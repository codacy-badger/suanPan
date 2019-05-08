////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2019 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "F31.h"
#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>
#include <Section/Section.h>
#include <Toolbox/IntegrationPlan.h>

const unsigned F31::b_node = 2;
const unsigned F31::b_dof = 6;
const unsigned F31::b_size = b_node * b_dof;
const span F31::b_span(0, 2);
const unsigned F31::max_iteration = 20;
const double F31::tolerance = 1E-13;

F31::IntegrationPoint::IntegrationPoint(const double C, const double W, unique_ptr<Section>&& M)
	: coor(C)
	, weight(W)
	, b_section(std::forward<unique_ptr<Section>>(M))
	, strain_mat(3, 6, fill::zeros) {}

F31::F31(const unsigned T, uvec&& N, const unsigned S, const unsigned O, const unsigned P, const bool F)
	: SectionElement(T, b_node, b_dof, std::forward<uvec>(N), uvec{S}, F)
	, int_pt_num(P > 20 ? 20 : P)
	, orientation_tag(O) {}

void F31::initialize(const shared_ptr<DomainBase>& D) {
	auto& sec_proto = D->get_section(unsigned(section_tag(0)));

	if(sec_proto->get_section_type() != SectionType::D3) {
		suanpan_warning("initialize() needs a 3D section.\n");
		D->disable_element(get_tag());
		return;
	}

	if(!D->find_orientation(orientation_tag)) {
		suanpan_warning("initialize() needs a valid transformation.\n");
		D->disable_element(get_tag());
		return;
	}

	b_trans = D->get_orientation(orientation_tag)->get_copy();
	b_trans->set_element_ptr(this);

	access::rw(length) = b_trans->get_length();

	const mat sec_stiff = sec_proto->get_initial_stiffness()(b_span, b_span);

	const IntegrationPlan plan(1, int_pt_num, IntegrationType::LOBATTO);

	initial_local_flexibility.zeros(6, 6);
	int_pt.clear(), int_pt.reserve(int_pt_num);
	for(unsigned I = 0; I < int_pt_num; ++I) {
		int_pt.emplace_back(plan(I, 0), .5 * plan(I, 1), sec_proto->get_copy());
		int_pt[I].strain_mat(0, 0) = 1.;
		int_pt[I].strain_mat(1, 1) = int_pt[I].strain_mat(2, 3) = .5 * plan(I, 0) - .5;
		int_pt[I].strain_mat(1, 2) = int_pt[I].strain_mat(2, 4) = .5 * plan(I, 0) + .5;
		// factor .5 moved to weight
		initial_local_flexibility += int_pt[I].strain_mat.t() * solve(sec_stiff, int_pt[I].strain_mat * int_pt[I].weight * length);
	}
	access::rw(torsion_stiff) = 1E-2 * vec(initial_local_flexibility.diag()).head(5).min();
	initial_local_flexibility(5, 5) = torsion_stiff;
	trial_local_flexibility = current_local_flexibility = initial_local_flexibility;

	trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(inv(initial_local_flexibility));

	// trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(section_proto->get_parameter(ParameterType::LINEARDENSITY));

	trial_local_deformation = current_local_deformation.zeros(6);
	trial_local_resistance = current_local_resistance.zeros(6);
}

int F31::update_status() {
	b_trans->update_status();

	vec residual_deformation = -trial_local_deformation;

	// transform global deformation to local one (remove rigid body motion)
	trial_local_deformation = b_trans->to_local_vec(get_trial_displacement());

	// initial residual be aware of how to compute it
	residual_deformation += trial_local_deformation;

	auto counter = 0;
	while(true) {
		trial_local_resistance += solve(trial_local_flexibility, residual_deformation);
		residual_deformation.zeros();
		trial_local_flexibility.zeros();
		for(const auto& I : int_pt) {
			const vec target_section_resistance = I.strain_mat * trial_local_resistance;
			// compute unbalanced deformation
			const vec incre_deformation = solve(I.b_section->get_trial_stiffness()(b_span, b_span), target_section_resistance - I.b_section->get_trial_resistance()(b_span));
			// update status
			if(I.b_section->update_trial_status(I.b_section->get_trial_deformation() + incre_deformation) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
			// collect new flexibility and deformation
			trial_local_flexibility += I.weight * length * I.strain_mat.t() * solve(I.b_section->get_trial_stiffness()(b_span, b_span), I.strain_mat);
			residual_deformation += I.weight * length * I.strain_mat.t() * solve(I.b_section->get_trial_stiffness()(b_span, b_span), I.b_section->get_trial_resistance()(b_span) - target_section_resistance);
		}
		trial_local_flexibility(5, 5) = torsion_stiff;
		// quit if converged
		if(norm(residual_deformation) < tolerance) break;
		// impose a relatively more strict rule
		if(++counter == max_iteration) {
			suanpan_extra_debug("iteration fails to converge at element level.\n");
			return SUANPAN_FAIL;
		}
	}

	trial_stiffness = b_trans->to_global_stiffness_mat(inv(trial_local_flexibility));
	trial_resistance = b_trans->to_global_vec(trial_local_resistance);

	if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(trial_local_resistance);

	return SUANPAN_SUCCESS;
}

int F31::clear_status() {
	trial_local_flexibility = current_local_flexibility = initial_local_flexibility;
	trial_local_deformation = current_local_deformation.zeros();
	trial_local_resistance = current_local_resistance.zeros();
	auto code = 0;
	for(const auto& I : int_pt) code += I.b_section->clear_status();
	return code;
}

int F31::commit_status() {
	current_local_flexibility = trial_local_flexibility;
	current_local_deformation = trial_local_deformation;
	current_local_resistance = trial_local_resistance;
	auto code = 0;
	for(const auto& I : int_pt) code += I.b_section->commit_status();
	return code;
}

int F31::reset_status() {
	trial_local_flexibility = current_local_flexibility;
	trial_local_deformation = current_local_deformation;
	trial_local_resistance = current_local_resistance;
	auto code = 0;
	for(const auto& I : int_pt) code += I.b_section->reset_status();
	return code;
}

vector<vec> F31::record(const OutputType P) {
	vector<vec> output;
	output.reserve(int_pt.size());

	if(P == OutputType::E) for(const auto& I : int_pt) output.emplace_back(I.b_section->get_current_deformation());
	else if(P == OutputType::S) for(const auto& I : int_pt) output.emplace_back(I.b_section->get_current_resistance());
	else if(P == OutputType::PE) for(const auto& I : int_pt) output.emplace_back(I.b_section->get_current_deformation() - solve(I.b_section->get_initial_stiffness(), I.b_section->get_current_resistance()));

	return output;
}

void F31::print() {
	suanpan_info("A 2D forced based beam element%s.\n", nlgeom ? " and corotational formulation" : "");
	suanpan_info("The element connects nodes:");
	node_encoding.t().print();
	suanpan_info("Section Model:\n");
	auto J = 1;
	for(const auto& I : int_pt) {
		suanpan_info("IP %d: ", J++);
		I.b_section->print();
	}
}
