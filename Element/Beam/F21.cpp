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

#include "F21.h"
#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>
#include <Section/Section.h>
#include <Toolbox/IntegrationPlan.h>

const unsigned F21::b_node = 2;
const unsigned F21::b_dof = 3;
const unsigned F21::b_size = b_node * b_dof;
const unsigned F21::max_iteration = 20;
const double F21::tolerance = 1E-14;

F21::IntegrationPoint::IntegrationPoint(const double C, const double W, unique_ptr<Section>&& M)
	: coor(C)
	, weight(W)
	, b_section(std::forward<unique_ptr<Section>>(M))
	, B(2, 3, fill::zeros) {}

F21::F21(const unsigned T, uvec&& N, const unsigned S, const unsigned P, const bool F)
	: SectionElement(T, b_node, b_dof, std::forward<uvec>(N), uvec{S}, F)
	, int_pt_num(P > 20 ? 20 : P)
	, b_trans(F ? make_unique<B2DC>() : make_unique<B2DL>()) {}

void F21::initialize(const shared_ptr<DomainBase>& D) {
	const auto& section_proto = D->get_section(unsigned(section_tag(0)));

	if(section_proto->get_section_type() != SectionType::D2) {
		suanpan_warning("initialize() needs a 2D section.\n");
		D->disable_element(get_tag());
		return;
	}

	b_trans->set_element_ptr(this);

	access::rw(length) = b_trans->get_length();

	const auto t_flexibility = inv(section_proto->get_initial_stiffness());

	const IntegrationPlan plan(1, int_pt_num, IntegrationType::LOBATTO);

	initial_local_flexibility.zeros(3, 3);
	int_pt.clear(), int_pt.reserve(int_pt_num);
	for(unsigned I = 0; I < int_pt_num; ++I) {
		int_pt.emplace_back(plan(I, 0), .5 * plan(I, 1), section_proto->get_copy());
		int_pt[I].B(0, 0) = 1.;
		int_pt[I].B(1, 1) = .5 * (plan(I, 0) - 1.);
		int_pt[I].B(1, 2) = .5 * (plan(I, 0) + 1.);
		// factor .5 moved to weight
		initial_local_flexibility += int_pt[I].B.t() * t_flexibility * int_pt[I].B * int_pt[I].weight * length;
	}
	trial_local_flexibility = current_local_flexibility = initial_local_flexibility;

	trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(inv(initial_local_flexibility));

	trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(section_proto->get_parameter(ParameterType::LINEARDENSITY));

	trial_local_deformation = current_local_deformation.zeros(3);
	trial_local_resistance = current_local_resistance.zeros(3);
}

int F21::update_status() {
	b_trans->update_status();

	auto residual_deformation = trial_local_deformation;

	// transform global deformation to local one (remove rigid body motion)
	trial_local_deformation = b_trans->to_local_vec(get_trial_displacement());

	// initial residual be aware of how to compute it
	residual_deformation -= trial_local_deformation;

	auto ref_error = 1.;

	unsigned counter = 0;
	while(++counter < max_iteration) {
		const vec incre_resistance = solve(trial_local_flexibility, residual_deformation, solve_opts::equilibrate);

		// quit if converged
		auto abs_error = norm(incre_resistance);
		if(1 == counter) ref_error = std::max(1., abs_error);
		suanpan_extra_debug("F21 local iteration error: %.4E.\n", abs_error /= ref_error);
		if(norm(residual_deformation) <= tolerance && abs_error <= tolerance) break;

		trial_local_resistance -= incre_resistance;
		residual_deformation.zeros();
		trial_local_flexibility.zeros();
		for(const auto& I : int_pt) {
			const vec target_section_resistance = I.B * trial_local_resistance;
			// compute unbalanced deformation
			const vec incre_deformation = solve(I.b_section->get_trial_stiffness(), target_section_resistance - I.b_section->get_trial_resistance());
			// update status
			if(I.b_section->update_trial_status(I.b_section->get_trial_deformation() + incre_deformation) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
			// collect new flexibility and deformation
			trial_local_flexibility += I.weight * length * I.B.t() * solve(I.b_section->get_trial_stiffness(), I.B);
			residual_deformation += I.weight * length * I.B.t() * solve(I.b_section->get_trial_stiffness(), target_section_resistance - I.b_section->get_trial_resistance());
		}
	}

	if(max_iteration == counter) {
		suanpan_extra_debug("iteration fails to converge at element level.\n");
		return SUANPAN_FAIL;
	}

	trial_stiffness = b_trans->to_global_stiffness_mat(inv(trial_local_flexibility));
	trial_resistance = b_trans->to_global_vec(trial_local_resistance);

	if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(trial_local_resistance);

	return SUANPAN_SUCCESS;
}

int F21::clear_status() {
	trial_local_flexibility = current_local_flexibility = initial_local_flexibility;
	current_local_deformation.zeros();
	trial_local_deformation.zeros();
	current_local_resistance.zeros();
	trial_local_resistance.zeros();
	auto code = 0;
	for(const auto& I : int_pt) code += I.b_section->clear_status();
	return code;
}

int F21::commit_status() {
	current_local_flexibility = trial_local_flexibility;
	current_local_deformation = trial_local_deformation;
	current_local_resistance = trial_local_resistance;
	auto code = 0;
	for(const auto& I : int_pt) code += I.b_section->commit_status();
	return code;
}

int F21::reset_status() {
	trial_local_flexibility = current_local_flexibility;
	trial_local_deformation = current_local_deformation;
	trial_local_resistance = current_local_resistance;
	auto code = 0;
	for(const auto& I : int_pt) code += I.b_section->reset_status();
	return code;
}

vector<vec> F21::record(const OutputType P) {
	vector<vec> output;
	output.reserve(int_pt.size());

	if(P == OutputType::E) for(const auto& I : int_pt) output.emplace_back(I.b_section->get_current_deformation());
	else if(P == OutputType::S) for(const auto& I : int_pt) output.emplace_back(I.b_section->get_current_resistance());
	else if(P == OutputType::PE) for(const auto& I : int_pt) output.emplace_back(I.b_section->get_current_deformation() - solve(I.b_section->get_initial_stiffness(), I.b_section->get_current_resistance()));

	return output;
}

void F21::print() {
	suanpan_info("A 2D forced based beam element%sReference: https://doi.org/10.1016/0045-7949(95)00103-N \n", nlgeom ? " and corotational formulation.\n" : ".\n");
	suanpan_info("The element connects nodes:");
	node_encoding.t().print();
	suanpan_info("Section Model:\n");
	auto J = 1;
	for(const auto& I : int_pt) {
		suanpan_info("IP %d: ", J++);
		I.b_section->print();
	}
}
