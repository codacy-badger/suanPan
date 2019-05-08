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

#include "B21.h"
#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>
#include <Section/Section.h>
#include <Toolbox/IntegrationPlan.h>

const unsigned B21::b_node = 2;
const unsigned B21::b_dof = 3;
const unsigned B21::b_size = b_node * b_dof;

B21::IntegrationPoint::IntegrationPoint(const double C, const double W, unique_ptr<Section>&& M)
	: coor(C)
	, weight(W)
	, b_section(std::forward<unique_ptr<Section>>(M))
	, strain_mat(2, 3, fill::zeros) {}

B21::B21(const unsigned T, uvec&& N, const unsigned S, const unsigned P, const bool F)
	: SectionElement(T, b_node, b_dof, std::forward<uvec>(N), uvec{S}, F)
	, int_pt_num(P > 20 ? 20 : P)
	, b_trans(F ? make_unique<B2DC>() : make_unique<B2DL>()) {}

void B21::initialize(const shared_ptr<DomainBase>& D) {
	auto& section_proto = D->get_section(unsigned(section_tag(0)));

	if(section_proto->get_section_type() != SectionType::D2) {
		suanpan_warning("initialize() needs a 2D section.\n");
		D->disable_element(get_tag());
		return;
	}

	b_trans->set_element_ptr(this);

	access::rw(length) = b_trans->get_length();

	const IntegrationPlan plan(1, int_pt_num, IntegrationType::LOBATTO);

	mat local_stiffness(3, 3, fill::zeros);
	int_pt.clear(), int_pt.reserve(int_pt_num);
	for(unsigned I = 0; I < int_pt_num; ++I) {
		int_pt.emplace_back(plan(I, 0), .5 * plan(I, 1), section_proto->get_copy());
		int_pt[I].strain_mat(0, 0) = 1.;
		int_pt[I].strain_mat(1, 1) = 3. * plan(I, 0) - 1.;
		int_pt[I].strain_mat(1, 2) = 3. * plan(I, 0) + 1.;
		local_stiffness += int_pt[I].strain_mat.t() * int_pt[I].b_section->get_initial_stiffness() * int_pt[I].strain_mat * int_pt[I].weight / length;
	}

	trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);

	trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(section_proto->get_parameter(ParameterType::LINEARDENSITY));
}

int B21::update_status() {
	b_trans->update_status();

	// const auto length = b_trans->get_length();

	const auto local_deformation = b_trans->to_local_vec(get_trial_displacement());

	mat local_stiffness(3, 3, fill::zeros);
	vec local_resistance(3, fill::zeros);
	for(const auto& I : int_pt) {
		if(I.b_section->update_trial_status(I.strain_mat * local_deformation / length) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		local_stiffness += I.strain_mat.t() * I.b_section->get_trial_stiffness() * I.strain_mat * I.weight / length;
		local_resistance += I.strain_mat.t() * I.b_section->get_trial_resistance() * I.weight;
	}

	trial_stiffness = b_trans->to_global_stiffness_mat(local_stiffness);
	trial_resistance = b_trans->to_global_vec(local_resistance);

	if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(local_resistance);

	return SUANPAN_SUCCESS;
}

int B21::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.b_section->commit_status();
	return code;
}

int B21::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.b_section->clear_status();
	return code;
}

int B21::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.b_section->reset_status();
	return code;
}

vector<vec> B21::record(const OutputType P) {
	vector<vec> output;
	output.reserve(int_pt.size());

	if(P == OutputType::E) for(const auto& I : int_pt) output.emplace_back(I.b_section->get_current_deformation());
	else if(P == OutputType::S) for(const auto& I : int_pt) output.emplace_back(I.b_section->get_current_resistance());
	else if(P == OutputType::PE) for(const auto& I : int_pt) output.emplace_back(I.b_section->get_current_deformation() - solve(I.b_section->get_initial_stiffness(), I.b_section->get_current_resistance()));

	return output;
}

void B21::print() {
	suanpan_info("A classic 2D displacement based beam element using Hermite interpolation functions%s", nlgeom ? " and corotational formulation.\n" : ".\n");
	suanpan_info("The element connects nodes:");
	node_encoding.t().print();
	suanpan_info("Section Model:\n");
	auto J = 1;
	for(const auto& I : int_pt) {
		suanpan_info("IP %d: ", J++);
		I.b_section->print();
	}
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

void B21::Setup() {
	vtk_cell = vtkSmartPointer<vtkLine>::New();
	auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < b_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void B21::GetDisplacement(vtkSmartPointer<vtkDoubleArray>& arrays) {
	mat t_disp(6, b_node, fill::zeros);
	t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_displacement(), b_dof, b_node);
	for(unsigned I = 0; I < b_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

void B21::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * mat(reshape(get_current_displacement(), b_dof, b_node).t()).cols(0, 1);
	for(unsigned I = 0; I < b_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
