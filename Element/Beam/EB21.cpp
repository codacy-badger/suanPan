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

#include "EB21.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Element/Utility/B2DC.h>

const unsigned EB21::b_node = 2;
const unsigned EB21::b_dof = 3;
const unsigned EB21::b_size = b_dof * b_node;

EB21::EB21(const unsigned T, uvec&& N, const double A, const double I, const unsigned M, const bool F)
	: MaterialElement(T, b_node, b_dof, std::forward<uvec>(N), uvec{M}, F)
	, area(A)
	, moment_inertia(I)
	, b_trans(F ? make_unique<B2DC>() : make_unique<B2DL>()) {}

void EB21::initialize(const shared_ptr<DomainBase>& D) {
	b_trans->set_element_ptr(this);

	b_material = D->get_material(unsigned(material_tag(0)))->get_copy();

	// stiffness
	const auto tmp_a = as_scalar(b_material->get_initial_stiffness()) / b_trans->get_length();

	local_stiff.zeros(3, 3);
	local_stiff(0, 0) = tmp_a * area;
	local_stiff(1, 1) = local_stiff(2, 2) = 2. * (local_stiff(1, 2) = local_stiff(2, 1) = 2. * tmp_a * moment_inertia);

	trial_stiffness = current_stiffness = initial_stiffness = b_trans->to_global_stiffness_mat(local_stiff);

	trial_mass = current_mass = initial_mass = b_trans->to_global_mass_mat(b_material->get_parameter(ParameterType::DENSITY) * area);
}

int EB21::update_status() {
	b_trans->update_status();

	const vec local_force = local_stiff * b_trans->to_local_vec(get_trial_displacement());

	trial_stiffness = b_trans->to_global_stiffness_mat(local_stiff);
	trial_resistance = b_trans->to_global_vec(local_force);

	if(nlgeom) trial_geometry = b_trans->to_global_geometry_mat(local_force);

	return SUANPAN_SUCCESS;
}

int EB21::commit_status() { return b_material->commit_status(); }

int EB21::clear_status() { return b_material->clear_status(); }

int EB21::reset_status() { return b_material->reset_status(); }

void EB21::print() { suanpan_info("An elastic B21 element%s", nlgeom ? " with corotational formulation.\n" : ".\n"); }

#ifdef SUANPAN_VTK
#include <vtkLine.h>

void EB21::Setup() {
	vtk_cell = vtkSmartPointer<vtkLine>::New();
	auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < b_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void EB21::GetDisplacement(vtkSmartPointer<vtkDoubleArray>& arrays) {
	mat t_disp(6, b_node, fill::zeros);
	t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_displacement(), b_dof, b_node);
	for(unsigned I = 0; I < b_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

void EB21::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * mat(reshape(get_current_displacement(), b_dof, b_node).t()).cols(0, 1);
	for(unsigned I = 0; I < b_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
