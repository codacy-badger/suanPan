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

#include "CAX3.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>

const unsigned CAX3::m_node = 3;
const unsigned CAX3::m_dof = 2;
const unsigned CAX3::m_size = m_dof * m_node;

CAX3::CAX3(const unsigned T, uvec&& NT, const unsigned MT, const bool R)
	: MaterialElement(T, m_node, m_dof, std::forward<uvec>(NT), uvec{MT}, R) {}

void CAX3::initialize(const shared_ptr<DomainBase>& D) {
	const auto material_proto = std::dynamic_pointer_cast<Material2D>(D->get_material(unsigned(material_tag(0))));

	if(material_proto == nullptr || material_proto->plane_type != PlaneType::A) {
		D->disable_element(get_tag());
		return;
	}

	m_material = material_proto->get_copy();

	mat ele_coor(m_node, m_node);
	ele_coor.col(0).fill(1.);
	ele_coor.cols(1, 2) = get_coordinate(2);

	const rowvec ele_centre = mean(ele_coor);
	const auto& r = ele_centre(1);
	const auto& z = ele_centre(2);

	inv_coor = inv(ele_coor);

	access::rw(area) = .5 * det(ele_coor);
	access::rw(weight) = 2. * datum::pi * area * r;

	strain_mat.zeros(4, m_size);
	for(unsigned J = 0, K = 0; J < m_node; ++J, K += m_dof) {
		strain_mat(0, K) = strain_mat(3, K + 1) = inv_coor(1, J);
		strain_mat(3, K) = strain_mat(1, K + 1) = inv_coor(2, J);
		strain_mat(2, K) = inv_coor(0, J) / r + inv_coor(1, J) + inv_coor(2, J) * z / r;
	}
	trial_stiffness = current_stiffness = initial_stiffness = weight * strain_mat.t() * m_material->get_initial_stiffness() * strain_mat;

	trial_mass = current_mass = initial_mass.zeros(m_size, m_size);
}

int CAX3::update_status() {
	vec t_strain(4, fill::zeros);
	for(unsigned I = 0, J = 0; I < m_node; ++I, J += m_dof) {
		const auto& t_disp = node_ptr[I].lock()->get_trial_displacement();
		t_strain(0) += t_disp(0) * inv_coor(1, I);
		t_strain(1) += t_disp(1) * inv_coor(2, I);
		t_strain(2) += t_disp(0) * strain_mat(2, J);
		t_strain(3) += t_disp(0) * inv_coor(2, I) + t_disp(1) * inv_coor(1, I);
	}

	if(m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

	trial_stiffness = weight * strain_mat.t() * m_material->get_trial_stiffness() * strain_mat;
	trial_resistance = weight * strain_mat.t() * m_material->get_trial_stress();

	return SUANPAN_SUCCESS;
}

int CAX3::commit_status() { return m_material->commit_status(); }

int CAX3::clear_status() { return m_material->clear_status(); }

int CAX3::reset_status() { return m_material->reset_status(); }

void CAX3::print() { suanpan_info("CAX3 element.\n"); }

#ifdef SUANPAN_VTK
#include <vtkTriangle.h>

void CAX3::Setup() {
	vtk_cell = vtkSmartPointer<vtkTriangle>::New();
	auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < m_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void CAX3::GetDisplacement(vtkSmartPointer<vtkDoubleArray>& arrays) {
	mat t_disp(6, m_node, fill::zeros);
	t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);
	for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

void CAX3::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
	for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
