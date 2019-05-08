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

#include "C3D4.h"
#include <Domain/DomainBase.h>
#include <Material/Material3D/Material3D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/tensorToolbox.h>

const unsigned C3D4::c_node = 4;
const unsigned C3D4::c_dof = 3;
const unsigned C3D4::c_size = c_dof * c_node;

C3D4::C3D4(const unsigned T, uvec&& N, const unsigned M, const bool F)
	: MaterialElement(T, c_node, c_dof, std::forward<uvec>(N), uvec{M}, F) {}

void C3D4::initialize(const shared_ptr<DomainBase>& D) {
	const auto material_proto = std::dynamic_pointer_cast<Material3D>(D->get_material(unsigned(material_tag(0))));

	if(material_proto == nullptr) {
		D->disable_element(get_tag());
		return;
	}

	c_material = material_proto->get_copy();

	mat ele_coor(c_node, c_node);
	ele_coor.col(0).fill(1.);
	ele_coor.cols(1, 3) = get_coordinate(3);

	access::rw(volume) = det(ele_coor) / 6.;

	const mat inv_coor = inv(ele_coor);

	pn_pxyz = inv_coor.rows(1, 3);

	strain_mat.zeros(6, c_size);
	for(unsigned J = 0, K = 0, L = 1, M = 2; J < c_node; ++J, K += c_dof, L += c_dof, M += c_dof) {
		strain_mat(0, K) = strain_mat(3, L) = strain_mat(5, M) = pn_pxyz(0, J);
		strain_mat(3, K) = strain_mat(1, L) = strain_mat(4, M) = pn_pxyz(1, J);
		strain_mat(5, K) = strain_mat(4, L) = strain_mat(2, M) = pn_pxyz(2, J);
	}

	trial_stiffness = current_stiffness = initial_stiffness = volume * strain_mat.t() * c_material->get_initial_stiffness() * strain_mat;

	initial_mass.zeros(c_size, c_size);
	const auto t_density = c_material->get_parameter() * volume;
	if(t_density != 0.) {
		const rowvec n = mean(ele_coor) * inv_coor;
		for(unsigned I = 0; I < c_node; ++I) for(auto J = I; J < c_node; ++J) initial_mass(c_dof * I, c_dof * J) += t_density * n(I) * n(J);
		for(unsigned I = 0; I < c_size; I += c_dof) {
			initial_mass(I + 1, I + 1) = initial_mass(I, I);
			for(auto J = I + c_dof; J < c_size; J += c_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
		}
	}
	trial_mass = current_mass = initial_mass;
}

int C3D4::update_status() {
	if(nlgeom) {
		const mat ele_disp = reshape(get_trial_displacement(), c_dof, c_node);

		mat BN(6, c_size);

		const mat gradient = ele_disp * pn_pxyz.t() + eye(c_dof, c_dof);
		for(unsigned I = 0, J = 0, K = 1, L = 2; I < c_node; ++I, J += c_dof, K += c_dof, L += c_dof) {
			BN(0, J) = pn_pxyz(0, I) * gradient(0, 0);
			BN(1, J) = pn_pxyz(1, I) * gradient(0, 1);
			BN(2, J) = pn_pxyz(2, I) * gradient(0, 2);
			BN(0, K) = pn_pxyz(0, I) * gradient(1, 0);
			BN(1, K) = pn_pxyz(1, I) * gradient(1, 1);
			BN(2, K) = pn_pxyz(2, I) * gradient(1, 2);
			BN(0, L) = pn_pxyz(0, I) * gradient(2, 0);
			BN(1, L) = pn_pxyz(1, I) * gradient(2, 1);
			BN(2, L) = pn_pxyz(2, I) * gradient(2, 2);
			BN(3, J) = pn_pxyz(0, I) * gradient(0, 1) + pn_pxyz(1, I) * gradient(0, 0);
			BN(4, J) = pn_pxyz(1, I) * gradient(0, 2) + pn_pxyz(2, I) * gradient(0, 1);
			BN(5, J) = pn_pxyz(2, I) * gradient(0, 0) + pn_pxyz(0, I) * gradient(0, 2);
			BN(3, K) = pn_pxyz(0, I) * gradient(1, 1) + pn_pxyz(1, I) * gradient(1, 0);
			BN(4, K) = pn_pxyz(1, I) * gradient(1, 2) + pn_pxyz(2, I) * gradient(1, 1);
			BN(5, K) = pn_pxyz(2, I) * gradient(1, 0) + pn_pxyz(0, I) * gradient(1, 2);
			BN(3, L) = pn_pxyz(0, I) * gradient(2, 1) + pn_pxyz(1, I) * gradient(2, 0);
			BN(4, L) = pn_pxyz(1, I) * gradient(2, 2) + pn_pxyz(2, I) * gradient(2, 1);
			BN(5, L) = pn_pxyz(2, I) * gradient(2, 0) + pn_pxyz(0, I) * gradient(2, 2);
		}

		if(c_material->update_trial_status(tensor::strain::to_voigt(tensor::strain::to_green(gradient))) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		auto& t_stress = c_material->get_trial_stress();

		trial_stiffness = volume * BN.t() * c_material->get_trial_stiffness() * BN;
		trial_resistance = volume * BN.t() * t_stress;

		trial_geometry.zeros(c_size, c_size);
		const auto sigma = tensor::stress::to_tensor(t_stress);

		for(unsigned I = 0; I < c_node; ++I) {
			const vec t_vec = sigma * pn_pxyz.col(I);
			auto t_factor = volume * dot(pn_pxyz.col(I), t_vec);
			const auto K = c_dof * I;
			trial_geometry(K, K) += t_factor;
			trial_geometry(K + 1, K + 1) += t_factor;
			trial_geometry(K + 2, K + 2) += t_factor;
			for(auto J = I + 1; J < c_node; ++J) {
				t_factor = volume * dot(pn_pxyz.col(J), t_vec);
				const auto L = c_dof * J;
				trial_geometry(L, K) += t_factor;
				trial_geometry(L + 1, K + 1) += t_factor;
				trial_geometry(L + 2, K + 2) += t_factor;
				trial_geometry(K, L) += t_factor;
				trial_geometry(K + 1, L + 1) += t_factor;
				trial_geometry(K + 2, L + 2) += t_factor;
			}
		}
	} else {
		if(c_material->update_trial_status(strain_mat * get_trial_displacement()) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		trial_stiffness = volume * strain_mat.t() * c_material->get_trial_stiffness() * strain_mat;
		trial_resistance = volume * strain_mat.t() * c_material->get_trial_stress();
	}

	return SUANPAN_SUCCESS;
}

int C3D4::commit_status() { return c_material->commit_status(); }

int C3D4::clear_status() { return c_material->clear_status(); }

int C3D4::reset_status() { return c_material->reset_status(); }

vector<vec> C3D4::record(const OutputType T) { return c_material->record(T); }

void C3D4::print() {
	suanpan_info("C3D4 element connects:");
	node_encoding.t().print();
	suanpan_info("Material model:\n");
	c_material->print();
	suanpan_info("strain:\t");
	c_material->get_trial_strain().t().print();
	suanpan_info("stress:\t");
	c_material->get_trial_stress().t().print();
}

#ifdef SUANPAN_VTK
#include <vtkTetra.h>

void C3D4::Setup() {
	vtk_cell = vtkSmartPointer<vtkTetra>::New();
	auto ele_coor = get_coordinate(3);
	for(unsigned I = 0; I < c_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
	}
}

void C3D4::GetDisplacement(vtkSmartPointer<vtkDoubleArray>& arrays) {
	mat t_disp(6, c_node, fill::zeros);
	t_disp.rows(0, 2) = reshape(get_current_displacement(), c_dof, c_node);
	for(unsigned I = 0; I < c_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

void C3D4::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(3) + amplifier * reshape(get_current_displacement(), c_dof, c_node).t();
	for(unsigned I = 0; I < c_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
