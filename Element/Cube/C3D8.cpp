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

#include "C3D8.h"
#include <Domain/DomainBase.h>
#include <Material/Material3D/Material3D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/tensorToolbox.h>

const unsigned C3D8::c_node = 8;
const unsigned C3D8::c_dof = 3;
const unsigned C3D8::c_size = c_dof * c_node;

C3D8::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& P)
	: coor(std::forward<vec>(C))
	, weight(W)
	, c_material(std::forward<unique_ptr<Material>>(M))
	, pn_pxyz(std::forward<mat>(P))
	, strain_mat(6, c_size, fill::zeros) {}

C3D8::C3D8(const unsigned T, uvec&& N, const unsigned M, const char R, const bool F)
	: MaterialElement(T, c_node, c_dof, std::forward<uvec>(N), uvec{M}, F)
	, int_scheme(R) {}

void C3D8::initialize(const shared_ptr<DomainBase>& D) {
	const auto ele_coor = get_coordinate(c_dof);

	auto& material_proto = D->get_material(unsigned(material_tag(0)));

	auto& ini_stiffness = material_proto->get_initial_stiffness();

	const IntegrationPlan plan(3, int_scheme == 'R' ? 1 : 2, int_scheme == 'I' ? IntegrationType::IRONS : IntegrationType::GAUSS);

	initial_stiffness.zeros(c_size, c_size);

	int_pt.clear(), int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		vec t_vec{plan(I, 0), plan(I, 1), plan(I, 2)};
		const auto pn = shape::cube(t_vec, 1);
		const mat jacob = pn * ele_coor;
		int_pt.emplace_back(std::move(t_vec), plan(I, c_dof) * det(jacob), material_proto->get_copy(), solve(jacob, pn));

		auto& c_pt = int_pt.back();
		for(unsigned J = 0; J < c_node; ++J) {
			const auto K = c_dof * J;
			c_pt.strain_mat(0, K) = c_pt.strain_mat(3, K + 1) = c_pt.strain_mat(5, K + 2) = c_pt.pn_pxyz(0, J);
			c_pt.strain_mat(3, K) = c_pt.strain_mat(1, K + 1) = c_pt.strain_mat(4, K + 2) = c_pt.pn_pxyz(1, J);
			c_pt.strain_mat(5, K) = c_pt.strain_mat(4, K + 1) = c_pt.strain_mat(2, K + 2) = c_pt.pn_pxyz(2, J);
		}
		initial_stiffness += c_pt.weight * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;
	}
	trial_stiffness = current_stiffness = initial_stiffness;

	initial_mass.zeros(c_size, c_size);
	const auto t_density = material_proto->get_parameter();
	if(t_density != 0.) {
		for(const auto& I : int_pt) {
			const auto n_int = shape::cube(I.coor, 0);
			const auto tmp_a = t_density * I.weight;
			for(unsigned J = 0; J < c_node; ++J) for(auto K = J; K < c_node; ++K) initial_mass(c_dof * J, c_dof * K) += tmp_a * n_int(J) * n_int(K);
		}
		for(unsigned I = 0, K = 1, L = 2; I < c_size; I += c_dof, K += c_dof, L += c_dof) {
			initial_mass(K, K) = initial_mass(L, L) = initial_mass(I, I);
			for(auto J = I + c_dof, M = J + 1, N = J + 2; J < c_size; J += c_dof, M += c_dof, N += c_dof) initial_mass(J, I) = initial_mass(K, M) = initial_mass(L, N) = initial_mass(M, K) = initial_mass(N, L) = initial_mass(I, J);
		}
	}
	trial_mass = current_mass = initial_mass;
}

int C3D8::update_status() {
	const auto t_disp = get_trial_displacement();

	trial_stiffness.zeros(c_size, c_size);
	trial_resistance.zeros(c_size);

	if(nlgeom) {
		const mat ele_disp = reshape(t_disp, c_dof, c_node);

		trial_geometry.zeros(c_size, c_size);

		mat BN(6, c_size);
		for(const auto& I : int_pt) {
			const mat gradient = ele_disp * I.pn_pxyz.t() + eye(c_dof, c_dof);
			for(unsigned J = 0, K = 0, L = 1, M = 2; J < c_node; ++J, K += c_dof, L += c_dof, M += c_dof) {
				BN(0, K) = I.pn_pxyz(0, J) * gradient(0, 0);
				BN(1, K) = I.pn_pxyz(1, J) * gradient(0, 1);
				BN(2, K) = I.pn_pxyz(2, J) * gradient(0, 2);
				BN(0, L) = I.pn_pxyz(0, J) * gradient(1, 0);
				BN(1, L) = I.pn_pxyz(1, J) * gradient(1, 1);
				BN(2, L) = I.pn_pxyz(2, J) * gradient(1, 2);
				BN(0, M) = I.pn_pxyz(0, J) * gradient(2, 0);
				BN(1, M) = I.pn_pxyz(1, J) * gradient(2, 1);
				BN(2, M) = I.pn_pxyz(2, J) * gradient(2, 2);
				BN(3, K) = I.pn_pxyz(0, J) * gradient(0, 1) + I.pn_pxyz(1, J) * gradient(0, 0);
				BN(4, K) = I.pn_pxyz(1, J) * gradient(0, 2) + I.pn_pxyz(2, J) * gradient(0, 1);
				BN(5, K) = I.pn_pxyz(2, J) * gradient(0, 0) + I.pn_pxyz(0, J) * gradient(0, 2);
				BN(3, L) = I.pn_pxyz(0, J) * gradient(1, 1) + I.pn_pxyz(1, J) * gradient(1, 0);
				BN(4, L) = I.pn_pxyz(1, J) * gradient(1, 2) + I.pn_pxyz(2, J) * gradient(1, 1);
				BN(5, L) = I.pn_pxyz(2, J) * gradient(1, 0) + I.pn_pxyz(0, J) * gradient(1, 2);
				BN(3, M) = I.pn_pxyz(0, J) * gradient(2, 1) + I.pn_pxyz(1, J) * gradient(2, 0);
				BN(4, M) = I.pn_pxyz(1, J) * gradient(2, 2) + I.pn_pxyz(2, J) * gradient(2, 1);
				BN(5, M) = I.pn_pxyz(2, J) * gradient(2, 0) + I.pn_pxyz(0, J) * gradient(2, 2);
			}

			if(I.c_material->update_trial_status(tensor::strain::to_voigt(tensor::strain::to_green(gradient))) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

			auto& t_stress = I.c_material->get_trial_stress();

			const auto sigma = tensor::stress::to_tensor(t_stress);

			for(unsigned J = 0, K = 0, L = 1, M = 2; J < c_node; ++J, K += c_dof, L += c_dof, M += c_dof) {
				const vec t_vec = I.weight * sigma * I.pn_pxyz.col(J);
				auto t_factor = dot(I.pn_pxyz.col(J), t_vec);
				trial_geometry(K, K) += t_factor;
				trial_geometry(L, L) += t_factor;
				trial_geometry(M, M) += t_factor;
				for(auto N = J + 1, O = c_dof * N, P = O + 1, Q = P + 1; N < c_node; ++N, O += c_dof, P += c_dof, Q += c_dof) {
					t_factor = dot(I.pn_pxyz.col(N), t_vec);
					trial_geometry(O, K) += t_factor;
					trial_geometry(P, L) += t_factor;
					trial_geometry(Q, M) += t_factor;
					trial_geometry(K, O) += t_factor;
					trial_geometry(L, P) += t_factor;
					trial_geometry(M, Q) += t_factor;
				}
			}

			trial_stiffness += I.weight * BN.t() * I.c_material->get_trial_stiffness() * BN;
			trial_resistance += I.weight * BN.t() * t_stress;
		}
	} else
		for(const auto& I : int_pt) {
			if(I.c_material->update_trial_status(I.strain_mat * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
			trial_stiffness += I.weight * I.strain_mat.t() * I.c_material->get_trial_stiffness() * I.strain_mat;
			trial_resistance += I.weight * I.strain_mat.t() * I.c_material->get_trial_stress();
		}

	return SUANPAN_SUCCESS;
}

int C3D8::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.c_material->commit_status();
	return code;
}

int C3D8::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.c_material->clear_status();
	return code;
}

int C3D8::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.c_material->reset_status();
	return code;
}

vector<vec> C3D8::record(const OutputType T) {
	vector<vec> data;
	switch(T) {
	case OutputType::E:
		for(const auto& I : int_pt) data.emplace_back(I.c_material->get_trial_strain());
		break;
	case OutputType::S:
		for(const auto& I : int_pt) data.emplace_back(I.c_material->get_trial_stress());
		break;
	default:
		for(const auto& I : int_pt) for(auto J : I.c_material->record(T)) data.emplace_back(J);
		break;
	}
	return data;
}

void C3D8::print() {
	suanpan_info("C3D8 element%s%s", int_scheme == 'R' ? " reduced integration" : int_scheme == 'I' ? " Iron's integration" : " full integration", nlgeom ? " nonlinear geometry" : "");
	suanpan_info(".\nThe element connects nodes:");
	node_encoding.t().print();
	suanpan_info("Material models:\n");
	for(const auto& t_pt : int_pt) {
		t_pt.c_material->print();
		suanpan_info("strain:\t");
		t_pt.c_material->get_trial_strain().t().print();
		suanpan_info("stress:\t");
		t_pt.c_material->get_trial_stress().t().print();
	}
}

#ifdef SUANPAN_VTK
#include <vtkHexahedron.h>

void C3D8::Setup() {
	vtk_cell = vtkSmartPointer<vtkHexahedron>::New();
	auto ele_coor = get_coordinate(3);
	for(unsigned I = 0; I < c_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
	}
}

void C3D8::GetDisplacement(vtkSmartPointer<vtkDoubleArray>& arrays) {
	mat t_disp(6, c_node, fill::zeros);
	t_disp.rows(0, 2) = reshape(get_current_displacement(), c_dof, c_node);
	for(unsigned I = 0; I < c_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

void C3D8::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(3) + amplifier * reshape(get_current_displacement(), c_dof, c_node).t();
	for(unsigned I = 0; I < c_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif
