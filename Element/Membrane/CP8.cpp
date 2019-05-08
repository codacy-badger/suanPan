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

#include "CP8.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/utility.h>
#include <Toolbox/tensorToolbox.h>

const unsigned CP8::m_node = 8;
const unsigned CP8::m_dof = 2;
const unsigned CP8::m_size = m_dof * m_node;

CP8::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& PNPXY)
	: coor(std::forward<vec>(C))
	, weight(W)
	, m_material(std::forward<unique_ptr<Material>>(M))
	, pn_pxy(std::forward<mat>(PNPXY))
	, strain_mat(3, m_size, fill::zeros) {}

CP8::CP8(const unsigned T, uvec&& N, const unsigned M, const double TH, const bool R, const bool F)
	: MaterialElement(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, F)
	, thickness(TH)
	, reduced_scheme(R) {}

void CP8::initialize(const shared_ptr<DomainBase>& D) {
	const auto material_proto = std::dynamic_pointer_cast<Material2D>(D->get_material(unsigned(material_tag(0))));

	if(material_proto == nullptr) {
		D->disable_element(get_tag());
		return;
	}

	if(material_proto->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

	const auto ele_coor = get_coordinate(m_dof);

	auto& ini_stiffness = material_proto->get_initial_stiffness();

	const IntegrationPlan plan(2, 2, reduced_scheme ? IntegrationType::GAUSS : IntegrationType::IRONS);

	initial_stiffness.zeros(m_size, m_size);

	int_pt.clear(), int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		vec t_vec{plan(I, 0), plan(I, 1)};
		const auto pn = shape::quad(t_vec, 1, m_node);
		const mat jacob = pn * ele_coor;
		int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob), material_proto->get_copy(), solve(jacob, pn));

		auto& c_int_pt = int_pt.back();

		for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
			c_int_pt.strain_mat(0, K) = c_int_pt.strain_mat(2, L) = c_int_pt.pn_pxy(0, J);
			c_int_pt.strain_mat(2, K) = c_int_pt.strain_mat(1, L) = c_int_pt.pn_pxy(1, J);
		}
		initial_stiffness += c_int_pt.weight * thickness * c_int_pt.strain_mat.t() * ini_stiffness * c_int_pt.strain_mat;
	}
	trial_stiffness = current_stiffness = initial_stiffness;

	initial_mass.zeros(m_size, m_size);
	const auto tmp_density = material_proto->get_parameter();
	if(tmp_density != 0.) {
		for(const auto& I : int_pt) {
			const auto n_int = shape::quad(I.coor, 0, m_node);
			const auto tmp_a = tmp_density * I.weight * thickness;
			for(unsigned J = 0; J < m_node; ++J) for(auto K = J; K < m_node; ++K) initial_mass(m_dof * J, m_dof * K) += tmp_a * n_int(J) * n_int(K);
		}

		for(unsigned I = 0; I < m_size; I += m_dof) {
			initial_mass(I + 1, I + 1) = initial_mass(I, I);
			for(auto J = I + m_dof; J < m_size; J += m_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
		}
	}
	trial_mass = current_mass = initial_mass;
}

int CP8::update_status() {
	trial_stiffness.zeros(m_size, m_size);
	trial_resistance.zeros(m_size);

	if(nlgeom) {
		trial_geometry.zeros(m_size, m_size);

		const mat ele_disp = reshape(get_trial_displacement(), m_dof, m_node);

		mat BN(3, m_size);
		for(const auto& I : int_pt) {
			const mat gradient = ele_disp * I.pn_pxy.t() + eye(m_dof, m_dof);
			for(unsigned J = 0; J < m_node; ++J) {
				const auto IDXA = m_dof * J;
				const auto IDXB = IDXA + 1;
				BN(0, IDXA) = I.pn_pxy(0, J) * gradient(0, 0);
				BN(1, IDXA) = I.pn_pxy(1, J) * gradient(0, 1);
				BN(0, IDXB) = I.pn_pxy(0, J) * gradient(1, 0);
				BN(1, IDXB) = I.pn_pxy(1, J) * gradient(1, 1);
				BN(2, IDXA) = I.pn_pxy(0, J) * gradient(0, 1) + I.pn_pxy(1, J) * gradient(0, 0);
				BN(2, IDXB) = I.pn_pxy(0, J) * gradient(1, 1) + I.pn_pxy(1, J) * gradient(1, 0);
			}

			if(I.m_material->update_trial_status(tensor::strain::to_voigt(tensor::strain::to_green(gradient))) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

			const auto t_weight = I.weight * thickness;

			auto& t_stress = I.m_material->get_trial_stress();

			const auto sigma = tensor::stress::to_tensor(t_stress);

			for(unsigned J = 0; J < m_node; ++J) {
				const vec t_vec = sigma * I.pn_pxy.col(J);
				auto t_factor = t_weight * dot(I.pn_pxy.col(J), t_vec);
				const auto IDXA = m_dof * J;
				trial_geometry(IDXA, IDXA) += t_factor;
				trial_geometry(IDXA + 1, IDXA + 1) += t_factor;
				for(auto K = J + 1; K < m_node; ++K) {
					t_factor = t_weight * dot(I.pn_pxy.col(K), t_vec);
					const auto IDXB = m_dof * K;
					trial_geometry(IDXB, IDXA) += t_factor;
					trial_geometry(IDXB + 1, IDXA + 1) += t_factor;
					trial_geometry(IDXA, IDXB) += t_factor;
					trial_geometry(IDXA + 1, IDXB + 1) += t_factor;
				}
			}

			trial_stiffness += t_weight * BN.t() * I.m_material->get_trial_stiffness() * BN;
			trial_resistance += t_weight * BN.t() * t_stress;
		}
	} else
		for(const auto& I : int_pt) {
			vec t_strain(3, fill::zeros);
			for(unsigned J = 0; J < m_node; ++J) {
				const auto& t_disp = node_ptr[J].lock()->get_trial_displacement();
				t_strain(0) += t_disp(0) * I.pn_pxy(0, J);
				t_strain(1) += t_disp(1) * I.pn_pxy(1, J);
				t_strain(2) += t_disp(0) * I.pn_pxy(1, J) + t_disp(1) * I.pn_pxy(0, J);
			}
			if(I.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

			trial_stiffness += I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat * I.weight * thickness;
			trial_resistance += I.strain_mat.t() * I.m_material->get_trial_stress() * I.weight * thickness;
		}

	return SUANPAN_SUCCESS;
}

int CP8::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	return code;
}

int CP8::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	return code;
}

int CP8::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	return code;
}

void CP8::print() {
	suanpan_info("A CP8%s element%s.\n", reduced_scheme ? "R" : "", nlgeom ? " with nonlinear geometry on" : "");
	suanpan_info("The nodes connected are: ");
	node_encoding.t().print();
	suanpan_info("Material models:\n");
	for(const auto& I : int_pt) I.m_material->print();
}

#ifdef SUANPAN_VTK
#include <vtkQuadraticQuad.h>

void CP8::Setup() {
	vtk_cell = vtkSmartPointer<vtkQuadraticQuad>::New();
	auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < m_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void CP8::GetDisplacement(vtkSmartPointer<vtkDoubleArray>& arrays) {
	mat t_disp(6, m_node, fill::zeros);
	t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);
	for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

mat CP8::GetData(const OutputType P) {
	if(P == OutputType::S) {
		mat B(3 * int_pt.size(), 1);
		mat A(3 * int_pt.size(), 11);
		for(unsigned I = 0, J = 0; I < unsigned(int_pt.size()); ++I, J = J + 3) {
			B.rows(J, J + 2) = int_pt[I].m_material->get_current_stress();
			A.rows(J, J + 2) = shape::stress11(int_pt[I].coor);
		}
		mat X;
		if(!solve(X, A, B)) return {};

		mat t_stress(6, m_node, fill::zeros);

		t_stress(uvec{0, 1, 3}, uvec{0}) = shape::stress11(-1., -1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{1}) = shape::stress11(1., -1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{2}) = shape::stress11(1., 1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{3}) = shape::stress11(-1., 1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{4}) = shape::stress11(0., -1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{5}) = shape::stress11(1., 0.) * X;
		t_stress(uvec{0, 1, 3}, uvec{6}) = shape::stress11(0., 1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{7}) = shape::stress11(-1., 0.) * X;

		return t_stress;
	}

	return {};
}

void CP8::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
	for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
