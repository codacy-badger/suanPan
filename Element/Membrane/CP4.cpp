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

#include "CP4.h"
#include <Domain/DomainBase.h>
#include <Element/Utility/MatrixModifier.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/utility.h>
#include <Toolbox/tensorToolbox.h>

const unsigned CP4::m_node = 4;
const unsigned CP4::m_dof = 2;
const unsigned CP4::m_size = m_dof * m_node;

CP4::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& PNPXY)
	: coor(std::forward<vec>(C))
	, weight(W)
	, m_material(std::forward<unique_ptr<Material>>(M))
	, pn_pxy(std::forward<mat>(PNPXY))
	, strain_mat(3, m_size, fill::zeros) {}

void CP4::stack_stiffness(mat& K, const mat& D, const mat& N, const double F) {
	const auto D11 = F * D(0, 0);
	const auto D12 = F * D(0, 1);
	const auto D13 = F * D(0, 2);
	const auto D21 = F * D(1, 0);
	const auto D22 = F * D(1, 1);
	const auto D23 = F * D(1, 2);
	const auto D31 = F * D(2, 0);
	const auto D32 = F * D(2, 1);
	const auto D33 = F * D(2, 2);

	const auto& NX1 = N(0, 0);
	const auto& NY1 = N(1, 1);
	const auto& NX2 = N(0, 2);
	const auto& NY2 = N(1, 3);
	const auto& NX3 = N(0, 4);
	const auto& NY3 = N(1, 5);
	const auto& NX4 = N(0, 6);
	const auto& NY4 = N(1, 7);

	const auto D11NX1 = D11 * NX1;
	const auto D11NX2 = D11 * NX2;
	const auto D11NX3 = D11 * NX3;
	const auto D11NX4 = D11 * NX4;

	const auto D12NX1 = D12 * NX1;
	const auto D12NX2 = D12 * NX2;
	const auto D12NX3 = D12 * NX3;
	const auto D12NX4 = D12 * NX4;

	const auto D13NX1 = D13 * NX1;
	const auto D13NX2 = D13 * NX2;
	const auto D13NX3 = D13 * NX3;
	const auto D13NX4 = D13 * NX4;

	const auto D21NY1 = D21 * NY1;
	const auto D21NY2 = D21 * NY2;
	const auto D21NY3 = D21 * NY3;
	const auto D21NY4 = D21 * NY4;

	const auto D22NY1 = D22 * NY1;
	const auto D22NY2 = D22 * NY2;
	const auto D22NY3 = D22 * NY3;
	const auto D22NY4 = D22 * NY4;

	const auto D23NY1 = D23 * NY1;
	const auto D23NY2 = D23 * NY2;
	const auto D23NY3 = D23 * NY3;
	const auto D23NY4 = D23 * NY4;

	const auto D31NX1 = D31 * NX1;
	const auto D31NX2 = D31 * NX2;
	const auto D31NX3 = D31 * NX3;
	const auto D31NX4 = D31 * NX4;
	const auto D31NY1 = D31 * NY1;
	const auto D31NY2 = D31 * NY2;
	const auto D31NY3 = D31 * NY3;
	const auto D31NY4 = D31 * NY4;

	const auto D32NX1 = D32 * NX1;
	const auto D32NX2 = D32 * NX2;
	const auto D32NX3 = D32 * NX3;
	const auto D32NX4 = D32 * NX4;
	const auto D32NY1 = D32 * NY1;
	const auto D32NY2 = D32 * NY2;
	const auto D32NY3 = D32 * NY3;
	const auto D32NY4 = D32 * NY4;

	const auto D33NX1 = D33 * NX1;
	const auto D33NX2 = D33 * NX2;
	const auto D33NX3 = D33 * NX3;
	const auto D33NX4 = D33 * NX4;
	const auto D33NY1 = D33 * NY1;
	const auto D33NY2 = D33 * NY2;
	const auto D33NY3 = D33 * NY3;
	const auto D33NY4 = D33 * NY4;

	const auto D11NX1D31NY1 = D11NX1 + D31NY1;
	const auto D13NX1D33NY1 = D13NX1 + D33NY1;
	const auto D12NX1D32NY1 = D12NX1 + D32NY1;
	const auto D31NX1D21NY1 = D31NX1 + D21NY1;
	const auto D33NX1D23NY1 = D33NX1 + D23NY1;
	const auto D32NX1D22NY1 = D32NX1 + D22NY1;
	const auto D11NX2D31NY2 = D11NX2 + D31NY2;
	const auto D13NX2D33NY2 = D13NX2 + D33NY2;
	const auto D12NX2D32NY2 = D12NX2 + D32NY2;
	const auto D31NX2D21NY2 = D31NX2 + D21NY2;
	const auto D33NX2D23NY2 = D33NX2 + D23NY2;
	const auto D32NX2D22NY2 = D32NX2 + D22NY2;
	const auto D11NX3D31NY3 = D11NX3 + D31NY3;
	const auto D13NX3D33NY3 = D13NX3 + D33NY3;
	const auto D12NX3D32NY3 = D12NX3 + D32NY3;
	const auto D31NX3D21NY3 = D31NX3 + D21NY3;
	const auto D33NX3D23NY3 = D33NX3 + D23NY3;
	const auto D32NX3D22NY3 = D32NX3 + D22NY3;
	const auto D11NX4D31NY4 = D11NX4 + D31NY4;
	const auto D13NX4D33NY4 = D13NX4 + D33NY4;
	const auto D12NX4D32NY4 = D12NX4 + D32NY4;
	const auto D31NX4D21NY4 = D31NX4 + D21NY4;
	const auto D33NX4D23NY4 = D33NX4 + D23NY4;
	const auto D32NX4D22NY4 = D32NX4 + D22NY4;

	K(0, 0) += NX1 * D11NX1D31NY1 + NY1 * D13NX1D33NY1;
	K(0, 1) += NX1 * D13NX1D33NY1 + NY1 * D12NX1D32NY1;
	K(0, 2) += NX2 * D11NX1D31NY1 + NY2 * D13NX1D33NY1;
	K(0, 3) += NX2 * D13NX1D33NY1 + NY2 * D12NX1D32NY1;
	K(0, 4) += NX3 * D11NX1D31NY1 + NY3 * D13NX1D33NY1;
	K(0, 5) += NX3 * D13NX1D33NY1 + NY3 * D12NX1D32NY1;
	K(0, 6) += NX4 * D11NX1D31NY1 + NY4 * D13NX1D33NY1;
	K(0, 7) += NX4 * D13NX1D33NY1 + NY4 * D12NX1D32NY1;
	K(1, 0) += NX1 * D31NX1D21NY1 + NY1 * D33NX1D23NY1;
	K(1, 1) += NX1 * D33NX1D23NY1 + NY1 * D32NX1D22NY1;
	K(1, 2) += NX2 * D31NX1D21NY1 + NY2 * D33NX1D23NY1;
	K(1, 3) += NX2 * D33NX1D23NY1 + NY2 * D32NX1D22NY1;
	K(1, 4) += NX3 * D31NX1D21NY1 + NY3 * D33NX1D23NY1;
	K(1, 5) += NX3 * D33NX1D23NY1 + NY3 * D32NX1D22NY1;
	K(1, 6) += NX4 * D31NX1D21NY1 + NY4 * D33NX1D23NY1;
	K(1, 7) += NX4 * D33NX1D23NY1 + NY4 * D32NX1D22NY1;
	K(2, 0) += NX1 * D11NX2D31NY2 + NY1 * D13NX2D33NY2;
	K(2, 1) += NX1 * D13NX2D33NY2 + NY1 * D12NX2D32NY2;
	K(2, 2) += NX2 * D11NX2D31NY2 + NY2 * D13NX2D33NY2;
	K(2, 3) += NX2 * D13NX2D33NY2 + NY2 * D12NX2D32NY2;
	K(2, 4) += NX3 * D11NX2D31NY2 + NY3 * D13NX2D33NY2;
	K(2, 5) += NX3 * D13NX2D33NY2 + NY3 * D12NX2D32NY2;
	K(2, 6) += NX4 * D11NX2D31NY2 + NY4 * D13NX2D33NY2;
	K(2, 7) += NX4 * D13NX2D33NY2 + NY4 * D12NX2D32NY2;
	K(3, 0) += NX1 * D31NX2D21NY2 + NY1 * D33NX2D23NY2;
	K(3, 1) += NX1 * D33NX2D23NY2 + NY1 * D32NX2D22NY2;
	K(3, 2) += NX2 * D31NX2D21NY2 + NY2 * D33NX2D23NY2;
	K(3, 3) += NX2 * D33NX2D23NY2 + NY2 * D32NX2D22NY2;
	K(3, 4) += NX3 * D31NX2D21NY2 + NY3 * D33NX2D23NY2;
	K(3, 5) += NX3 * D33NX2D23NY2 + NY3 * D32NX2D22NY2;
	K(3, 6) += NX4 * D31NX2D21NY2 + NY4 * D33NX2D23NY2;
	K(3, 7) += NX4 * D33NX2D23NY2 + NY4 * D32NX2D22NY2;
	K(4, 0) += NX1 * D11NX3D31NY3 + NY1 * D13NX3D33NY3;
	K(4, 1) += NX1 * D13NX3D33NY3 + NY1 * D12NX3D32NY3;
	K(4, 2) += NX2 * D11NX3D31NY3 + NY2 * D13NX3D33NY3;
	K(4, 3) += NX2 * D13NX3D33NY3 + NY2 * D12NX3D32NY3;
	K(4, 4) += NX3 * D11NX3D31NY3 + NY3 * D13NX3D33NY3;
	K(4, 5) += NX3 * D13NX3D33NY3 + NY3 * D12NX3D32NY3;
	K(4, 6) += NX4 * D11NX3D31NY3 + NY4 * D13NX3D33NY3;
	K(4, 7) += NX4 * D13NX3D33NY3 + NY4 * D12NX3D32NY3;
	K(5, 0) += NX1 * D31NX3D21NY3 + NY1 * D33NX3D23NY3;
	K(5, 1) += NX1 * D33NX3D23NY3 + NY1 * D32NX3D22NY3;
	K(5, 2) += NX2 * D31NX3D21NY3 + NY2 * D33NX3D23NY3;
	K(5, 3) += NX2 * D33NX3D23NY3 + NY2 * D32NX3D22NY3;
	K(5, 4) += NX3 * D31NX3D21NY3 + NY3 * D33NX3D23NY3;
	K(5, 5) += NX3 * D33NX3D23NY3 + NY3 * D32NX3D22NY3;
	K(5, 6) += NX4 * D31NX3D21NY3 + NY4 * D33NX3D23NY3;
	K(5, 7) += NX4 * D33NX3D23NY3 + NY4 * D32NX3D22NY3;
	K(6, 0) += NX1 * D11NX4D31NY4 + NY1 * D13NX4D33NY4;
	K(6, 1) += NX1 * D13NX4D33NY4 + NY1 * D12NX4D32NY4;
	K(6, 2) += NX2 * D11NX4D31NY4 + NY2 * D13NX4D33NY4;
	K(6, 3) += NX2 * D13NX4D33NY4 + NY2 * D12NX4D32NY4;
	K(6, 4) += NX3 * D11NX4D31NY4 + NY3 * D13NX4D33NY4;
	K(6, 5) += NX3 * D13NX4D33NY4 + NY3 * D12NX4D32NY4;
	K(6, 6) += NX4 * D11NX4D31NY4 + NY4 * D13NX4D33NY4;
	K(6, 7) += NX4 * D13NX4D33NY4 + NY4 * D12NX4D32NY4;
	K(7, 0) += NX1 * D31NX4D21NY4 + NY1 * D33NX4D23NY4;
	K(7, 1) += NX1 * D33NX4D23NY4 + NY1 * D32NX4D22NY4;
	K(7, 2) += NX2 * D31NX4D21NY4 + NY2 * D33NX4D23NY4;
	K(7, 3) += NX2 * D33NX4D23NY4 + NY2 * D32NX4D22NY4;
	K(7, 4) += NX3 * D31NX4D21NY4 + NY3 * D33NX4D23NY4;
	K(7, 5) += NX3 * D33NX4D23NY4 + NY3 * D32NX4D22NY4;
	K(7, 6) += NX4 * D31NX4D21NY4 + NY4 * D33NX4D23NY4;
	K(7, 7) += NX4 * D33NX4D23NY4 + NY4 * D32NX4D22NY4;
}

CP4::CP4(const unsigned T, uvec&& N, const unsigned M, const double TH, const bool R, const bool F)
	: MaterialElement(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, F)
	, thickness(TH)
	, reduced_scheme(R) {}

void CP4::initialize(const shared_ptr<DomainBase>& D) {
	const auto material_proto = std::dynamic_pointer_cast<Material2D>(D->get_material(unsigned(material_tag(0))));

	if(material_proto == nullptr) {
		D->disable_element(get_tag());
		return;
	}

	if(material_proto->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

	auto ele_coor = get_coordinate(2);

	if(reduced_scheme) {
		hourglassing.zeros(m_size, m_size);
		const auto area = .5 * ((ele_coor(2, 0) - ele_coor(0, 0)) * (ele_coor(3, 1) - ele_coor(1, 1)) + (ele_coor(1, 0) - ele_coor(3, 0)) * (ele_coor(2, 1) - ele_coor(0, 1)));
		vec b1(4), b2(4);
		b1(0) = ele_coor(1, 1) - ele_coor(3, 1);
		b1(1) = ele_coor(2, 1) - ele_coor(0, 1);
		b1(2) = ele_coor(3, 1) - ele_coor(1, 1);
		b1(3) = ele_coor(0, 1) - ele_coor(2, 1);
		b2(0) = ele_coor(3, 0) - ele_coor(1, 0);
		b2(1) = ele_coor(0, 0) - ele_coor(2, 0);
		b2(2) = ele_coor(1, 0) - ele_coor(3, 0);
		b2(3) = ele_coor(2, 0) - ele_coor(0, 0);
		const vec h{1., -1., 1., -1.};
		vec gamma = 2. * area * h - dot(h, ele_coor.col(0)) * b1 - dot(h, ele_coor.col(1)) * b2;
		mat t_hourglassing = gamma * gamma.t();
		for(unsigned I = 0, K = 0; I < m_node; ++I, K += m_dof) for(unsigned J = 0, L = 0; J < m_node; ++J, L += m_dof) hourglassing(K + 1, L + 1) = hourglassing(K, L) = t_hourglassing(I, J);
	}

	auto& ini_stiffness = material_proto->get_initial_stiffness();

	const IntegrationPlan plan(2, reduced_scheme ? 1 : 2, IntegrationType::GAUSS);

	initial_stiffness.zeros(m_size, m_size);

	int_pt.clear(), int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		vec t_vec{plan(I, 0), plan(I, 1)};
		const auto pn = shape::quad(t_vec, 1);
		const mat jacob = pn * ele_coor;
		int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob), material_proto->get_copy(), solve(jacob, pn));

		auto& c_pt = int_pt.back();

		for(unsigned J = 0; J < m_node; ++J) {
			const auto K = m_dof * J;
			c_pt.strain_mat(0, K) = c_pt.strain_mat(2, K + 1) = c_pt.pn_pxy(0, J);
			c_pt.strain_mat(2, K) = c_pt.strain_mat(1, K + 1) = c_pt.pn_pxy(1, J);
		}
		initial_stiffness += c_pt.weight * thickness * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;
	}
	trial_stiffness = current_stiffness = initial_stiffness;

	initial_mass.zeros(m_size, m_size);
	const auto tmp_density = material_proto->get_parameter();
	if(tmp_density != 0.) {
		for(const auto& I : int_pt) {
			const auto n_int = shape::quad(I.coor, 0);
			const auto tmp_a = tmp_density * I.weight * thickness;
			for(unsigned J = 0; J < m_node; ++J) for(auto K = J; K < m_node; ++K) initial_mass(m_dof * J, m_dof * K) += tmp_a * n_int(J) * n_int(K);
		}
		for(unsigned I = 0; I < m_node * m_dof; I += m_dof) {
			initial_mass(I + 1, I + 1) = initial_mass(I, I);
			for(auto J = I + m_dof; J < m_size; J += m_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
		}
	}
	trial_mass = current_mass = initial_mass;
}

int CP4::update_status() {
	const auto t_disp = get_trial_displacement();

	trial_stiffness.zeros(m_size, m_size);
	trial_resistance.zeros(m_size);

	if(nlgeom) {
		const mat ele_disp = reshape(t_disp, m_dof, m_node);

		trial_geometry.zeros(m_size, m_size);

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
				const auto K = J * m_dof;
				t_strain(0) += t_disp(K) * I.pn_pxy(0, J);
				t_strain(1) += t_disp(K + 1) * I.pn_pxy(1, J);
				t_strain(2) += t_disp(K) * I.pn_pxy(1, J) + t_disp(K + 1) * I.pn_pxy(0, J);
			}

			if(I.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

			const auto t_factor = I.weight * thickness;

			stack_stiffness(trial_stiffness, I.m_material->get_trial_stiffness(), I.strain_mat, t_factor);
			trial_resistance += t_factor * I.strain_mat.t() * I.m_material->get_trial_stress();
		}

	if(reduced_scheme) {
		trial_stiffness += hourglassing;
		trial_resistance += hourglassing * t_disp;
	}

	return SUANPAN_SUCCESS;
}

int CP4::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	return code;
}

int CP4::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	return code;
}

int CP4::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	return code;
}

vector<vec> CP4::record(const OutputType P) {
	vector<vec> output;
	output.reserve(int_pt.size());

	for(const auto& I : int_pt) for(const auto& J : I.m_material->record(P)) output.emplace_back(J);

	return output;
}

void CP4::print() {
	suanpan_info("Element %u is a four-node membrane element (CP4)%s.\n", get_tag(), nlgeom ? " with nonlinear geomotry (TL formulation)" : "");
	suanpan_info("The nodes connected are:\n");
	node_encoding.t().print();
	suanpan_info("Material model response:\n");
	for(size_t I = 0; I < int_pt.size(); ++I) {
		suanpan_info("Integration Point %lu:\t", I + 1);
		int_pt[I].coor.t().print();
		int_pt[I].m_material->print();
	}
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void CP4::Setup() {
	vtk_cell = vtkSmartPointer<vtkQuad>::New();
	auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < m_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void CP4::GetDisplacement(vtkSmartPointer<vtkDoubleArray>& arrays) {
	mat t_disp(6, m_node, fill::zeros);
	t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);
	for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

mat CP4::GetData(const OutputType P) {
	if(P == OutputType::S) {
		mat B(3 * int_pt.size(), 1);
		mat A(3 * int_pt.size(), 7);
		for(unsigned I = 0, J = 0; I < unsigned(int_pt.size()); ++I, J = J + 3) {
			B.rows(J, J + 2) = int_pt[I].m_material->get_current_stress();
			A.rows(J, J + 2) = shape::stress7(int_pt[I].coor);
		}
		mat X;
		if(!solve(X, A, B)) return {};

		mat t_stress(6, m_node, fill::zeros);

		t_stress(uvec{0, 1, 3}, uvec{0}) = shape::stress7(-1., -1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{1}) = shape::stress7(1., -1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{2}) = shape::stress7(1., 1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{3}) = shape::stress7(-1., 1.) * X;

		return t_stress;
	}

	return {};
}

void CP4::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
	for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
