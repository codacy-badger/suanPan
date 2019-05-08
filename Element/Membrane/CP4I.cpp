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

#include "CP4I.h"
#include <Domain/DomainBase.h>
#include <Element/Utility/MatrixModifier.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>

const unsigned CP4I::m_node = 4;
const unsigned CP4I::m_dof = 2;
const unsigned CP4I::m_size = m_dof * m_node;

CP4I::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& PNPXY)
	: coor(std::forward<vec>(C))
	, weight(W)
	, m_material(std::forward<unique_ptr<Material>>(M))
	, pn_pxy(std::forward<mat>(PNPXY))
	, B1(3, m_size, fill::zeros)
	, B2(3, 4, fill::zeros) {}

void CP4I::stack_stiffness(mat& K, const mat& D, const mat& N, const double F) {
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

void CP4I::stack_stiffness_incompatible(mat& K, const mat& D, const mat& N, const double F) {
	const auto& D11 = D(0, 0);
	const auto& D12 = D(0, 1);
	const auto& D13 = D(0, 2);
	const auto& D21 = D(1, 0);
	const auto& D22 = D(1, 1);
	const auto& D23 = D(1, 2);
	const auto& D31 = D(2, 0);
	const auto& D32 = D(2, 1);
	const auto& D33 = D(2, 2);

	const auto& NX = N(0, 0);
	const auto& NY = N(1, 3);

	const auto NXNX = F * NX * NX;
	const auto NXNY = F * NX * NY;
	const auto NYNY = F * NY * NY;

	K(0) += D11 * NXNX;
	K(1) += D31 * NXNX;
	K(2) += D31 * NXNY;
	K(3) += D21 * NXNY;
	K(4) += D13 * NXNX;
	K(5) += D33 * NXNX;
	K(6) += D33 * NXNY;
	K(7) += D23 * NXNY;
	K(8) += D13 * NXNY;
	K(9) += D33 * NXNY;
	K(10) += D33 * NYNY;
	K(11) += D23 * NYNY;
	K(12) += D12 * NXNY;
	K(13) += D32 * NXNY;
	K(14) += D32 * NYNY;
	K(15) += D22 * NYNY;
}

CP4I::CP4I(const unsigned T, uvec&& N, const unsigned M, const double TH)
	: MaterialElement(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, false)
	, thickness(TH) {}

void CP4I::initialize(const shared_ptr<DomainBase>& D) {
	const auto material_proto = std::dynamic_pointer_cast<Material2D>(D->get_material(unsigned(material_tag(0))));

	if(material_proto == nullptr) {
		D->disable_element(get_tag());
		return;
	}

	if(material_proto->plane_type == PlaneType::E) access::rw(thickness) = 1.;

	auto& ini_stiffness = material_proto->get_initial_stiffness();

	const auto ele_coor = get_coordinate(2);

	const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

	initial_stiffness.zeros(m_size, m_size);
	mat stiff_a(4, 4, fill::zeros), stiff_b(m_size, 4, fill::zeros);

	int_pt.clear(), int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		vec t_vec{plan(I, 0), plan(I, 1)};
		const auto pn = shape::quad(t_vec, 1);
		const mat jacob = pn * ele_coor;
		int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob) * thickness, material_proto->get_copy(), solve(jacob, pn));

		auto& c_pt = int_pt.back();

		for(unsigned J = 0; J < m_node; ++J) {
			const auto K = m_dof * J;
			c_pt.B1(0, K) = c_pt.B1(2, K + 1) = c_pt.pn_pxy(0, J);
			c_pt.B1(2, K) = c_pt.B1(1, K + 1) = c_pt.pn_pxy(1, J);
		}

		const vec pbn_pxy = solve(jacob, -2. * c_pt.coor);
		c_pt.B2(0, 0) = c_pt.B2(2, 1) = pbn_pxy(0);
		c_pt.B2(1, 3) = c_pt.B2(2, 2) = pbn_pxy(1);

		stack_stiffness(initial_stiffness, ini_stiffness, c_pt.B1, c_pt.weight);
		stack_stiffness_incompatible(stiff_a, ini_stiffness, c_pt.B2, c_pt.weight);
		stiff_b += c_pt.B1.t() * ini_stiffness * c_pt.B2 * c_pt.weight;
	}
	initial_stiffness -= stiff_b * solve(stiff_a, stiff_b.t());
	trial_stiffness = current_stiffness = initial_stiffness;

	initial_mass.zeros(m_size, m_size);
	const auto tmp_density = material_proto->get_parameter();
	if(tmp_density != 0.) {
		for(const auto& I : int_pt) {
			const auto n_int = shape::quad(I.coor, 0);
			const auto tmp_a = tmp_density * I.weight;
			for(unsigned J = 0; J < m_node; ++J) for(auto K = J; K < m_node; ++K) initial_mass(m_dof * J, m_dof * K) += tmp_a * n_int(J) * n_int(K);
		}
		for(unsigned I = 0; I < m_node * m_dof; I += m_dof) {
			initial_mass(I + 1, I + 1) = initial_mass(I, I);
			for(auto J = I + m_dof; J < m_size; J += m_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
		}
	}
	trial_mass = current_mass = initial_mass;
}

int CP4I::update_status() {
	const auto t_disp = get_trial_displacement();

	trial_stiffness.zeros(m_size, m_size);
	trial_resistance.zeros(m_size);
	mat stiff_a(4, 4, fill::zeros), stiff_b(m_size, 4, fill::zeros);
	vec resistance_a(4, fill::zeros);

	for(const auto& I : int_pt) {
		vec t_strain(3, fill::zeros);
		for(unsigned J = 0; J < m_node; ++J) {
			const auto K = J * m_dof;
			t_strain(0) += t_disp(K) * I.pn_pxy(0, J);
			t_strain(1) += t_disp(K + 1) * I.pn_pxy(1, J);
			t_strain(2) += t_disp(K) * I.pn_pxy(1, J) + t_disp(K + 1) * I.pn_pxy(0, J);
		}

		if(I.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		stack_stiffness(trial_stiffness, I.m_material->get_trial_stiffness(), I.B1, I.weight);
		trial_resistance += I.weight * I.B1.t() * I.m_material->get_trial_stress();
		stack_stiffness_incompatible(stiff_a, I.m_material->get_trial_stiffness(), I.B2, I.weight);
		stiff_b += I.B1.t() * I.m_material->get_trial_stiffness() * I.B2 * I.weight;
		resistance_a += I.weight * I.B2.t() * I.m_material->get_trial_stress();
	}
	trial_stiffness -= stiff_b * solve(stiff_a, stiff_b.t());
	trial_resistance -= stiff_b * solve(stiff_a, resistance_a);

	return SUANPAN_SUCCESS;
}

int CP4I::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	return code;
}

int CP4I::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	return code;
}

int CP4I::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	return code;
}

vector<vec> CP4I::record(const OutputType P) {
	vector<vec> output;
	output.reserve(int_pt.size());

	for(const auto& I : int_pt) for(const auto& J : I.m_material->record(P)) output.emplace_back(J);

	return output;
}

void CP4I::print() {
	suanpan_info("Element %u is a four-node membrane element (CP4I).\n", get_tag());
	suanpan_info("The nodes connected are:\n");
	node_encoding.t().print();
	suanpan_info("Material model response:\n");
	for(size_t I = 0; I < int_pt.size(); ++I) {
		suanpan_info("Integration Point %lu:\t", I + 1);
		int_pt[I].coor.t().print();
		int_pt[I].m_material->print();
	}
}
