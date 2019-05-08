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

#include "C3D8I.h"
#include <Domain/DomainBase.h>
#include <Material/Material3D/Material3D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/tensorToolbox.h>

const unsigned C3D8I::c_node = 8;
const unsigned C3D8I::c_dof = 3;
const unsigned C3D8I::c_size = c_dof * c_node;

C3D8I::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& P)
	: coor(std::forward<vec>(C))
	, weight(W)
	, c_material(std::forward<unique_ptr<Material>>(M))
	, pn_pxyz(std::forward<mat>(P))
	, B1(6, c_size, fill::zeros)
	, B2(6, 9, fill::zeros) {}

C3D8I::C3D8I(const unsigned T, uvec&& N, const unsigned M)
	: MaterialElement(T, c_node, c_dof, std::forward<uvec>(N), uvec{M}, false) {}

void C3D8I::initialize(const shared_ptr<DomainBase>& D) {
	const auto ele_coor = get_coordinate(c_dof);

	auto& mat_proto = D->get_material(unsigned(material_tag(0)));
	auto& mat_stiff = mat_proto->get_initial_stiffness();

	const IntegrationPlan plan(3, 2, IntegrationType::IRONS);

	initial_stiffness.zeros(c_size, c_size);
	mat stiff_a(9, 9, fill::zeros), stiff_b(c_size, 9, fill::zeros);

	int_pt.clear(), int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		vec t_vec{plan(I, 0), plan(I, 1), plan(I, 2)};
		const auto pn = shape::cube(t_vec, 1);
		const mat jacob = pn * ele_coor;
		int_pt.emplace_back(std::move(t_vec), plan(I, c_dof) * det(jacob), mat_proto->get_copy(), solve(jacob, pn));

		auto& c_pt = int_pt.back();
		for(unsigned J = 0; J < c_node; ++J) {
			const auto K = c_dof * J;
			c_pt.B1(0, K) = c_pt.B1(3, K + 1) = c_pt.B1(5, K + 2) = c_pt.pn_pxyz(0, J);
			c_pt.B1(3, K) = c_pt.B1(1, K + 1) = c_pt.B1(4, K + 2) = c_pt.pn_pxyz(1, J);
			c_pt.B1(5, K) = c_pt.B1(4, K + 1) = c_pt.B1(2, K + 2) = c_pt.pn_pxyz(2, J);
		}

		const vec pbn_pxyz = solve(jacob, -2. * c_pt.coor);
		c_pt.B2(0, 0) = c_pt.B2(3, 1) = c_pt.B2(5, 2) = pbn_pxyz(0);
		c_pt.B2(1, 4) = c_pt.B2(3, 3) = c_pt.B2(4, 5) = pbn_pxyz(1);
		c_pt.B2(2, 8) = c_pt.B2(4, 7) = c_pt.B2(5, 6) = pbn_pxyz(2);

		initial_stiffness += c_pt.weight * c_pt.B1.t() * mat_stiff * c_pt.B1;

		const auto t_stiff = mat_stiff * c_pt.B2 * c_pt.weight;
		stiff_a += c_pt.B2.t() * t_stiff;
		stiff_b += c_pt.B1.t() * t_stiff;
	}
	initial_stiffness -= stiff_b * solve(stiff_a, stiff_b.t());
	trial_stiffness = current_stiffness = initial_stiffness;

	initial_mass.zeros(c_size, c_size);
	const auto t_density = mat_proto->get_parameter();
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

int C3D8I::update_status() {
	const auto t_disp = get_trial_displacement();

	trial_stiffness.zeros(c_size, c_size);
	trial_resistance.zeros(c_size);
	mat stiff_a(9, 9, fill::zeros), stiff_b(c_size, 9, fill::zeros);
	vec resistance_a(9, fill::zeros);

	for(const auto& I : int_pt) {
		if(I.c_material->update_trial_status(I.B1 * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		trial_stiffness += I.weight * I.B1.t() * I.c_material->get_trial_stiffness() * I.B1;
		trial_resistance += I.weight * I.B1.t() * I.c_material->get_trial_stress();
		const auto t_stiff = I.c_material->get_trial_stiffness() * I.B2 * I.weight;
		stiff_a += I.B2.t() * t_stiff;
		stiff_b += I.B1.t() * t_stiff;
		resistance_a += I.weight * I.B2.t() * I.c_material->get_trial_stress();
	}
	trial_stiffness -= stiff_b * solve(stiff_a, stiff_b.t());
	trial_resistance -= stiff_b * solve(stiff_a, resistance_a);

	return SUANPAN_SUCCESS;
}

int C3D8I::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.c_material->commit_status();
	return code;
}

int C3D8I::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.c_material->clear_status();
	return code;
}

int C3D8I::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.c_material->reset_status();
	return code;
}

vector<vec> C3D8I::record(const OutputType T) {
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

void C3D8I::print() {
	suanpan_info("C3D8I element connects nodes:");
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
