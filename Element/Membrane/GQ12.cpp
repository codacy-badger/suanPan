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

#include "GQ12.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>

const unsigned GQ12::m_node = 4;
const unsigned GQ12::m_dof = 3;
const unsigned GQ12::m_size = m_dof * m_node;

GQ12::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
	: coor(std::forward<vec>(C))
	, weight(W)
	, m_material(std::forward<unique_ptr<Material>>(M))
	, strain_mat(3, m_size, fill::zeros) {}

GQ12::GQ12(const unsigned T, uvec&& N, const unsigned M, const double TH)
	: MaterialElement(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, false)
	, thickness(TH) {}

void GQ12::initialize(const shared_ptr<DomainBase>& D) {
	const auto mat_proto = std::dynamic_pointer_cast<Material2D>(D->get_material(unsigned(material_tag(0))));

	if(mat_proto == nullptr) {
		D->disable_element(get_tag());
		return;
	}

	if(mat_proto->plane_type == PlaneType::E) access::rw(thickness) = 1.;

	auto& mat_stiff = mat_proto->get_initial_stiffness();

	const auto ele_coor = get_coordinate(2);

	const auto LX1 = ele_coor(1, 1) - ele_coor(0, 1);
	const auto LX2 = ele_coor(2, 1) - ele_coor(1, 1);
	const auto LX3 = ele_coor(3, 1) - ele_coor(2, 1);
	const auto LX4 = ele_coor(0, 1) - ele_coor(3, 1);
	const auto LY1 = ele_coor(0, 0) - ele_coor(1, 0);
	const auto LY2 = ele_coor(1, 0) - ele_coor(2, 0);
	const auto LY3 = ele_coor(2, 0) - ele_coor(3, 0);
	const auto LY4 = ele_coor(3, 0) - ele_coor(0, 0);

	const IntegrationPlan plan(2, 2, IntegrationType::IRONS);

	mat pnt(2, 8);

	initial_stiffness.zeros(m_size, m_size);

	int_pt.clear(), int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		vec t_vec{plan(I, 0), plan(I, 1)};
		const auto pn = shape::quad(t_vec, 1);
		const mat jacob = pn * ele_coor;
		const mat pn_pxy = solve(jacob, pn);
		int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob) * thickness, mat_proto->get_copy());

		const auto TX = 2. * plan(I, 0);
		const auto TY = 2. * plan(I, 1);

		const auto AA = plan(I, 0) + 1.;
		const auto BB = plan(I, 0) - 1.;
		const auto CC = plan(I, 1) + 1.;
		const auto DD = plan(I, 1) - 1.;

		pnt(0, 0) = +DD * (LX4 * CC - LX1 * TX);
		pnt(0, 1) = +DD * (LX2 * CC + LX1 * TX);
		pnt(0, 2) = -CC * (LX2 * DD - LX3 * TX);
		pnt(0, 3) = -CC * (LX4 * DD + LX3 * TX);
		pnt(0, 4) = +DD * (LY4 * CC - LY1 * TX);
		pnt(0, 5) = +DD * (LY2 * CC + LY1 * TX);
		pnt(0, 6) = -CC * (LY2 * DD - LY3 * TX);
		pnt(0, 7) = -CC * (LY4 * DD + LY3 * TX);
		pnt(1, 0) = -BB * (LX1 * AA - LX4 * TY);
		pnt(1, 1) = +AA * (LX1 * BB + LX2 * TY);
		pnt(1, 2) = +AA * (LX3 * BB - LX2 * TY);
		pnt(1, 3) = -BB * (LX3 * AA + LX4 * TY);
		pnt(1, 4) = -BB * (LY1 * AA - LY4 * TY);
		pnt(1, 5) = +AA * (LY1 * BB + LY2 * TY);
		pnt(1, 6) = +AA * (LY3 * BB - LY2 * TY);
		pnt(1, 7) = -BB * (LY3 * AA + LY4 * TY);

		const mat pnt_pxy = solve(jacob, 625E-4 * pnt);

		auto& c_pt = int_pt.back();

		for(unsigned J = 0, K = 0, L = 1, M = 2, N = 4; J < m_node; ++J, K += m_dof, L += m_dof, M += m_dof, ++N) {
			c_pt.strain_mat(0, K) = c_pt.strain_mat(2, L) = pn_pxy(0, J);
			c_pt.strain_mat(2, K) = c_pt.strain_mat(1, L) = pn_pxy(1, J);
			c_pt.strain_mat(0, M) = pnt_pxy(0, J);
			c_pt.strain_mat(1, M) = pnt_pxy(1, N);
			c_pt.strain_mat(2, M) = pnt_pxy(0, N) + pnt_pxy(1, J);
		}

		initial_stiffness += c_pt.strain_mat.t() * mat_stiff * c_pt.strain_mat * c_pt.weight;
	}
	trial_stiffness = current_stiffness = initial_stiffness;
}

int GQ12::update_status() {
	const auto t_disp = get_trial_displacement();

	trial_stiffness.zeros(m_size, m_size);
	trial_resistance.zeros(m_size);
	for(const auto& I : int_pt) {
		if(I.m_material->update_trial_status(I.strain_mat * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		trial_stiffness += I.weight * I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat;
		trial_resistance += I.weight * I.strain_mat.t() * I.m_material->get_trial_stress();
	}

	return SUANPAN_SUCCESS;
}

int GQ12::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	return code;
}

int GQ12::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	return code;
}

int GQ12::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	return code;
}

void GQ12::print() {
	suanpan_info("material model response:\n");
	for(size_t I = 0; I < int_pt.size(); ++I) {
		suanpan_info("Integration Point %lu:\n", I + 1);
		int_pt[I].m_material->print();
	}
}
