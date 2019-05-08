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

#include "SRGCMQ.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>

const unsigned SRGCMQ::m_node = 4;
const unsigned SRGCMQ::m_dof = 3;
const unsigned SRGCMQ::m_size = m_dof * m_node;
const mat SRGCMQ::iso_mapping;

SRGCMQ::IntegrationPoint::IntegrationPoint(vec&& C, const double F, unique_ptr<Material>&& M)
	: coor(std::forward<vec>(C))
	, factor(F)
	, m_material(std::forward<unique_ptr<Material>>(M)) {}

mat SRGCMQ::form_transformation(const mat& jacobian) {
	mat trans_mat(3, 3);

	trans_mat(0, 0) = jacobian(0, 0) * jacobian(0, 0);
	trans_mat(1, 0) = jacobian(0, 1) * jacobian(0, 1);
	trans_mat(2, 0) = jacobian(0, 0) * jacobian(0, 1);

	trans_mat(0, 1) = jacobian(1, 0) * jacobian(1, 0);
	trans_mat(1, 1) = jacobian(1, 1) * jacobian(1, 1);
	trans_mat(2, 1) = jacobian(1, 0) * jacobian(1, 1);

	trans_mat(0, 2) = 2. * jacobian(0, 0) * jacobian(1, 0);
	trans_mat(1, 2) = 2. * jacobian(1, 0) * jacobian(1, 1);
	trans_mat(2, 2) = jacobian(0, 0) * jacobian(1, 1) + jacobian(0, 1) * jacobian(1, 0);

	return trans_mat / accu(square(trans_mat));
}

mat SRGCMQ::form_drilling_mass(const vec& coor, const vec& lxy) {
	mat poly_mass(2, m_size, fill::zeros);

	auto &X = coor(0), &Y = coor(1);

	auto &LX1 = lxy(0), &LX2 = lxy(1), &LX3 = lxy(2), &LX4 = lxy(3);
	auto &LY1 = lxy(4), &LY2 = lxy(5), &LY3 = lxy(6), &LY4 = lxy(7);

	const auto XX = X * X, YY = Y * Y, XXY = XX * Y, YYX = YY * X;

	poly_mass(0, 2) = -LX1 * XXY + LX1 * XX + LX4 * YYX - LX4 * X - LX4 * YY + LX1 * Y - LX1 + LX4;
	poly_mass(0, 5) = LX1 * XXY - LX1 * XX + LX2 * YYX - LX2 * X + LX2 * YY - LX1 * Y + LX1 - LX2;
	poly_mass(0, 8) = LX3 * XXY + LX3 * XX - LX2 * YYX + LX2 * X - LX2 * YY - LX3 * Y + LX2 - LX3;
	poly_mass(0, 11) = -LX3 * XXY - LX3 * XX - LX4 * YYX + LX4 * X + LX4 * YY + LX3 * Y + LX3 - LX4;
	poly_mass(1, 2) = -LY1 * XXY + LY1 * XX + LY4 * YYX - LY4 * X - LY4 * YY + LY1 * Y - LY1 + LY4;
	poly_mass(1, 5) = LY1 * XXY - LY1 * XX + LY2 * YYX - LY2 * X + LY2 * YY - LY1 * Y + LY1 - LY2;
	poly_mass(1, 8) = LY3 * XXY + LY3 * XX - LY2 * YYX + LY2 * X - LY2 * YY - LY3 * Y + LY2 - LY3;
	poly_mass(1, 11) = -LY3 * XXY - LY3 * XX - LY4 * YYX + LY4 * X + LY4 * YY + LY3 * Y + LY3 - LY4;

	poly_mass /= 16.;

	return poly_mass;
}

mat SRGCMQ::form_drilling_displacement(const vec& coor, const vec& lxy) {
	mat poly_drilling(2, 8);

	auto &X = coor(0), &Y = coor(1);

	auto &LX1 = lxy(0), &LX2 = lxy(1), &LX3 = lxy(2), &LX4 = lxy(3);
	auto &LY1 = lxy(4), &LY2 = lxy(5), &LY3 = lxy(6), &LY4 = lxy(7);

	const auto X2 = 2. * X, Y2 = 2. * Y, XP = X + 1., XM = X - 1., YP = Y + 1., YM = Y - 1.;

	poly_drilling(0, 0) = YM * (LX4 * YP - LX1 * X2);
	poly_drilling(0, 1) = YM * (LX2 * YP + LX1 * X2);
	poly_drilling(0, 2) = YP * (LX3 * X2 - LX2 * YM);
	poly_drilling(0, 3) = -YP * (LX3 * X2 + LX4 * YM);
	poly_drilling(0, 4) = YM * (LY4 * YP - LY1 * X2);
	poly_drilling(0, 5) = YM * (LY2 * YP + LY1 * X2);
	poly_drilling(0, 6) = YP * (LY3 * X2 - LY2 * YM);
	poly_drilling(0, 7) = -YP * (LY3 * X2 + LY4 * YM);
	poly_drilling(1, 0) = XM * (LX4 * Y2 - LX1 * XP);
	poly_drilling(1, 1) = XP * (LX1 * XM + LX2 * Y2);
	poly_drilling(1, 2) = XP * (LX3 * XM - LX2 * Y2);
	poly_drilling(1, 3) = -XM * (LX3 * XP + LX4 * Y2);
	poly_drilling(1, 4) = XM * (LY4 * Y2 - LY1 * XP);
	poly_drilling(1, 5) = XP * (LY1 * XM + LY2 * Y2);
	poly_drilling(1, 6) = XP * (LY3 * XM - LY2 * Y2);
	poly_drilling(1, 7) = -XM * (LY3 * XP + LY4 * Y2);

	poly_drilling /= 16.;

	return poly_drilling;
}

mat SRGCMQ::form_displacement(const mat& pn_pxy, const mat& pnt_pxy) {
	mat poly_disp(3, m_size, fill::zeros);

	for(unsigned J = 0, K = 0, L = 1, M = 2; J < m_node; ++J, K += m_dof, L += m_dof, M += m_dof) {
		poly_disp(0, K) = poly_disp(2, L) = pn_pxy(0, J);
		poly_disp(2, K) = poly_disp(1, L) = pn_pxy(1, J);
		poly_disp(0, M) = pnt_pxy(0, J);
		poly_disp(1, M) = pnt_pxy(1, J + 4);
		poly_disp(2, M) = pnt_pxy(0, J + 4) + pnt_pxy(1, J);
	}

	return poly_disp;
}

vec SRGCMQ::form_stress_mode(const double X, const double Y) { return vec{0., X, Y, X * Y}; }

SRGCMQ::SRGCMQ(const unsigned T, uvec&& N, const unsigned M, const double TH, const char IP, const double RHOX, const unsigned MATX, const double RHOY, const unsigned MATY)
	: MaterialElement(T, m_node, m_dof, std::forward<uvec>(N), uvec{M})
	, thickness(TH)
	, int_scheme(IP)
	, rho_x(RHOX)
	, rho_y(RHOY)
	, mat_x(MATX)
	, mat_y(MATY) {
	if(iso_mapping.is_empty()) {
		mat t_mapping(4, 4);
		t_mapping.fill(.25);
		t_mapping(1, 0) = t_mapping(1, 3) = t_mapping(2, 0) = t_mapping(2, 1) = t_mapping(3, 1) = t_mapping(3, 3) = -.25;
		access::rw(iso_mapping) = t_mapping;
	}
}

void SRGCMQ::initialize(const shared_ptr<DomainBase>& D) {
	const auto ele_coor = get_coordinate(2);

	const mat trans_mat = trans(iso_mapping * ele_coor);

	const auto jacob_trans = form_transformation(shape::quad(vec{0., 0.}, 1) * ele_coor);

	vec diff_coor(8);
	diff_coor(0) = ele_coor(1, 1) - ele_coor(0, 1);
	diff_coor(1) = ele_coor(2, 1) - ele_coor(1, 1);
	diff_coor(2) = ele_coor(3, 1) - ele_coor(2, 1);
	diff_coor(3) = ele_coor(0, 1) - ele_coor(3, 1);
	diff_coor(4) = ele_coor(0, 0) - ele_coor(1, 0);
	diff_coor(5) = ele_coor(1, 0) - ele_coor(2, 0);
	diff_coor(6) = ele_coor(2, 0) - ele_coor(3, 0);
	diff_coor(7) = ele_coor(3, 0) - ele_coor(0, 0);

	const auto mat_proto = std::dynamic_pointer_cast<Material2D>(D->get_material(unsigned(material_tag(0))));

	if(mat_proto == nullptr) {
		D->disable_element(get_tag());
		return;
	}

	if(mat_proto->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

	auto mat_stiff = mat_proto->get_initial_stiffness();

	if(reinforced_x) {
		if(D->find_material(mat_x) && std::dynamic_pointer_cast<Material1D>(D->get_material(mat_x)) != nullptr) mat_stiff(0, 0) += D->get_material(mat_x)->get_initial_stiffness().at(0) * rho_x;
		else {
			suanpan_error("cannot find the reinforcement material.\n");
			D->disable_element(get_tag());
			return;
		}
	}
	if(reinforced_y) {
		if(D->find_material(mat_y) && std::dynamic_pointer_cast<Material1D>(D->get_material(mat_y)) != nullptr) mat_stiff(1, 1) += D->get_material(mat_y)->get_initial_stiffness().at(0) * rho_y;
		else {
			suanpan_error("cannot find the reinforcement material.\n");
			D->disable_element(get_tag());
			return;
		}
	}

	const IntegrationPlan plan(2, int_scheme == 'I' ? 2 : 3, int_scheme == 'I' ? IntegrationType::IRONS : int_scheme == 'L' ? IntegrationType::LOBATTO : IntegrationType::GAUSS);

	mat N(11, 12, fill::zeros), H(11, 11, fill::zeros), HT(11, 11, fill::zeros);

	int_pt.clear(), int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		const auto &X = plan(I, 0), &Y = plan(I, 1);

		vec t_vec{X, Y};
		const auto pn = shape::quad(t_vec, 1);
		const mat jacob = pn * ele_coor;
		int_pt.emplace_back(std::move(t_vec), det(jacob) * plan(I, 2) * thickness, mat_proto->get_copy());

		auto& c_pt = int_pt.back();

		if(reinforced_x) c_pt.rebar_x = D->get_material(mat_x)->get_copy();
		if(reinforced_y) c_pt.rebar_y = D->get_material(mat_y)->get_copy();

		c_pt.poly_strain = shape::stress11(trans_mat * form_stress_mode(X, Y));
		const auto& poly_stress = c_pt.poly_strain;

		N += c_pt.factor * poly_stress.t() * form_displacement(solve(jacob, pn), solve(jacob, form_drilling_displacement(c_pt.coor, diff_coor)));
		H += c_pt.factor * poly_stress.t() * c_pt.poly_strain;
		HT += c_pt.factor * c_pt.poly_strain.t() * mat_stiff * c_pt.poly_strain;
	}

	const mat NT = solve(H, N);

	trial_stiffness = current_stiffness = initial_stiffness = NT.t() * HT * NT;

	for(auto& I : int_pt) I.poly_strain *= NT;

	initial_mass.zeros(m_size, m_size);
	const auto t_density = mat_proto->get_parameter();
	if(t_density != 0.) {
		for(const auto& I : int_pt) {
			const auto n_int = shape::quad(I.coor, 0);
			const auto tmp_a = t_density * I.factor;
			for(unsigned J = 0; J < m_node; ++J) for(auto K = J; K < m_node; ++K) initial_mass(m_dof * J, m_dof * K) += tmp_a * n_int(J) * n_int(K);
		}
		for(unsigned I = 0; I < m_size; I += m_dof) {
			initial_mass(I + 1, I + 1) = initial_mass(I, I);
			for(auto J = I + m_dof; J < m_size; J += m_dof) initial_mass(J, I) = initial_mass(I + 1, J + 1) = initial_mass(J + 1, I + 1) = initial_mass(I, J);
		}
		for(const auto& I : int_pt) {
			const auto n_int = form_drilling_mass(I.coor, diff_coor);
			initial_mass += n_int.t() * n_int * t_density * I.factor;
		}
	}
	trial_mass = current_mass = initial_mass;
}

int SRGCMQ::update_status() {
	const auto trial_disp = get_trial_displacement();

	trial_resistance.zeros(m_size);
	trial_stiffness.zeros(m_size, m_size);
	for(const auto& t_pt : int_pt) {
		const vec t_strain = t_pt.poly_strain * trial_disp;
		if(t_pt.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		if(!reinforced_x && !reinforced_y) {
			trial_resistance += t_pt.poly_strain.t() * t_pt.factor * t_pt.m_material->get_trial_stress();
			trial_stiffness += t_pt.poly_strain.t() * t_pt.factor * t_pt.m_material->get_trial_stiffness() * t_pt.poly_strain;
		} else {
			auto t_stress = t_pt.m_material->get_trial_stress();
			auto t_stiffness = t_pt.m_material->get_trial_stiffness();
			if(reinforced_x) {
				t_pt.rebar_x->update_trial_status(t_strain(0));
				t_stress(0) += t_pt.rebar_x->get_trial_stress().at(0) * rho_x;
				t_stiffness(0, 0) += t_pt.rebar_x->get_trial_stiffness().at(0) * rho_x;
			}
			if(reinforced_y) {
				t_pt.rebar_y->update_trial_status(t_strain(1));
				t_stress(1) += t_pt.rebar_y->get_trial_stress().at(0) * rho_y;
				t_stiffness(1, 1) += t_pt.rebar_y->get_trial_stiffness().at(0) * rho_y;
			}
			trial_resistance += t_pt.poly_strain.t() * t_pt.factor * t_stress;
			trial_stiffness += t_pt.poly_strain.t() * t_pt.factor * t_stiffness * t_pt.poly_strain;
		}
	}

	return SUANPAN_SUCCESS;
}

int SRGCMQ::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	if(reinforced_x) for(const auto& I : int_pt) code += I.rebar_x->commit_status();
	if(reinforced_y) for(const auto& I : int_pt) code += I.rebar_y->commit_status();
	return code;
}

int SRGCMQ::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	if(reinforced_x) for(const auto& I : int_pt) code += I.rebar_x->clear_status();
	if(reinforced_y) for(const auto& I : int_pt) code += I.rebar_y->clear_status();
	return code;
}

int SRGCMQ::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	if(reinforced_x) for(const auto& I : int_pt) code += I.rebar_x->reset_status();
	if(reinforced_y) for(const auto& I : int_pt) code += I.rebar_y->reset_status();
	return code;
}

vector<vec> SRGCMQ::record(const OutputType T) {
	vector<vec> data;
	for(const auto& I : int_pt) for(const auto& J : I.m_material->record(T)) data.emplace_back(J);
	return data;
}

void SRGCMQ::print() {
	suanpan_info("SRGCMQ mixed quad element %u connects nodes:\n", get_tag());
	node_encoding.t().print();
	suanpan_info("\nMaterial model response:");
	for(size_t I = 0; I < int_pt.size(); ++I) {
		suanpan_info("\nIntegration Point %lu:\t", I + 1);
		int_pt[I].coor.t().print();
		int_pt[I].m_material->print();
	}
}
