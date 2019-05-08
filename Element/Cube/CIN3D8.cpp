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

#include "CIN3D8.h"
#include <Domain/DomainBase.h>
#include <Material/Material3D/Material3D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/tensorToolbox.h>

const unsigned CIN3D8::c_node = 8;
const unsigned CIN3D8::c_dof = 3;
const unsigned CIN3D8::c_size = c_dof * c_node;

CIN3D8::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& P)
	: coor(std::forward<vec>(C))
	, weight(W)
	, c_material(std::forward<unique_ptr<Material>>(M))
	, pn_pxyz(std::forward<mat>(P))
	, strain_mat(6, c_size, fill::zeros) {}

mat CIN3D8::compute_mapping(const vec& C) {
	const auto& X = C(0);
	const auto& Y = C(1);
	const auto& Z = C(2);

	mat m(c_node, 3);

	const auto XP = X + 1.;
	const auto XM = X - 1.;
	const auto YP = Y + 1.;
	const auto YM = Y - 1.;
	const auto ZP = Z + 1.;
	const auto ZM = Z - 1.;
	const auto ZM2 = 2. * ZM;
	const auto ZMZM = ZM2 * ZM;
	const auto ZM4 = 2. * ZM2;

	m(0) = Z * YM / ZM2;
	m(1) = -Z * YM / ZM2;
	m(2) = Z * YP / ZM2;
	m(3) = -Z * YP / ZM2;
	m(4) = -YM * ZP / ZM4;
	m(5) = YM * ZP / ZM4;
	m(6) = -YP * ZP / ZM4;
	m(7) = YP * ZP / ZM4;
	m(8) = Z * XM / ZM2;
	m(9) = -Z * XP / ZM2;
	m(10) = Z * XP / ZM2;
	m(11) = -Z * XM / ZM2;
	m(12) = -XM * ZP / ZM4;
	m(13) = XP * ZP / ZM4;
	m(14) = -XP * ZP / ZM4;
	m(15) = XM * ZP / ZM4;
	m(16) = -XM * YM / ZMZM;
	m(17) = XP * YM / ZMZM;
	m(18) = -XP * YP / ZMZM;
	m(19) = XM * YP / ZMZM;
	m(20) = XM * YM / ZMZM;
	m(21) = -XP * YM / ZMZM;
	m(22) = XP * YP / ZMZM;
	m(23) = -XM * YP / ZMZM;

	return m.t();
}

mat CIN3D8::compute_n(const vec& C) {
	const auto& X = C(0);
	const auto& Y = C(1);
	const auto& Z = C(2);

	mat n(c_node, 1);

	const auto XP = 1. + X;
	const auto XM = 1. - X;
	const auto YP = 1. + Y;
	const auto YM = 1. - Y;
	const auto ZZMZ = Z * Z - Z;
	const auto ZZM = 1. - Z * Z;

	n(0) = XM * YM * ZZMZ * .125;
	n(1) = XP * YM * ZZMZ * .125;
	n(2) = XP * YP * ZZMZ * .125;
	n(3) = XM * YP * ZZMZ * .125;
	n(4) = XM * YM * ZZM * .25;
	n(5) = XP * YM * ZZM * .25;
	n(6) = XP * YP * ZZM * .25;
	n(7) = XM * YP * ZZM * .25;

	return n;
}

mat CIN3D8::compute_dn(const vec& C) {
	const auto& X = C(0);
	const auto& Y = C(1);
	const auto& Z = C(2);

	mat pn(c_node, 3);

	const auto XP = X + 1.;
	const auto XM = X - 1.;
	const auto YP = Y + 1.;
	const auto YM = Y - 1.;
	const auto ZP = Z + 1.;
	const auto ZM = Z - 1.;
	const auto ZM2 = 2. * Z - 1.;

	pn(0) = Z * YM * ZM * .125;
	pn(1) = -Z * YM * ZM * .125;
	pn(2) = Z * YP * ZM * .125;
	pn(3) = -Z * YP * ZM * .125;
	pn(4) = -ZP * ZM * YM * .25;
	pn(5) = ZP * ZM * YM * .25;
	pn(6) = -ZP * ZM * YP * .25;
	pn(7) = ZP * ZM * YP * .25;
	pn(8) = Z * XM * ZM * .125;
	pn(9) = -Z * XP * ZM * .125;
	pn(10) = Z * XP * ZM * .125;
	pn(11) = -Z * XM * ZM * .125;
	pn(12) = -ZP * ZM * XM * .25;
	pn(13) = ZP * ZM * XP * .25;
	pn(14) = -ZP * ZM * XP * .25;
	pn(15) = ZP * ZM * XM * .25;
	pn(16) = ZM2 * XM * YM * .125;
	pn(17) = -ZM2 * XP * YM * .125;
	pn(18) = ZM2 * XP * YP * .125;
	pn(19) = -ZM2 * XM * YP * .125;
	pn(20) = -Z * XM * YM * .5;
	pn(21) = Z * XP * YM * .5;
	pn(22) = -Z * XP * YP * .5;
	pn(23) = Z * XM * YP * .5;

	return pn.t();
}

CIN3D8::CIN3D8(const unsigned T, uvec&& N, const unsigned M)
	: MaterialElement(T, c_node, c_dof, std::forward<uvec>(N), uvec{M}, false) {}

void CIN3D8::initialize(const shared_ptr<DomainBase>& D) {
	const auto ele_coor = get_coordinate(c_dof);

	auto& material_proto = D->get_material(unsigned(material_tag(0)));

	auto& ini_stiffness = material_proto->get_initial_stiffness();

	const IntegrationPlan plan(3, 2, IntegrationType::GAUSS);

	initial_stiffness.zeros(c_size, c_size);

	int_pt.clear(), int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		vec t_vec{plan(I, 0), plan(I, 1), plan(I, 2)};
		const mat jacob = compute_mapping(t_vec) * ele_coor;
		int_pt.emplace_back(std::move(t_vec), plan(I, c_dof) * det(jacob), material_proto->get_copy(), solve(jacob, compute_dn(t_vec)));

		auto& c_pt = int_pt.back();
		for(unsigned J = 0, K = 0, L = 1, M = 2; J < c_node; ++J, K += c_dof, L += c_dof, M += c_dof) {
			c_pt.strain_mat(0, K) = c_pt.strain_mat(3, L) = c_pt.strain_mat(5, M) = c_pt.pn_pxyz(0, J);
			c_pt.strain_mat(3, K) = c_pt.strain_mat(1, L) = c_pt.strain_mat(4, M) = c_pt.pn_pxyz(1, J);
			c_pt.strain_mat(5, K) = c_pt.strain_mat(4, L) = c_pt.strain_mat(2, M) = c_pt.pn_pxyz(2, J);
		}
		initial_stiffness += c_pt.weight * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;
	}
	trial_stiffness = current_stiffness = initial_stiffness;

	initial_mass.zeros(c_size, c_size);
	const auto t_density = material_proto->get_parameter();
	if(t_density != 0.) {
		for(const auto& I : int_pt) {
			const auto n_int = compute_n(I.coor);
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

int CIN3D8::update_status() {
	const auto t_disp = get_trial_displacement();

	trial_stiffness.zeros(c_size, c_size);
	trial_resistance.zeros(c_size);

	for(const auto& I : int_pt) {
		if(I.c_material->update_trial_status(I.strain_mat * t_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		trial_stiffness += I.weight * I.strain_mat.t() * I.c_material->get_trial_stiffness() * I.strain_mat;
		trial_resistance += I.weight * I.strain_mat.t() * I.c_material->get_trial_stress();
	}

	return SUANPAN_SUCCESS;
}

int CIN3D8::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.c_material->commit_status();
	return code;
}

int CIN3D8::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.c_material->clear_status();
	return code;
}

int CIN3D8::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.c_material->reset_status();
	return code;
}

vector<vec> CIN3D8::record(const OutputType T) {
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

void CIN3D8::print() {
	suanpan_info("CIN3D8 element.\n");
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

void CIN3D8::Setup() {
	vtk_cell = vtkSmartPointer<vtkHexahedron>::New();
	auto ele_coor = get_coordinate(3);
	for(unsigned I = 0; I < c_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), ele_coor(I, 2));
	}
}

void CIN3D8::GetDisplacement(vtkSmartPointer<vtkDoubleArray>& arrays) {
	mat t_disp(6, c_node, fill::zeros);
	t_disp.rows(0, 2) = reshape(get_current_displacement(), c_dof, c_node);
	for(unsigned I = 0; I < c_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

void CIN3D8::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(3) + amplifier * reshape(get_current_displacement(), c_dof, c_node).t();
	for(unsigned I = 0; I < c_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), ele_disp(I, 2));
}

#endif