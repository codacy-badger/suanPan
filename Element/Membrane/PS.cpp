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

#include "PS.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/utility.h>

const unsigned PS::m_node = 4;
const unsigned PS::m_dof = 2;
const unsigned PS::m_size = m_dof * m_node;

PS::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
	: coor(std::forward<vec>(C))
	, weight(W)
	, m_material(std::forward<unique_ptr<Material>>(M)) {}

mat PS::form_transformation(const mat& jacobian) {
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

	return trans_mat;
}

PS::PS(const unsigned T, uvec&& N, const unsigned M, const double TH)
	: MaterialElement(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, false)
	, thickness(TH) {}

void PS::initialize(const shared_ptr<DomainBase>& D) {
	const auto material_proto = std::dynamic_pointer_cast<Material2D>(D->get_material(unsigned(material_tag(0))));

	if(material_proto == nullptr) {
		D->disable_element(get_tag());
		return;
	}

	if(material_proto->plane_type == PlaneType::E) suanpan::hacker(thickness) = 1.;

	auto& ini_stiffness = material_proto->get_initial_stiffness();

	const auto ele_coor = get_coordinate(2);

	const auto jacob_trans = form_transformation(shape::quad(vec{0., 0.}, 1) * ele_coor);

	const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

	mat poly_disp(3, 8, fill::zeros), poly_stress(3, 5, fill::zeros);
	for(auto J = 0; J < 3; ++J) poly_stress(J, J) = 1.;

	mat H(5, 5, fill::zeros), N(5, 8, fill::zeros);
	int_pt.clear(), int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		vec t_vec{plan(I, 0), plan(I, 1)};
		const auto pn = shape::quad(t_vec, 1);
		const mat jacob = pn * ele_coor;
		int_pt.emplace_back(std::move(t_vec), thickness * plan(I, 2) * det(jacob), material_proto->get_copy());

		auto& c_pt = int_pt.back();

		const mat pn_pxy = solve(jacob, pn);
		for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
			poly_disp(2, L) = poly_disp(0, K) = pn_pxy(0, J);
			poly_disp(2, K) = poly_disp(1, L) = pn_pxy(1, J);
		}

		poly_stress.col(3) = jacob_trans.col(0) * c_pt.coor(1);
		poly_stress.col(4) = jacob_trans.col(1) * c_pt.coor(0);

		c_pt.poly_strain = solve(ini_stiffness, poly_stress);

		N += c_pt.weight * poly_stress.t() * poly_disp;
		H += c_pt.weight * poly_stress.t() * c_pt.poly_strain;
	}

	const mat NT = solve(H, N);

	trial_stiffness = current_stiffness = initial_stiffness = N.t() * NT;

	for(auto& I : int_pt) I.poly_strain *= NT;

	trial_mass.zeros(m_size, m_size);
	const auto tmp_density = material_proto->get_parameter(ParameterType::DENSITY);
	if(tmp_density != 0.) {
		for(const auto& I : int_pt) {
			const auto n_int = shape::quad(I.coor, 0);
			const auto tmp_a = tmp_density * I.weight;
			for(unsigned J = 0; J < m_node; ++J) for(auto K = J; K < m_node; ++K) trial_mass(m_dof * J, m_dof * K) += tmp_a * n_int(J) * n_int(K);
		}
		for(unsigned I = 0; I < m_size; I += m_dof) {
			trial_mass(I + 1, I + 1) = trial_mass(I, I);
			for(auto J = I + m_dof; J < m_size; J += m_dof) trial_mass(J, I) = trial_mass(I + 1, J + 1) = trial_mass(J + 1, I + 1) = trial_mass(I, J);
		}
	}
	initial_mass = current_mass = trial_mass;
}

int PS::update_status() {
	const auto trial_disp = get_trial_displacement();

	trial_resistance.zeros(m_size);
	trial_stiffness.zeros(m_size, m_size);
	for(const auto& I : int_pt) {
		if(I.m_material->update_trial_status(I.poly_strain * trial_disp) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		trial_resistance += I.weight * I.poly_strain.t() * I.m_material->get_trial_stress();
		trial_stiffness += I.weight * I.poly_strain.t() * I.m_material->get_trial_stiffness() * I.poly_strain;
	}

	return SUANPAN_SUCCESS;
}

int PS::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	return code;
}

int PS::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	return code;
}

int PS::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	return code;
}

vector<vec> PS::record(const OutputType P) {
	vector<vec> output;
	output.reserve(int_pt.size());

	for(const auto& I : int_pt) for(const auto& J : I.m_material->record(P)) output.emplace_back(J);

	return output;
}

void PS::print() {
	suanpan_info("Element %u is a four-node membrane element (Pian-Sumihara).\n", get_tag());
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

void PS::Setup() {
	vtk_cell = vtkSmartPointer<vtkQuad>::New();
	auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < m_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void PS::GetDisplacement(vtkSmartPointer<vtkDoubleArray>& arrays) {
	mat t_disp(6, m_node, fill::zeros);
	t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);
	for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

void PS::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
	for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
