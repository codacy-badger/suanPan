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

#include "T2D2.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

const unsigned T2D2::t_node = 2;
const unsigned T2D2::t_dof = 2;
const unsigned T2D2::t_size = t_dof * t_node;

T2D2::T2D2(const unsigned T, uvec&& N, const unsigned M, const double A, const bool F, const bool UA, const bool LS)
	: MaterialElement(T, t_node, t_dof, std::forward<uvec>(N), uvec{M}, F)
	, area(A)
	, t_trans(F ? make_unique<T2DC>() : make_unique<T2DL>())
	, update_area(UA)
	, log_strain(LS) {}

void T2D2::initialize(const shared_ptr<DomainBase>& D) {
	t_trans->set_element_ptr(this);

	access::rw(length) = t_trans->get_length();

	if(!D->find_material(unsigned(material_tag(0)))) {
		D->disable_element(get_tag());
		return;
	}

	t_material = D->get_material(unsigned(material_tag(0)))->get_copy();

	if(t_material->get_material_type() != MaterialType::D1) {
		D->disable_element(get_tag());
		return;
	}

	trial_stiffness = current_stiffness = initial_stiffness = t_trans->to_global_stiffness_mat(area / length * t_material->get_initial_stiffness());

	trial_mass = current_mass = initial_mass = t_trans->to_global_mass_mat(t_material->get_parameter() * area);
}

int T2D2::update_status() {
	auto new_area = area;

	if(nlgeom) {
		t_trans->update_status();

		const auto new_length = t_trans->get_length();

		double t_strain, d_strain;
		if(log_strain) {
			t_strain = log(new_length / length);
			d_strain = new_length;
		} else {
			t_strain = new_length / length - 1.;
			d_strain = length;
		}

		if(t_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		const auto t_stiff = update_area ? -t_material->get_trial_stress().at(0) * (new_area *= length / new_length) / new_length : 0.;

		trial_stiffness = t_trans->to_global_stiffness_mat(t_material->get_trial_stiffness() * new_area / d_strain + t_stiff);
		trial_geometry = t_trans->to_global_geometry_mat(new_area / new_length * t_material->get_trial_stress());
	} else {
		if(t_material->update_trial_status(t_trans->to_local_vec(get_trial_displacement()) / length, t_trans->to_local_vec(get_trial_velocity()) / length) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		trial_stiffness = t_trans->to_global_stiffness_mat(area / length * t_material->get_trial_stiffness());
	}

	const vec t_strain = t_trans->to_local_vec(get_trial_displacement()) / length;
	const vec t_strain_rate = t_trans->to_local_vec(get_trial_velocity()) / length;
	const auto t_stiff = t_material->get_trial_stiffness();

	trial_resistance = t_trans->to_global_vec(new_area * t_material->get_trial_stress());

	suanpan_debug([&]() { if(!trial_stiffness.is_finite() || !trial_resistance.is_finite()) throw invalid_argument("infinite number detected"); });

	return SUANPAN_SUCCESS;
}

int T2D2::commit_status() { return t_material->commit_status(); }

int T2D2::clear_status() { return t_material->clear_status(); }

int T2D2::reset_status() { return t_material->reset_status(); }

vector<vec> T2D2::record(const OutputType P) { return t_material->record(P); }

void T2D2::print() {
	suanpan_info("2-D truss element with ");
	if(nlgeom) suanpan_info("corotational formulation, assuming constant %s and %s strain. ", update_area ? "volume" : "area", log_strain ? "logarithmic" : "engineering");
	else suanpan_info("linear formulation. ");
	suanpan_info("The nodes connected are\n");
	node_encoding.t().print();
	suanpan_info("The area is %.4E. The initial element length is %.4E.\n", area, length);
	suanpan_info("Material Model: ");
	t_material->print();
}

#ifdef SUANPAN_VTK
#include <vtkLine.h>

void T2D2::Setup() {
	vtk_cell = vtkSmartPointer<vtkLine>::New();
	auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < t_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void T2D2::GetDisplacement(vtkSmartPointer<vtkDoubleArray>& arrays) {
	mat t_disp(6, t_node, fill::zeros);
	t_disp.rows(0, 1) = reshape(get_current_displacement(), t_dof, t_node);
	for(unsigned I = 0; I < t_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

void T2D2::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), t_dof, t_node).t();
	for(unsigned I = 0; I < t_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
