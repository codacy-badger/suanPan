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

#include "T2D2S.h"
#include <Domain/DomainBase.h>
#include <Section/Section.h>

const unsigned T2D2S::t_node = 2;
const unsigned T2D2S::t_dof = 2;
const unsigned T2D2S::t_size = t_dof * t_node;

T2D2S::T2D2S(const unsigned T, uvec&& N, const unsigned M, const bool F, const bool LS)
	: SectionElement(T, t_node, t_dof, std::forward<uvec>(N), uvec{M}, F)
	, t_trans(F ? make_unique<T2DC>() : make_unique<T2DL>())
	, log_strain(LS) {}

void T2D2S::initialize(const shared_ptr<DomainBase>& D) {
	t_trans->set_element_ptr(this);

	access::rw(length) = t_trans->get_length();

	t_section = D->get_section(unsigned(section_tag(0)))->get_copy();

	trial_stiffness = current_stiffness = initial_stiffness = t_trans->to_global_stiffness_mat(t_section->get_initial_stiffness() / length);

	trial_mass = current_mass = initial_mass = t_trans->to_global_mass_mat(t_section->get_parameter(ParameterType::LINEARDENSITY));
}

int T2D2S::update_status() {
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

		if(t_section->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		trial_stiffness = t_trans->to_global_stiffness_mat(t_section->get_trial_stiffness() / d_strain);
		trial_geometry = t_trans->to_global_geometry_mat(t_section->get_trial_resistance() / new_length);
	} else {
		if(t_section->update_trial_status(t_trans->to_local_vec(get_trial_displacement()) / length, t_trans->to_local_vec(get_trial_velocity()) / length) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		trial_stiffness = t_trans->to_global_stiffness_mat(t_section->get_trial_stiffness() / length);
	}

	trial_resistance = t_trans->to_global_vec(t_section->get_trial_resistance());

	return SUANPAN_SUCCESS;
}

int T2D2S::commit_status() { return t_section->commit_status(); }

int T2D2S::clear_status() { return t_section->clear_status(); }

int T2D2S::reset_status() { return t_section->reset_status(); }

vector<vec> T2D2S::record(const OutputType P) { return t_section->record(P); }

void T2D2S::print() {
	suanpan_info("2-D truss element with ");
	if(nlgeom) suanpan_info("corotational formulation, assuming constant area and %s strain. ", log_strain ? "logarithmic" : "engineering");
	else suanpan_info("linear formulation. ");
	suanpan_info("The nodes connected are\n");
	node_encoding.t().print();
	suanpan_info("Section model: ");
	t_section->print();
}
