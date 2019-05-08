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

#include "SingleSection3D.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Section/Section.h>

const unsigned SingleSection3D::s_node = 1;
const unsigned SingleSection3D::s_dof = 3;

SingleSection3D::SingleSection3D(const unsigned T, const unsigned NT, const unsigned ST)
	: SectionElement(T, s_node, s_dof, uvec{NT}, uvec{ST}, false) {}

void SingleSection3D::initialize(const shared_ptr<DomainBase>& D) {
	const auto s_tag = unsigned(section_tag(0));

	if(!D->find_section(s_tag) || D->get_section(s_tag)->get_section_type() != SectionType::D3) {
		D->disable_element(get_tag());
		return;
	}

	s_section = suanpan::make_copy(D->get_section(s_tag));

	initial_stiffness = s_section->get_initial_stiffness();
}

int SingleSection3D::update_status() {
	s_section->update_trial_status(node_ptr[0].lock()->get_trial_displacement());

	trial_stiffness = s_section->get_trial_stiffness();

	trial_resistance = s_section->get_trial_resistance();

	return SUANPAN_SUCCESS;
}

int SingleSection3D::commit_status() { return s_section->commit_status(); }

int SingleSection3D::clear_status() { return s_section->clear_status(); }

int SingleSection3D::reset_status() { return s_section->reset_status(); }

void SingleSection3D::print() {
	suanpan_info("A SingleSection3D element that represents a section which can be used for section analysis.\n");
	s_section->print();
}
