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

#include "Element.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material.h>
#include <Section/Section.h>

/**
 * \brief generate a matrix that contains coordinates of connected nodes
 * \param num_dim number of dimension required
 * \return a matrix of following form
 *         | x_1  y_1  z_1  ... |
 *         | x_2  y_2  z_2  ... |
 *         | x_3  y_3  z_3  ... |
 */
mat Element::get_coordinate(const unsigned num_dim) const {
	mat ele_coor(num_node, num_dim, fill::zeros);

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_coor = node_ptr[I].lock()->get_coordinate();
		for(uword J = 0; J < std::min(uword(num_dim), t_coor.n_elem); ++J) ele_coor(I, J) = t_coor(J);
	}

	return ele_coor;
}

vec Element::get_incre_displacement() const {
	vec incre_displacement(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_disp = node_ptr[I].lock()->get_incre_displacement();
		for(unsigned J = 0; J < num_dof; ++J) incre_displacement(idx++) = t_disp(J);
	}

	return incre_displacement;
}

vec Element::get_incre_velocity() const {
	vec incre_velocity(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_vec = node_ptr[I].lock()->get_incre_velocity();
		for(unsigned J = 0; J < num_dof; ++J) incre_velocity(idx++) = t_vec(J);
	}

	return incre_velocity;
}

vec Element::get_incre_acceleration() const {
	vec incre_acceleration(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_acc = node_ptr[I].lock()->get_incre_acceleration();
		for(unsigned J = 0; J < num_dof; ++J) incre_acceleration(idx++) = t_acc(J);
	}

	return incre_acceleration;
}

vec Element::get_trial_displacement() const {
	vec trial_displacement(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_disp = node_ptr[I].lock()->get_trial_displacement();
		for(unsigned J = 0; J < num_dof; ++J) trial_displacement(idx++) = t_disp(J);
	}

	return trial_displacement;
}

vec Element::get_trial_velocity() const {
	vec trial_velocity(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_vec = node_ptr[I].lock()->get_trial_velocity();
		for(unsigned J = 0; J < num_dof; ++J) trial_velocity(idx++) = t_vec(J);
	}

	return trial_velocity;
}

vec Element::get_trial_acceleration() const {
	vec trial_acceleration(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_acc = node_ptr[I].lock()->get_trial_acceleration();
		for(unsigned J = 0; J < num_dof; ++J) trial_acceleration(idx++) = t_acc(J);
	}

	return trial_acceleration;
}

vec Element::get_current_displacement() const {
	vec current_displacement(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_disp = node_ptr[I].lock()->get_current_displacement();
		for(unsigned J = 0; J < num_dof; ++J) current_displacement(idx++) = t_disp(J);
	}

	return current_displacement;
}

vec Element::get_current_velocity() const {
	vec current_velocity(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_vec = node_ptr[I].lock()->get_current_velocity();
		for(unsigned J = 0; J < num_dof; ++J) current_velocity(idx++) = t_vec(J);
	}

	return current_velocity;
}

vec Element::get_current_acceleration() const {
	vec current_acceleration(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_acc = node_ptr[I].lock()->get_current_acceleration();
		for(unsigned J = 0; J < num_dof; ++J) current_acceleration(idx++) = t_acc(J);
	}

	return current_acceleration;
}

vec Element::get_node_incre_resistance() const {
	vec node_incre_resistance(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_acc = node_ptr[I].lock()->get_incre_resistance();
		for(unsigned J = 0; J < num_dof; ++J) node_incre_resistance(idx++) = t_acc(J);
	}

	return node_incre_resistance;
}

vec Element::get_node_trial_resistance() const {
	vec node_trial_resistance(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_acc = node_ptr[I].lock()->get_trial_resistance();
		for(unsigned J = 0; J < num_dof; ++J) node_trial_resistance(idx++) = t_acc(J);
	}

	return node_trial_resistance;
}

vec Element::get_node_current_resistance() const {
	vec node_current_resistance(num_size);

	auto idx = 0;

	for(unsigned I = 0; I < num_node; ++I) {
		auto& t_acc = node_ptr[I].lock()->get_current_resistance();
		for(unsigned J = 0; J < num_dof; ++J) node_current_resistance(idx++) = t_acc(J);
	}

	return node_current_resistance;
}

vector<shared_ptr<Material>> Element::get_material(const shared_ptr<DomainBase>& D) const {
	vector<shared_ptr<Material>> material_pool;
	for(const auto& I : material_tag) material_pool.emplace_back(D->find_material(unsigned(I)) ? D->get_material(unsigned(I)) : nullptr);
	return material_pool;
}

vector<shared_ptr<Section>> Element::get_section(const shared_ptr<DomainBase>& D) const {
	vector<shared_ptr<Section>> section_pool;
	for(const auto& I : section_tag) section_pool.emplace_back(D->find_section(unsigned(I)) ? D->get_section(unsigned(I)) : nullptr);
	return section_pool;
}

/**
 * \brief the basic ctor of a typical element
 * \param T element tag
 * \param NN  number of nodes
 * \param ND  number of dofs
 * \param NT node tags
 * \param MT material tags
 * \param ST section tags
 * \param F nlgeom switch
 */
Element::Element(const unsigned T, const unsigned NN, const unsigned ND, uvec&& NT, uvec&& MT, uvec&& ST, const bool F)
	: ElementData{std::forward<uvec>(NT), std::forward<uvec>(MT), std::forward<uvec>(ST), F, true, true, true, true, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}
	, Tag(T)
	, num_node(NN)
	, num_dof(ND) { suanpan_debug("Element %u ctor() called.\n", T); }

Element::~Element() { suanpan_debug("Element %u dtor() called.\n", get_tag()); }

void Element::initialize(const shared_ptr<DomainBase>& D) {
	// initialized before check node vadality
	if(node_ptr.size() == num_node) {
		for(const auto& I : node_ptr) {
			const auto& t_node = I.lock();
			if(t_node == nullptr || !t_node->is_active()) {
				D->disable_element(get_tag());
				return;
			}
		}
		return;
	}

	// first initiliazation
	const auto total_dof = num_node * num_dof;

	if(total_dof == 0) {
		D->disable_element(get_tag());
		return;
	}

	dof_encoding.set_size(total_dof);

	// check if nodes are still valid
	node_ptr.clear(), node_ptr.reserve(num_node);
	for(const auto& t_tag : node_encoding) {
		if(!D->find_node(unsigned(t_tag)) || !D->get_node(unsigned(t_tag))->is_active()) {
			suanpan_debug("Element %u finds an invalid node %u, now disable it.\n", get_tag(), t_tag);
			D->disable_element(get_tag());
			return;
		}
		auto& t_node = D->get_node(unsigned(t_tag));
		if(t_node->get_dof_number() < num_dof) t_node->set_dof_number(num_dof);
		node_ptr.emplace_back(t_node);
	}

	// check if material models are valid
	for(const auto& t_material : material_tag) {
		const auto t_tag = unsigned(t_material);
		if(!D->find_material(t_tag) || !D->get_material(t_tag)->is_active()) {
			suanpan_debug("Element %u cannot find valid material %u, now disable it.\n", get_tag(), t_tag);
			D->disable_element(get_tag());
			return;
		}
	}

	// check if section models are valid
	for(const auto& t_section : section_tag) {
		const auto t_tag = unsigned(t_section);
		if(!D->find_section(t_tag) || !D->get_section(t_tag)->is_active()) {
			suanpan_debug("Element %u cannot find valid section %u, now disable it.\n", get_tag(), t_tag);
			D->disable_element(get_tag());
			return;
		}
	}

#ifdef SUANPAN_VTK
	// vtk visualization setup
	Setup();
#endif
}

void Element::set_initialized(const bool F) const { access::rw(initialized) = F; }

void Element::set_symmetric(const bool F) const { access::rw(symmetric) = F; }

bool Element::is_initialized() const { return initialized; }

bool Element::is_symmetric() const { return symmetric; }

bool Element::is_nlgeom() const { return nlgeom; }

void Element::update_dof_encoding() {
	unsigned idx = 0;
	for(const auto& tmp_ptr : node_ptr) {
		auto& node_dof = tmp_ptr.lock()->get_reordered_dof();
		for(unsigned i = 0; i < num_dof; ++i) dof_encoding(idx++) = node_dof(i);
	}
}

const uvec& Element::get_dof_encoding() const { return dof_encoding; }

const uvec& Element::get_node_encoding() const { return node_encoding; }

unsigned Element::get_dof_number() const { return num_dof; }

unsigned Element::get_node_number() const { return num_node; }

unsigned Element::get_total_number() const { return num_size; }

const vector<weak_ptr<Node>>& Element::get_node_ptr() const { return node_ptr; }

const vec& Element::get_trial_resistance() const { return trial_resistance; }

const vec& Element::get_current_resistance() const { return current_resistance; }

const mat& Element::get_trial_mass() const { return trial_mass; }

const mat& Element::get_trial_damping() const { return trial_damping; }

const mat& Element::get_trial_stiffness() const { return trial_stiffness; }

const mat& Element::get_trial_geometry() const { return trial_geometry; }

const mat& Element::get_trial_secant() const { return get_trial_stiffness(); }

const mat& Element::get_current_mass() const { return current_mass; }

const mat& Element::get_current_damping() const { return current_damping; }

const mat& Element::get_current_stiffness() const { return current_stiffness; }

const mat& Element::get_current_geometry() const { return current_geometry; }

const mat& Element::get_current_secant() const { return get_current_stiffness(); }

const mat& Element::get_initial_mass() const { return initial_mass; }

const mat& Element::get_initial_damping() const { return initial_damping; }

const mat& Element::get_initial_stiffness() const { return initial_stiffness; }

const mat& Element::get_initial_geometry() const { return initial_geometry; }

const mat& Element::get_initial_secant() const { return get_initial_stiffness(); }

int Element::update_status() { throw invalid_argument("hidden method called.\n"); }

int Element::clear_status() {
	if(update_mass && !initial_mass.is_empty()) trial_mass = current_mass = initial_mass;
	if(update_damping && !initial_damping.is_empty()) trial_damping = current_damping = initial_damping;
	if(update_stiffness && !initial_stiffness.is_empty()) trial_stiffness = current_stiffness = initial_stiffness;
	if(update_geometry && !initial_geometry.is_empty()) trial_geometry = current_geometry = initial_geometry;

	if(!trial_resistance.is_empty()) trial_resistance.zeros();
	if(!current_resistance.is_empty()) current_resistance.zeros();

	return 0;
}

int Element::commit_status() {
	if(update_mass && !trial_mass.is_empty()) current_mass = trial_mass;
	if(update_damping && !trial_damping.is_empty()) current_damping = trial_damping;
	if(update_stiffness && !trial_stiffness.is_empty()) current_stiffness = trial_stiffness;
	if(update_geometry && !trial_geometry.is_empty()) current_geometry = trial_geometry;
	if(!trial_resistance.is_empty()) current_resistance = trial_resistance;

	return 0;
}

int Element::reset_status() {
	if(update_mass && !trial_mass.is_empty()) trial_mass = current_mass;
	if(update_damping && !trial_damping.is_empty()) trial_damping = current_damping;
	if(update_stiffness && !trial_stiffness.is_empty()) trial_stiffness = current_stiffness;
	if(update_geometry && !trial_geometry.is_empty()) trial_geometry = current_geometry;
	if(!trial_resistance.is_empty()) trial_resistance = current_resistance;

	return 0;
}

vector<vec> Element::record(const OutputType) { return {}; }

mat get_coordinate(const Element* const E, const unsigned N) { return E->get_coordinate(N); }

vec get_incre_displacement(const Element* const E) { return E->get_incre_displacement(); }

vec get_incre_velocity(const Element* const E) { return E->get_incre_velocity(); }

vec get_incre_acceleration(const Element* const E) { return E->get_incre_acceleration(); }

vec get_trial_displacement(const Element* const E) { return E->get_trial_displacement(); }

vec get_trial_velocity(const Element* const E) { return E->get_trial_velocity(); }

vec get_trial_acceleration(const Element* const E) { return E->get_trial_acceleration(); }

vec get_node_incre_resistance(const Element* const E) { return E->get_node_incre_resistance(); }

vec get_node_trial_resistance(const Element* const E) { return E->get_node_trial_resistance(); }
