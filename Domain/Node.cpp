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

#include "Node.h"
#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>

Node::Node(const unsigned T)
	: Tag(T) { suanpan_debug("Node %u ctor() called.\n", T); }

Node::Node(const unsigned T, vec&& C)
	: Tag(T) {
	num_dof = unsigned(C.n_elem);
	coordinate = std::forward<vec>(C);
	suanpan_debug("Node %u ctor() called.\n", T);
}

/**
 * \brief initialize `num_dof` and set the size of `coordinate` to `num_dof`.
 * \param T `unique_tag`
 * \param D `num_dof`
 */
Node::Node(const unsigned T, const unsigned D)
	: Tag(T) {
	num_dof = D;
	coordinate.zeros(D);
	suanpan_debug("Node %u ctor() called.\n", T);
}

/**
 * \brief initialize `num_dof` and `coordinate`.
 * \param T `unique_tag`
 * \param D `num_dof`
 * \param C `coordinate`
 */
Node::Node(const unsigned T, const unsigned D, vec&& C)
	: Tag(T) {
	num_dof = D;
	coordinate = std::forward<vec>(C);
	suanpan_debug("Node %u ctor() called.\n", T);
}

/**
 * \brief default destructor.
 */
Node::~Node() { suanpan_debug("Node %u dtor() called.\n", get_tag()); }

/**
 * \brief This method should be called after Element objects are set. Element
 * objects will set the minimum number of DoFs for all related Node objects.
 * This method initialize all member variables with the size of `num_dof` and
 * fill `original_dof` with `-1` to indicated it should be omitted from the
 * system. Finally check if the size of `coordinate` is the same of `num_dof`,
 * if not, resize it to `num_dof`. This will be necessary for beam/plate/shell
 * problems which have more DoFs than coordinates.
 */
void Node::initialize(const shared_ptr<DomainBase>& D) {
	if(initialized || !is_active()) return;

	if(num_dof != 0) {
		original_dof.zeros(num_dof);
		original_dof.fill(static_cast<uword>(-1));

		reordered_dof.reset();

		current_resistance.resize(num_dof);
		current_displacement.resize(num_dof);
		current_velocity.resize(num_dof);
		current_acceleration.resize(num_dof);

		incre_resistance.resize(num_dof);
		incre_displacement.resize(num_dof);
		incre_velocity.resize(num_dof);
		incre_acceleration.resize(num_dof);

		trial_resistance.resize(num_dof);
		trial_displacement.resize(num_dof);
		trial_velocity.resize(num_dof);
		trial_acceleration.resize(num_dof);
	} else {
		suanpan_debug("Node %u is not used in the problem, now disable it.\n", get_tag());
		D->disable_node(get_tag());
	}

	initialized = true;
}

void Node::set_dof_number(const unsigned& D) {
	if(num_dof == D) return;

	num_dof = D;
	initialized = false;
}

const unsigned& Node::get_dof_number() const { return num_dof; }

void Node::set_original_dof(unsigned& F) {
	if(!is_active()) return;

	for(unsigned I = 0; I < num_dof; ++I, ++F)
		if(original_dof(I) != F) {
			original_dof(I) = F;
			initialized = false;
		}
}

void Node::set_original_dof(const uvec& D) {
	if(original_dof.size() != D.size()) {
		original_dof = D;
		initialized = false;
	} else
		for(uword I = 0; I < D.size(); ++I)
			if(original_dof(I) != D(I)) {
				original_dof = D;
				initialized = false;
				break;
			}
}

const uvec& Node::get_original_dof() const { return original_dof; }

void Node::set_reordered_dof(const uvec& R) { reordered_dof = R; }

const uvec& Node::get_reordered_dof() const { return reordered_dof.is_empty() ? original_dof : reordered_dof; }

void Node::set_coordinate(const vec& C) { coordinate = C; }

const vec& Node::get_coordinate() const { return coordinate; }

void Node::set_current_resistance(const vec& R) { current_resistance = R; }

void Node::set_current_displacement(const vec& D) { current_displacement = D; }

void Node::set_current_velocity(const vec& V) { current_velocity = V; }

void Node::set_current_acceleration(const vec& A) { current_acceleration = A; }

void Node::set_incre_resistance(const vec& R) { incre_resistance = R; }

void Node::set_incre_displacement(const vec& D) { incre_displacement = D; }

void Node::set_incre_velocity(const vec& V) { incre_velocity = V; }

void Node::set_incre_acceleration(const vec& A) { incre_acceleration = A; }

void Node::set_trial_resistance(const vec& R) { trial_resistance = R; }

void Node::set_trial_displacement(const vec& D) { trial_displacement = D; }

void Node::set_trial_velocity(const vec& V) { trial_velocity = V; }

void Node::set_trial_acceleration(const vec& A) { trial_acceleration = A; }

const vec& Node::get_current_resistance() const { return current_resistance; }

const vec& Node::get_current_displacement() const { return current_displacement; }

const vec& Node::get_current_velocity() const { return current_velocity; }

const vec& Node::get_current_acceleration() const { return current_acceleration; }

const vec& Node::get_incre_resistance() const { return incre_resistance; }

const vec& Node::get_incre_displacement() const { return incre_displacement; }

const vec& Node::get_incre_velocity() const { return incre_velocity; }

const vec& Node::get_incre_acceleration() const { return incre_acceleration; }

const vec& Node::get_trial_resistance() const { return trial_resistance; }

const vec& Node::get_trial_displacement() const { return trial_displacement; }

const vec& Node::get_trial_velocity() const { return trial_velocity; }

const vec& Node::get_trial_acceleration() const { return trial_acceleration; }

void Node::update_current_resistance(const vec& R) {
	trial_resistance = current_resistance = R(reordered_dof);
	incre_resistance.zeros();
}

void Node::update_incre_resistance(const vec& R) {
	incre_resistance = R(reordered_dof);
	trial_resistance = current_resistance + incre_resistance;
}

void Node::update_trial_resistance(const vec& R) {
	trial_resistance = R(reordered_dof);
	incre_resistance = trial_resistance - current_resistance;
}

void Node::update_current_status(const vec& D) {
	trial_displacement = current_displacement = D(reordered_dof);
	incre_displacement.zeros();
}

void Node::update_current_status(const vec& D, const vec& V) {
	trial_velocity = current_velocity = V(reordered_dof);
	incre_velocity.zeros();
	update_current_status(D);
}

void Node::update_current_status(const vec& D, const vec& V, const vec& A) {
	trial_acceleration = current_acceleration = A(reordered_dof);
	incre_acceleration.zeros();
	update_current_status(D, V);
}

void Node::update_incre_status(const vec& D) {
	incre_displacement = D(reordered_dof);
	trial_displacement = current_displacement + incre_displacement;
}

void Node::update_incre_status(const vec& D, const vec& V) {
	incre_velocity = V(reordered_dof);
	trial_velocity = current_velocity + incre_velocity;
	update_incre_status(D);
}

void Node::update_incre_status(const vec& D, const vec& V, const vec& A) {
	incre_acceleration = A(reordered_dof);
	trial_acceleration = current_acceleration + incre_acceleration;
	update_incre_status(D, V);
}

void Node::update_trial_status(const vec& D) {
	trial_displacement = D(reordered_dof);
	incre_displacement = trial_displacement - current_displacement;
}

void Node::update_trial_status(const vec& D, const vec& V) {
	trial_velocity = V(reordered_dof);
	incre_velocity = trial_velocity - current_velocity;
	update_trial_status(D);
}

void Node::update_trial_status(const vec& D, const vec& V, const vec& A) {
	trial_acceleration = A(reordered_dof);
	incre_acceleration = trial_acceleration - current_acceleration;
	update_trial_status(D, V);
}

void Node::commit_status() {
	if(!trial_resistance.is_empty()) {
		current_resistance = trial_resistance;
		incre_resistance.zeros();
	}
	if(!trial_displacement.is_empty()) {
		current_displacement = trial_displacement;
		incre_displacement.zeros();
	}
	if(!trial_velocity.is_empty()) {
		current_velocity = trial_velocity;
		incre_displacement.zeros();
	}
	if(!trial_acceleration.is_empty()) {
		current_acceleration = trial_acceleration;
		incre_acceleration.zeros();
	}
}

void Node::reset_status() {
	if(!current_resistance.is_empty()) {
		trial_resistance = current_resistance;
		incre_resistance.zeros();
	}
	if(!current_displacement.is_empty()) {
		trial_displacement = current_displacement;
		incre_displacement.zeros();
	}
	if(!current_velocity.is_empty()) {
		trial_velocity = current_velocity;
		incre_velocity.zeros();
	}
	if(!current_acceleration.is_empty()) {
		trial_acceleration = current_acceleration;
		incre_acceleration.zeros();
	}
}

void Node::clear_status() {
	if(!current_resistance.is_empty()) {
		current_resistance.zeros();
		incre_resistance.zeros();
		trial_resistance.zeros();
	}
	if(!current_displacement.is_empty()) {
		current_displacement.zeros();
		incre_displacement.zeros();
		trial_displacement.zeros();
	}
	if(!current_velocity.is_empty()) {
		current_velocity.zeros();
		incre_velocity.zeros();
		trial_velocity.zeros();
	}
	if(!current_acceleration.is_empty()) {
		current_acceleration.zeros();
		incre_acceleration.zeros();
		trial_acceleration.zeros();
	}
}

vector<vec> Node::record(const OutputType L) const {
	vector<vec> data;

	if(L == OutputType::RF) data.push_back(current_resistance);
	else if(L == OutputType::U) data.push_back(current_displacement);
	else if(L == OutputType::V) data.push_back(current_velocity);
	else if(L == OutputType::A) data.push_back(current_acceleration);
	else if(L == OutputType::U1 && current_displacement.n_elem >= 1) data.emplace_back(vec{current_displacement(0)});
	else if(L == OutputType::U2 && current_displacement.n_elem >= 2) data.emplace_back(vec{current_displacement(1)});
	else if(L == OutputType::U3 && current_displacement.n_elem >= 3) data.emplace_back(vec{current_displacement(2)});
	else if(L == OutputType::UR1 && current_displacement.n_elem >= 4) data.emplace_back(vec{current_displacement(3)});
	else if(L == OutputType::UR2 && current_displacement.n_elem >= 5) data.emplace_back(vec{current_displacement(4)});
	else if(L == OutputType::UR3 && current_displacement.n_elem >= 6) data.emplace_back(vec{current_displacement(5)});
	else if(L == OutputType::RF1 && current_resistance.n_elem >= 1) data.emplace_back(vec{current_resistance(0)});
	else if(L == OutputType::RF2 && current_resistance.n_elem >= 2) data.emplace_back(vec{current_resistance(1)});
	else if(L == OutputType::RF3 && current_resistance.n_elem >= 3) data.emplace_back(vec{current_resistance(2)});
	else if(L == OutputType::RM1 && current_resistance.n_elem >= 4) data.emplace_back(vec{current_resistance(3)});
	else if(L == OutputType::RM2 && current_resistance.n_elem >= 5) data.emplace_back(vec{current_resistance(4)});
	else if(L == OutputType::RM3 && current_resistance.n_elem >= 6) data.emplace_back(vec{current_resistance(5)});
	else if(L == OutputType::V1 && current_velocity.n_elem >= 1) data.emplace_back(vec{current_velocity(0)});
	else if(L == OutputType::V2 && current_velocity.n_elem >= 2) data.emplace_back(vec{current_velocity(1)});
	else if(L == OutputType::V3 && current_velocity.n_elem >= 3) data.emplace_back(vec{current_velocity(2)});
	else if(L == OutputType::VR1 && current_velocity.n_elem >= 4) data.emplace_back(vec{current_velocity(3)});
	else if(L == OutputType::VR2 && current_velocity.n_elem >= 5) data.emplace_back(vec{current_velocity(4)});
	else if(L == OutputType::VR3 && current_velocity.n_elem >= 6) data.emplace_back(vec{current_velocity(5)});
	else if(L == OutputType::A1 && current_acceleration.n_elem >= 1) data.emplace_back(vec{current_acceleration(0)});
	else if(L == OutputType::A2 && current_acceleration.n_elem >= 2) data.emplace_back(vec{current_acceleration(1)});
	else if(L == OutputType::A3 && current_acceleration.n_elem >= 3) data.emplace_back(vec{current_acceleration(2)});
	else if(L == OutputType::AR1 && current_acceleration.n_elem >= 4) data.emplace_back(vec{current_acceleration(3)});
	else if(L == OutputType::AR2 && current_acceleration.n_elem >= 5) data.emplace_back(vec{current_acceleration(4)});
	else if(L == OutputType::AR3 && current_acceleration.n_elem >= 6) data.emplace_back(vec{current_acceleration(5)});

	return data;
}

void Node::print() {
	suanpan_info("Node %u:\n", get_tag(), is_active() ? "" : " is currently inactive");
	coordinate.t().print();
	suanpan_info("Displacement:\n");
	current_displacement.t().print();
	suanpan_info("Resistance:\n");
	current_resistance.t().print();
	if(accu(current_velocity) != 0.) {
		suanpan_info("Velocity:\n");
		current_velocity.t().print();
	}
	if(accu(current_acceleration) != 0.) {
		suanpan_info("Acceleration:\n");
		current_acceleration.t().print();
	}
}
