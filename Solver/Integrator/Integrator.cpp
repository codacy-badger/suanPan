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

#include "Integrator.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

Integrator::Integrator(const unsigned T)
	: Tag(T) { suanpan_debug("Integrator %u ctor() called.\n", T); }

Integrator::~Integrator() { suanpan_debug("Integrator %u dtor() called.\n", get_tag()); }

void Integrator::set_domain(const weak_ptr<DomainBase>& D) { if(database.lock() != D.lock()) database = D; }

const weak_ptr<DomainBase>& Integrator::get_domain() const { return database; }

int Integrator::initialize() {
	if(database.lock() == nullptr) {
		suanpan_error("initialize() needs a valid domain.\n");
		return SUANPAN_FAIL;
	}

	return SUANPAN_SUCCESS;
}

int Integrator::process_load() const { return database.lock()->process_load(); }

int Integrator::process_constraint() const { return database.lock()->process_constraint(); }

int Integrator::process_criterion() const { return database.lock()->process_criterion(); }

void Integrator::record() const { database.lock()->record(); }

void Integrator::assemble_resistance() {
	// done in update_trail_status()
	// D->assemble_resistance();
}

void Integrator::assemble_matrix() {
	const auto& D = database.lock();
	D->assemble_stiffness();
	const auto& W = D->get_factory();
	if(W->get_nlgeom()) {
		D->assemble_geometry();
		get_stiffness(W) += get_geometry(W);
	}
}

void Integrator::update_trial_time(const double T) const { database.lock()->get_factory()->update_trial_time(T); }

void Integrator::update_incre_time(const double T) const { database.lock()->get_factory()->update_incre_time(T); }

int Integrator::update_trial_status() { return database.lock()->update_trial_status(); }

int Integrator::update_incre_status() { return database.lock()->update_incre_status(); }

void Integrator::erase_machine_error() const { database.lock()->erase_machine_error(); }

void Integrator::commit_status() const { database.lock()->commit_status(); }

void Integrator::clear_status() const { database.lock()->clear_status(); }

void Integrator::reset_status() const { database.lock()->reset_status(); }
