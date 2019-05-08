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

#include "Step.h"
#include <Converger/RelIncreDisp.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Solver.h>

Step::Step(const unsigned T, const double P)
	: Tag(T)
	, time_period(P) { suanpan_debug("Step %u ctor() called.\n", T); }

Step::~Step() { suanpan_debug("Step %u dtor() called.\n", get_tag()); }

bool Step::is_updated() const { return updated; }

int Step::initialize() {
	const auto& t_domain = database.lock();

	if(converger_tag != 0 && t_domain->find_converger(converger_tag)) tester = t_domain->get_converger(converger_tag);
	else if(t_domain->get_current_converger_tag() != 0) tester = t_domain->get_current_converger();

	if(integrator_tag != 0 && t_domain->find_integrator(integrator_tag)) modifier = t_domain->get_integrator(integrator_tag);
	else if(t_domain->get_current_integrator_tag() != 0) modifier = t_domain->get_current_integrator();

	if(solver_tag != 0 && t_domain->find_solver(solver_tag)) solver = t_domain->get_solver(solver_tag);
	else if(t_domain->get_current_solver_tag() != 0) solver = t_domain->get_current_solver();

	if(tester == nullptr) tester = make_shared<RelIncreDisp>();

	factory = t_domain->get_factory();

	return 0;
}

void Step::set_domain(const weak_ptr<DomainBase>& D) {
	if(database.lock() == D.lock()) return;
	database = D;
	updated = false;
}

const weak_ptr<DomainBase>& Step::get_domain() const { return database; }

void Step::set_factory(const shared_ptr<Factory<double>>& F) {
	if(factory == F) return;
	factory = F;
	updated = false;
}

const shared_ptr<Factory<double>>& Step::get_factory() const { return factory; }

void Step::set_solver_tag(const unsigned T) { solver_tag = T; }

void Step::set_solver(const shared_ptr<Solver>& S) {
	if(solver == S) return;
	solver = S;
	updated = false;
}

const shared_ptr<Solver>& Step::get_solver() const { return solver; }

void Step::set_converger_tag(const unsigned T) { converger_tag = T; }

void Step::set_converger(const shared_ptr<Converger>& C) {
	if(tester == C) return;
	tester = C;
	updated = false;
}

const shared_ptr<Converger>& Step::get_converger() const { return tester; }

void Step::set_integrator_tag(const unsigned T) { integrator_tag = T; }

void Step::set_integrator(const shared_ptr<Integrator>& G) {
	if(modifier == G) return;
	modifier = G;
	updated = false;
}

const shared_ptr<Integrator>& Step::get_integrator() const { return modifier; }

void Step::set_time_perid(const double T) {
	if(time_period == T) return;
	time_period = T;
	updated = false;
	const auto tmp_iteration = static_cast<int>(floor(time_period / ini_step_size)) + 1;
	if(tmp_iteration > static_cast<int>(max_substep) && max_substep != 0) {
		if(tmp_iteration > static_cast<int>(std::numeric_limits<unsigned>::max())) {
			suanpan_warning("set_ini_step_size() exceeds limits.\n");
			set_max_substep(std::numeric_limits<unsigned>::max());
		} else set_max_substep(tmp_iteration);
	}
}

double Step::get_time_period() const { return time_period; }

void Step::set_ini_step_size(const double T) {
	if(ini_step_size == T) return;
	ini_step_size = T > time_period ? time_period : T;
	updated = false;
	const auto tmp_iteration = static_cast<int>(floor(time_period / ini_step_size)) + 1;
	if(tmp_iteration > static_cast<int>(max_substep) && max_substep != 0) set_max_substep(tmp_iteration);
}

void Step::set_min_step_size(const double T) {
	if(min_step_size == T) return;
	min_step_size = T;
	updated = false;
}

void Step::set_max_step_size(const double T) {
	if(max_step_size == T) return;
	max_step_size = T;
	updated = false;
}

void Step::set_max_substep(const unsigned M) {
	if(max_substep == M) return;
	max_substep = M;
	updated = false;
}

double Step::get_ini_step_size() const { return ini_step_size; }

double Step::get_min_step_size() const { return min_step_size; }

double Step::get_max_step_size() const { return max_step_size; }

unsigned Step::get_max_substep() const { return max_substep; }

bool Step::is_fixed_step_size() const { return fixed_step_size; }

void Step::set_fixed_step_size(const bool B) {
	if(fixed_step_size == B) return;
	fixed_step_size = B;
	updated = false;
}

bool Step::is_symm() const { return symm_mat; }

bool Step::is_band() const { return band_mat; }

bool Step::is_sparse() const { return sparse_mat; }

void Step::set_symm(const bool B) {
	if(symm_mat == B) return;
	access::rw(symm_mat) = B;
	updated = false;
}

void Step::set_band(const bool B) {
	if(band_mat == B) return;
	access::rw(band_mat) = B;
	updated = false;
}

void Step::set_sparse(const bool B) {
	if(sparse_mat == B) return;
	access::rw(sparse_mat) = B;
	updated = false;
}
