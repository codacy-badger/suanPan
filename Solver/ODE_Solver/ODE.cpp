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

#include "ODE.h"
#include <Domain/DomainBase.h>
#include <Solver/ODE_Solver/ODE_Solver.h>

const unsigned ODE::max_iteration = 2000;

ODE::ODE(const unsigned T, const unsigned D)
	: ODE(T, D, nullptr) {}

ODE::ODE(const unsigned T, const unsigned D, shared_ptr<ODE_Solver> O)
	: Tag(T)
	, backup_variable(D, fill::zeros)
	, current_variable(D, fill::zeros)
	, trial_variable(D, fill::zeros)
	, incre_variable(D, fill::zeros)
	, ode_solver(std::move(O)) { suanpan_debug("ODE %u ctor called.\n", T); }

ODE::~ODE() { suanpan_debug("ODE %u dtor called.\n", get_tag()); }

void ODE::initialize(const shared_ptr<DomainBase>&) { if(ode_solver == nullptr) disable(); }

vec ODE::operator()(const double t_time, const vec& t_variable) { return eval(t_time, t_variable); }

int ODE::update_incre_status(const double i_time) { return update_trial_status(i_time + current_time); }

int ODE::update_trial_status(const double t_time) {
	current_time = backup_time;
	current_variable = backup_variable;

	set_trial_time(t_time);

	auto time_left = incre_time;
	auto step = .1 * time_left;

	unsigned counter = 0;
	while(++counter < max_iteration) {
		if(ode_solver->analyze(this) == SUANPAN_SUCCESS && (fabs(error) <= abs_tolerance || fabs(error) <= rel_tolerance * norm(trial_variable))) {
			current_time = trial_time;
			current_variable = trial_variable;
			if((time_left -= step) <= 0.) break;
		}
		set_incre_time(step = std::min(step * std::min(4., std::max(.1, .95 * pow(fabs(.5 * abs_tolerance / error) + datum::eps, ode_solver->get_scale()))), time_left));
	}

	suanpan_extra_debug("ODE local iteration counter: %u.\n", counter);

	if(max_iteration == counter) return SUANPAN_FAIL;

	return SUANPAN_SUCCESS;
}

int ODE::commit_status() {
	backup_time = current_time = trial_time;
	incre_time = 0.;
	backup_variable = current_variable = trial_variable;
	incre_variable.zeros();
	return SUANPAN_SUCCESS;
}

int ODE::clear_status() {
	trial_time = backup_time = current_time = incre_time = 0.;
	trial_variable = backup_variable = current_variable = incre_variable.zeros();
	return SUANPAN_SUCCESS;
}

int ODE::reset_status() {
	trial_time = current_time = backup_time;
	incre_time = 0.;
	trial_variable = current_variable = backup_variable;
	incre_variable.zeros();
	return SUANPAN_SUCCESS;
}

void ODE::set_error(const double E) { error = E; }

double ODE::get_error() const { return error; }

void ODE::set_abs_tolerance(const double T) const { access::rw(abs_tolerance) = T; }

double ODE::get_abs_tolerance() const { return abs_tolerance; }

void ODE::set_rel_tolerance(const double T) const { access::rw(rel_tolerance) = T; }

double ODE::get_rel_tolerance() const { return rel_tolerance; }

void ODE::set_ode_solver(shared_ptr<ODE_Solver> solver) { ode_solver = std::move(solver); }

const shared_ptr<ODE_Solver>& ODE::get_ode_solver() const { return ode_solver; }

void ODE::set_current_time(const double c_time) {
	trial_time = current_time = c_time;
	incre_time = 0.;
}

void ODE::set_trial_time(const double t_time) { incre_time = (trial_time = t_time) - current_time; }

void ODE::set_incre_time(const double i_time) { trial_time = current_time + (incre_time = i_time); }

double ODE::get_current_time() const { return current_time; }

double ODE::get_trial_time() const { return trial_time; }

double ODE::get_incre_time() const { return incre_time; }

void ODE::set_current_variable(const vec& c_variable) {
	trial_variable = current_variable = backup_variable = c_variable;
	incre_variable.zeros(c_variable.size());
}

void ODE::set_trial_variable(const vec& t_variable) { incre_variable = (trial_variable = t_variable) - current_variable; }

void ODE::set_incre_variable(const vec& i_variable) { trial_variable = current_variable + (incre_variable = i_variable); }

const vec& ODE::get_current_variable() const { return current_variable; }

const vec& ODE::get_trial_variable() const { return trial_variable; }

const vec& ODE::get_incre_variable() const { return incre_variable; }
