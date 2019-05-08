/*******************************************************************************
 * Copyright (C) 2017-2019 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
/**
 * @class ODE
 * @brief An ODE class in charge of defining an ODE system.
 * 
 * An ordinary differential equation (ODE) (system) is normally expressed as
 * \f{gather}{
 * y'=f\left(t,y\right),
 * \f}
 * where \f$t\f$ is the generalized time which is always a scalar, \f$y\f$ are
 * independent variables and \f$y'\f$ are the corresponding derivatives. There
 * could be only one independent variable, or more often, multiple variables.
 * 
 * @author tlc
 * @date 22/02/2019
 * @version 0.1.1
 * @file ODE.h
 */

#ifndef ODE_H
#define ODE_H

#include <Domain/Tag.h>

class ODE_Solver;
class DomainBase;

class ODE : public Tag {
	static const unsigned max_iteration;

	const double abs_tolerance = 1E-10;
	const double rel_tolerance = 1E-6;

	double error = 0.;

	double backup_time = 0.;
	double current_time = 0.;
	double trial_time = 0.;
	double incre_time = 0.;

	vec backup_variable;
	vec current_variable;
	vec trial_variable;
	vec incre_variable;

	shared_ptr<ODE_Solver> ode_solver;
public:
	ODE(unsigned, unsigned);
	ODE(unsigned, unsigned, shared_ptr<ODE_Solver>);
	ODE(const ODE&) = default;
	ODE(ODE&&) = default;
	ODE& operator=(const ODE&) = delete;
	ODE& operator=(ODE&&) = delete;
	virtual ~ODE();

	void initialize(const shared_ptr<DomainBase>&);

	virtual unique_ptr<ODE> get_copy() = 0;

	virtual vec eval(double, const vec&) = 0;

	vec operator()(double, const vec&);

	int update_incre_status(double);
	int update_trial_status(double);

	int commit_status();
	int clear_status();
	int reset_status();

	void set_error(double);
	double get_error() const;

	void set_abs_tolerance(double) const;
	double get_abs_tolerance() const;
	void set_rel_tolerance(double) const;
	double get_rel_tolerance() const;

	void set_ode_solver(shared_ptr<ODE_Solver>);
	const shared_ptr<ODE_Solver>& get_ode_solver() const;

	void set_current_time(double);
	void set_trial_time(double);
	void set_incre_time(double);

	double get_current_time() const;
	double get_trial_time() const;
	double get_incre_time() const;

	void set_current_variable(const vec&);
	void set_trial_variable(const vec&);
	void set_incre_variable(const vec&);

	const vec& get_current_variable() const;
	const vec& get_trial_variable() const;
	const vec& get_incre_variable() const;
};

#endif
