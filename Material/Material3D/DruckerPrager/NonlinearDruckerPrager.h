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
 * @class NonlinearDruckerPrager
 * @brief The NonlinearDruckerPrager class.
 * 
 * algorithm verified at 24 April 2019 by tlc
 * 
 * @author tlc
 * @date 24/04/2019
 * @version 1.0.0
 * @file NonlinearDruckerPrager.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NONLINEARDRUCKERPRAGER_H
#define NONLINEARDRUCKERPRAGER_H

#include <Material/Material3D/Material3D.h>

struct DataNonlinearDruckerPrager {
	const double elastic_modulus; // elastic modulus
	const double poissons_ratio;  // poisson's ratio
	const double shear;           // shear modulus
	const double bulk;            // bulk modulus
	const double eta_yield;
	const double eta_flow;
	const double xi;

	const double double_shear; // double shear modulus

	const double factor_a, factor_b, factor_c, factor_d;

	const bool associated = eta_yield == eta_flow;
};

class NonlinearDruckerPrager : protected DataNonlinearDruckerPrager, public Material3D {
	static const unsigned max_iteration;
	static const mat unit_dev_tensor;
	static const mat unit_x_unit;

	virtual double compute_c(double) const = 0;
	virtual double compute_dc(double) const = 0;
public:
	NonlinearDruckerPrager(unsigned,   // tag
	                       double,     // elastic modulus
	                       double,     // poisson's ratio
	                       double,     // eta_yield (hydrostatic stress related)
	                       double,     // eta_flow (dilatancy angle related)
	                       double,     // xi (cohesion related)
	                       double = 0. // density
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	double get_parameter(ParameterType) const override;

	int update_trial_status(const vec&) override;

	int clear_status() override;
	int commit_status() override;
	int reset_status() override;

	void print() override;
};

#endif

//! @}
