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
 * @class NonlinearJ2
 * @brief The NonlinearJ2 class.
 * 
 * algorithm verified at 24 April 2019 by tlc
 * 
 * @author tlc
 * @date 24/04/2019
 * @version 1.0.0
 * @file NonlinearJ2.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef NONLINEARJ2_H
#define NONLINEARJ2_H

#include <Material/Material3D/Material3D.h>

struct DataNonlinearJ2 {
	const double elastic_modulus; // elastic modulus
	const double poissons_ratio;  // poisson's ratio

	const double shear_modulus;       // shear modulus
	const double double_shear;        // double shear modulus
	const double square_double_shear; // square shear modulus
};

class NonlinearJ2 : protected DataNonlinearJ2, public Material3D {
	static const double two_third;
	static const double root_two_third;
	static const mat unit_dev_tensor;
	static const unsigned max_iteration;

	virtual double compute_k(double) const = 0;
	virtual double compute_dk(double) const = 0;
	virtual double compute_h(double) const = 0;
	virtual double compute_dh(double) const = 0;
public:
	NonlinearJ2(unsigned,   // tag
	            double,     // elastic modulus
	            double,     // poisson's ratio
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