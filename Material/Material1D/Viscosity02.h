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
 * @class Viscosity02
 * @brief A 1D Viscosity class.
 * @author tlc
 * @date 01/03/2019
 * @version 0.2.0
 * @file Viscosity02.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef VISCOSITY02_H
#define VISCOSITY02_H

#include <Material/Material1D/Material1D.h>

class Viscosity02 final : public Material1D {
	const double alpha;
	const double damping_a, damping_b, damping_c, damping_d;
	const double gap_a, gap_b;

	const double limit;                                     // y=ax+b
	const double a = (alpha - 1.) * pow(limit, alpha - 2.); // a
	const double b = (2. - alpha) * pow(limit, alpha - 1.); // b

	double compute_dstrain(double, double) const;
	double compute_dstrainrate(double, double) const;
	double compute_damping_coefficient(double, double) const;
public:
	Viscosity02(unsigned, // tag
	            double,   // alpha
	            double,   // damp_a
	            double,   // damp_b
	            double,   // damp_c
	            double,   // damp_d
	            double,   // gap_a
	            double,   // gap_b
	            double    // cut-off
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	unique_ptr<Material> get_copy() override;

	int update_trial_status(const vec&, const vec&) override;

	int clear_status() override;
	int commit_status() override;
	int reset_status() override;

	vector<vec> record(OutputType) override;

	void print() override;
};

#endif

//! @}
