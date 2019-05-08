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
 * @class Viscosity01
 * @brief A 1D Elastic class.
 * @author tlc
 * @date 01/03/2019
 * @version 0.3.0
 * @file Viscosity01.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef VISCOSITY01_H
#define VISCOSITY01_H

#include <Material/Material1D/Material1D.h>

struct DataViscosity01 {
	const double alpha, damping, limit;

	const double a = (alpha - 1.) * pow(limit, alpha - 2.); // a
	const double b = (2. - alpha) * pow(limit, alpha - 1.); // b
};

class Viscosity01 final : DataViscosity01, public Material1D {
public:
	Viscosity01(unsigned, // tag
	            double,   // alpha
	            double,   // damp coefficient
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
