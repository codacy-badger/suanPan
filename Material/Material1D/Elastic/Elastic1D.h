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
 * @class Elastic1D
 * @brief A 1-D Elastic class.
 * @author tlc
 * @date 27/07/2017
 * @file Elastic1D.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef ELASTIC1D_H
#define ELASTIC1D_H

#include <Material/Material1D/Material1D.h>

struct DataElastic1D {
	const double elastic_modulus; // elastic modulus
};

class Elastic1D final : DataElastic1D, public Material1D {
public:
	Elastic1D(unsigned,   // tag 
	          double,     // elastic modulus
	          double = 0. // density
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	unique_ptr<Material> get_copy() override;

	int update_trial_status(const vec&) override;

	int clear_status() override;
	int commit_status() override;
	int reset_status() override;

	void print() override;
};

#endif

//! @}