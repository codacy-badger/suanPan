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
 * @class BilinearPerzyna
 * @brief The BilinearPerzyna class.
 * @author tlc
 * @date 21/02/2019
 * @version 0.1.0
 * @file BilinearPerzyna.h
 * @addtogroup Material-3D
 * @{
 */

#ifndef BILINEARPERZYNA_H
#define BILINEARPERZYNA_H

#include "NonlinearPerzyna.h"

struct DataBilinearPerzyna {
	const double yield_stress, hardening_modulus;
};

class BilinearPerzyna final : DataBilinearPerzyna, public NonlinearPerzyna {
	double compute_k(double) const override;
	double compute_dk(double) const override;
public:
	BilinearPerzyna(unsigned,    // tag
	                double,      // elastic modulus
	                double,      // poisson's ratio
	                double,      // initial yield stress
	                double = 0., // hardening modulus
	                double = 0., // mu
	                double = 0., // epsilon
	                double = 0.  // density
	);

	unique_ptr<Material> get_copy() override;
};

#endif

//! @}
