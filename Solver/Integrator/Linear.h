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
 * @class Linear
 * @brief A Linear class defines a solver using Linear algorithm.
 *
 * @author tlc
 * @date 25/08/2017
 * @version 0.1.1
 * @file Linear.h
 * @addtogroup Integrator
 * @{
 */

#ifndef LINEAR_H
#define LINEAR_H

#include "Integrator.h"

class Linear final : public Integrator {
public:
	using Integrator::Integrator;

	void assemble_matrix() override;
};

#endif

//! @}
