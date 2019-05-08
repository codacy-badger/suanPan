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
 * @class ODE_Explicit
 * @brief A ODE_Explicit class.
 * @author tlc
 * @date 22/02/2019
 * @version 0.1.0
 * @file ODE_Explicit.h
 * @addtogroup ODE_Explicit
 * @ingroup ODE_Solver
 * @{
 */

#ifndef ODE_EXPLICIT_H
#define ODE_EXPLICIT_H

#include "ODE_Solver.h"

class ODE_Explicit : public ODE_Solver {
public:
	using ODE_Solver::ODE_Solver;
};

#endif

//! @}
