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
 * @class ODE_Implicit
 * @brief A ODE_Implicit class.
 * @author tlc
 * @date 22/10/2017
 * @version 0.1.0
 * @file ODE_Implicit.h
 * @addtogroup ODE_Implicit
 * @ingroup ODE_Solver
 * @{
 */

#ifndef ODE_IMPLICIT_H
#define ODE_IMPLICIT_H

#include "ODE_Solver.h"
#include <deque>

using std::deque;

class ODE_Implicit : public ODE_Solver {
protected:
	const unsigned num_step;

	const bool use_corrector;

	deque<vec> history_step;
public:
	ODE_Implicit(unsigned, unsigned, double, bool);
};

inline ODE_Implicit::ODE_Implicit(const unsigned T, const unsigned N, const double S, const bool C)
	: ODE_Solver(T, S)
	, num_step(N)
	, use_corrector(C) {}

#endif

//! @}
