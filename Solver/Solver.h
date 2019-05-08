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
 * @class Solver
 * @brief A Solver class defines solvers used in analysis.
 *
 * @author tlc
 * @date 27/07/2017
 * @version 0.2.1
 * @file Solver.h
 * @addtogroup Solver
 * @ingroup Analysis
 * @{
 */

#ifndef SOLVER_H
#define SOLVER_H

#include <Domain/Tag.h>

class Converger;
class Integrator;

class Solver : public Tag {
	shared_ptr<Converger> converger = nullptr;
	shared_ptr<Integrator> modifier = nullptr;
public:
	explicit Solver(unsigned = 0);
	virtual ~Solver();

	virtual int initialize();

	virtual int analyze() = 0;

	void set_converger(const shared_ptr<Converger>&);
	const shared_ptr<Converger>& get_converger() const;

	void set_integrator(const shared_ptr<Integrator>&);
	const shared_ptr<Integrator>& get_integrator() const;
};

#endif

//! @}