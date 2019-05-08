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
 * @class MinDisplacement
 * @brief A MinDisplacement class.
 *
 * The MinDisplacement class.
 *
 * @author tlc
 * @date 27/09/2017
 * @version 0.1.0
 * @file MinDisplacement.h
 * @addtogroup Criterion
 * @{
 */

#ifndef MINDISPLACEMENT_H
#define MINDISPLACEMENT_H

#include <Constraint/Criterion/Criterion.h>

class MinDisplacement : public Criterion {
	unsigned node, dof;
	double limit;
public:
	explicit MinDisplacement(unsigned = 0, // tag
	                         unsigned = 0, // step tag
	                         unsigned = 0, // node tag
	                         unsigned = 0, // dof tag
	                         double = 0.   // limit
	);

	int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
