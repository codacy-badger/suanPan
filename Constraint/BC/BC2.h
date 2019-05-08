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
 * @class BC2
 * @brief A BC2 class handles boundary conditions.
 *
 * The BC2 class is in charge of applying boundary conditions to the system. The
 * BC2 class only takes care of homogeneous Dirichlet conditions. Non-homogeneous
 * displacement boundary conditions are treated as Load so that can be solved
 * iteratively. Others are handled by general constraint class such as MPC. The
 * BC2 class stores the boundary condition category, type, node(s) and
 * corresponding DoF(s). The Domain invokes `process(const shared_ptr<Domain>&)`
 * method to modify the global stiffness matrix.
 *
 * @author tlc
 * @date 23/07/2017
 * @version 0.1.0
 * @file BC2.h
 * @addtogroup Constraint
 * @{
 */

#ifndef BC2_H
#define BC2_H

#include <Constraint/BC/BC.h>

class BC2 final : public BC {
public:
	using BC::BC;

	int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
