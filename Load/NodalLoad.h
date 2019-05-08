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
 * @class NodalLoad
 * @brief A NodalLoad class.
 *
 * The NodalLoad class is in charge of handling concentrated load.
 *
 * @author tlc
 * @date 23/07/2017
 * @version 0.1.0
 * @file NodalLoad.h
 * @addtogroup Load
 * @{
 */

#ifndef NODALLOAD_H
#define NODALLOAD_H

#include <Load/Load.h>

class NodalLoad final : public Load {
public:
	explicit NodalLoad(unsigned = 0, // tag
	                   unsigned = 0, // start step tag
	                   double = 0.,  // maginitude
	                   uvec&& = {},  // node tags
	                   unsigned = 0, // dof tag
	                   unsigned = 0  // amplitude tag
	);
	NodalLoad(unsigned,    // tag
	          unsigned,    // start step tag
	          double,      // maginitude
	          uvec&&,      // node tags
	          uvec&&,      // dof tags
	          unsigned = 0 // amplitude tag
	);

	int process(const shared_ptr<DomainBase>&) override;
};

#endif

//! @}
