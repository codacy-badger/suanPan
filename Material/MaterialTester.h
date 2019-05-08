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
 * @fn MaterialTester
 * @brief A MaterialTester function.
 * @author tlc
 * @date 12/01/2019
 * @version 0.2.0
 * @file MaterialTester.h
 * @addtogroup Utility
 * @{
 */

#ifndef MATERIALTESTER_H
#define MATERIALTESTER_H

#include <Material/Material.h>

mat material_tester(const shared_ptr<Material>&, const vector<unsigned>&, const vec&);

#endif

//! @}
