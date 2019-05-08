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
 * @class T3D2S
 * @brief The T3D2S class handles both linear and nonlinear problems by using a corotational transformation.
 * @author tlc
 * @date 14/08/2017
 * @version 0.1.1
 * @file T3D2S.h
 * @addtogroup Truss
 * @ingroup Element
 * @{
 */

#ifndef T3D2S_H
#define T3D2S_H

#include <Element/SectionElement.h>
#include <Element/Utility/T3DC.h>

class T3D2S final : public SectionElement {
	static const unsigned t_node, t_dof, t_size;

	const double length = 0.; // length of the element

	unique_ptr<Section> t_section;   // section model
	unique_ptr<Orientation> t_trans; // transformation

	const bool log_strain;
public:
	T3D2S(unsigned,     // tag
	      uvec&&,       // node tag
	      unsigned,     // material tag
	      bool = false, // nonlinear geometry switch
	      bool = true   // log strain swicth
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	int update_status() override;

	int commit_status() override;
	int clear_status() override;
	int reset_status() override;

	void print() override;
};

#endif

//! @}
