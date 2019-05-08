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
 * @class Damper01
 * @brief A Damper01 class.
 * 
 * Using quandrant damper and displacement and velocity as basic input.
 * 
 * @author tlc
 * @date 23/07/2018
 * @file Damper01.h
 * @addtogroup Special
 * @ingroup Element
 * @{
 */

#ifndef DAMPER01_H
#define DAMPER01_H

#include <Element/Element.h>

class Damper01 final : public Element {
	static unsigned d_node, d_dof, d_size;
	static uvec IS, JS;

	const unsigned damper_tag;

	const vec direction_cosine;

	unique_ptr<Material> damper;
public:
	Damper01(unsigned, // tag
	         uvec&&,   // node tag
	         unsigned  // damper tag
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	int update_status() override;

	int commit_status() override;
	int clear_status() override;
	int reset_status() override;

	vector<vec> record(OutputType) override;

	void print() override;
};

#endif

//! @}
