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
 * @class CP4I
 * @brief The CP4I class handles CPS4I, CPE4I.
 * 
 * @author tlc
 * @date 25/07/2018
 * @version 0.1.0
 * @file CP4I.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CP4I_H
#define CP4I_H

#include <Element/MaterialElement.h>

class CP4I final : public MaterialElement {
	static const unsigned m_node, m_dof, m_size;

	struct IntegrationPoint {
		vec coor;
		double weight;
		unique_ptr<Material> m_material;
		mat pn_pxy, B1, B2;
		IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
	};

	const double thickness;

	vector<IntegrationPoint> int_pt;

	static void stack_stiffness(mat&, const mat&, const mat&, double);
	static void stack_stiffness_incompatible(mat&, const mat&, const mat&, double);
public:
	CP4I(unsigned,   // tag
	     uvec&&,     // node tag
	     unsigned,   // material tag
	     double = 1. // thickness
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
