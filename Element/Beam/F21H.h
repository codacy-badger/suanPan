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
 * @class F21H
 * @brief The F21H class.
 * @author tlc
 * @date 11/10/2017
 * @version 0.2.1
 * @file F21H.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef F21H_H
#define F21H_H

#include <Element/SectionElement.h>
#include <Element/Utility/B2DC.h>

class F21H final : public SectionElement {
	static const unsigned b_node, b_dof, b_size;
	static const unsigned max_iteration;
	static const double tolerance;

	struct IntegrationPoint final {
		double coor, weight;
		unique_ptr<Section> b_section;
		mat B;
		IntegrationPoint(double, double, unique_ptr<Section>&&);
	};

	const double hinge_length;

	const double length = 0.;

	vector<IntegrationPoint> int_pt, elastic_int_pt;

	unique_ptr<Orientation> b_trans;

	mat initial_local_flexibility, elastic_local_flexibility, elastic_section_flexibility;

	mat current_local_flexibility, trial_local_flexibility;
	vec current_local_deformation, trial_local_deformation;
	vec current_local_resistance, trial_local_resistance;
public:
	F21H(unsigned,    // tag
	     uvec&&,      // node tag
	     unsigned,    // section tags
	     double = .2, // hinge length
	     bool = false // nonliear geometry switch
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	int update_status() override;

	int clear_status() override;
	int commit_status() override;
	int reset_status() override;

	vector<vec> record(OutputType) override;

	void print() override;
};

#endif

//! @}
