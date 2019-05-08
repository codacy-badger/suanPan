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
 * @class F31
 * @brief The F31 class.
 * 
 * Reference: https://doi.org/10.1016/0045-7949(95)00103-N
 * 
 * @author tlc
 * @date 27/07/2018
 * @version 0.1.0
 * @file F31.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef F31_H
#define F31_H

#include <Element/SectionElement.h>
#include <Element/Utility/B2DC.h>

class F31 final : public SectionElement {
	static const unsigned b_node, b_dof, b_size;
	static const span b_span;
	static const unsigned max_iteration;
	static const double tolerance;

	struct IntegrationPoint {
		double coor, weight;
		unique_ptr<Section> b_section;
		mat strain_mat;
		IntegrationPoint(double, double, unique_ptr<Section>&&);
	};

	const unsigned int_pt_num, orientation_tag;

	const double length = 0., torsion_stiff = 0.;

	vector<IntegrationPoint> int_pt;

	unique_ptr<Orientation> b_trans;

	mat initial_local_flexibility;
	mat current_local_flexibility, trial_local_flexibility;
	vec current_local_deformation, trial_local_deformation;
	vec current_local_resistance, trial_local_resistance;
public:
	F31(unsigned,     // tag
	    uvec&&,       // node tag
	    unsigned,     // section tag
	    unsigned,     // orientation tag
	    unsigned = 6, // integration points
	    bool = false  // nonliear geometry switch
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
