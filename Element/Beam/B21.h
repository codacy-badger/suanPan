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
 * @class B21
 * @brief The B21 class.
 * @author tlc
 * @date 11/10/2017
 * @version 0.1.1
 * @file B21.h
 * @addtogroup Beam
 * @ingroup Element
 * @{
 */

#ifndef B21_H
#define B21_H

#include <Element/SectionElement.h>
#include <Element/Utility/B2DC.h>

class B21 final : public SectionElement {
	static const unsigned b_node, b_dof, b_size;

	struct IntegrationPoint {
		double coor, weight;
		unique_ptr<Section> b_section;
		mat strain_mat;
		IntegrationPoint(double, double, unique_ptr<Section>&&);
	};

	const unsigned int_pt_num;

	const double length = 0.;

	vector<IntegrationPoint> int_pt;

	unique_ptr<Orientation> b_trans;
public:
	B21(unsigned,     // tag
	    uvec&&,       // node tags
	    unsigned,     // section tag
	    unsigned = 6, // integration points
	    bool = false  // nonliear geometry switch
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	int update_status() override;

	int commit_status() override;
	int clear_status() override;
	int reset_status() override;

	vector<vec> record(OutputType) override;

	void print() override;

#ifdef SUANPAN_VTK
	void Setup() override;
	void GetDisplacement(vtkSmartPointer<vtkDoubleArray>&) override;
	void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
