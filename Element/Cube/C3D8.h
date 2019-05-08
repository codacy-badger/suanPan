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
 * @class C3D8
 * @brief The C3D8 class defines C3D8 C3D8R elements.
 * @author tlc
 * @date 13/07/2018
 * @version 0.2.0
 * @file C3D8.h
 * @addtogroup Cube
 * @ingroup Element
 * @{
 */

#ifndef C3D8_H
#define C3D8_H

#include <Element/MaterialElement.h>

class C3D8 final : public MaterialElement {
	struct IntegrationPoint {
		vec coor;
		double weight;
		unique_ptr<Material> c_material;
		mat pn_pxyz, strain_mat;
		IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
	};

	static const unsigned c_node, c_dof, c_size;

	const char int_scheme;

	vector<IntegrationPoint> int_pt;
public:
	C3D8(unsigned,    // tag
	     uvec&&,      // node tag
	     unsigned,    // material tag
	     char = 'I',  // reduced integration
	     bool = false // nonlinear geometry switch
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
