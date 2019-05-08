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
 * @class T2D2
 * @brief A T2D2 class
 * 
 * The T2D2 class handles both linear and nonlinear problems by using a optional corotational transformation.
 *  
 * @author tlc
 * @date 28/06/2018
 * @version 0.2.0
 * @file T2D2.h
 * @addtogroup Truss
 * @ingroup Element
 * @{
 */

#ifndef T2D2_H
#define T2D2_H

#include <Element/MaterialElement.h>
#include <Element/Utility/T2DC.h>

class T2D2 final : public MaterialElement {
	static const unsigned t_node, t_dof, t_size;

	const double length = 0.; // length of the element
	const double area;        // area of the cross section

	unique_ptr<Material> t_material; // material model
	unique_ptr<Orientation> t_trans; // transformation

	const bool update_area; // flag to indicate if update section area
	const bool log_strain;  // flag to indicate if use log strain
public:
	T2D2(unsigned,     // tag
	     uvec&&,       // node tag
	     unsigned,     // material tag
	     double,       // area
	     bool = false, // nonlinear geometry switch
	     bool = true,  // update area swicth
	     bool = true   // log strain swicth
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