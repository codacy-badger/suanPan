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
 * @class PS
 * @brief The PS class defines a Pian--Sumihara membrane element.
 * @author tlc
 * @date 30/07/2017
 * @version 0.1.0
 * @file PS.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef PS_H
#define PS_H

#include <Element/MaterialElement.h>

class PS final : public MaterialElement {
	struct IntegrationPoint {
		vec coor;
		double weight;
		unique_ptr<Material> m_material;
		mat poly_strain;
		IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
	};

	static const unsigned m_node, m_dof, m_size;

	const double thickness;

	vector<IntegrationPoint> int_pt;

	static mat form_transformation(const mat&);
public:
	PS(unsigned,   // unique tag
	   uvec&&,     // node tags
	   unsigned,   // material tag
	   double = 1. // thickness
	);

	void initialize(const shared_ptr<DomainBase>& = nullptr) override;

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
