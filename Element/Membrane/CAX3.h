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
 * @class CAX3
 * @brief The CAX3 class defines CAX3 element.
 * @author tlc
 * @date 13/02/2019
 * @version 0.1.0
 * @file CAX3.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CAX3_H
#define CAX3_H

#include <Element/MaterialElement.h>

class CAX3 final : public MaterialElement {
	static const unsigned m_node, m_dof, m_size;

	const double area = 0.; // area
	const double weight = 0.;

	mat inv_coor, strain_mat;

	unique_ptr<Material> m_material; // store material model
public:
	CAX3(unsigned, uvec&&, unsigned, bool = false);

	void initialize(const shared_ptr<DomainBase>&) override;

	int update_status() override;

	int commit_status() override;
	int clear_status() override;
	int reset_status() override;

	void print() override;

#ifdef SUANPAN_VTK
	void Setup() override;
	void GetDisplacement(vtkSmartPointer<vtkDoubleArray>&) override;
	void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
