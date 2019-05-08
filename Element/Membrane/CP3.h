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
 * @class CP3
 * @brief The CP3 class defines CPS3 CPE3 elements.
 * @author tlc
 * @date 13/07/2018
 * @version 0.3.0
 * @file CP3.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CP3_H
#define CP3_H

#include <Element/MaterialElement.h>

class CP3 final : public MaterialElement {
	static const unsigned m_node, m_dof, m_size;

	const double thickness; // thickness
	const double area = 0.; // area

	mat pn_pxy, strain_mat;

	unique_ptr<Material> m_material; // store material model

	static void stack_stiffness(mat&, const mat&, const mat&, double);
public:
	CP3(unsigned, uvec&&, unsigned, double = 1., bool = false);

	void initialize(const shared_ptr<DomainBase>&) override;

	int update_status() override;

	int commit_status() override;
	int clear_status() override;
	int reset_status() override;

	void print() override;

#ifdef SUANPAN_VTK
	void Setup() override;
	void GetDisplacement(vtkSmartPointer<vtkDoubleArray>&) override;
	mat GetData(OutputType) override;
	void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
