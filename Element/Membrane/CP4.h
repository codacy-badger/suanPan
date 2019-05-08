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
 * @class CP4
 * @brief The CP4 class handles CPS4, CPE4, CPS4R and CPE4R elements. It is a
 * four node constant strain membrane element with optional reduced integration
 * for both plane stress and plane strain problems and optional switch for TL
 * nonlinear geometry formulation.
 * 
 * @author tlc
 * @date 10/06/2018
 * @version 0.1.3
 * @file CP4.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef CP4_H
#define CP4_H

#include <Element/MaterialElement.h>

class CP4 final : public MaterialElement {
	struct IntegrationPoint {
		vec coor;
		double weight;
		unique_ptr<Material> m_material;
		mat pn_pxy, strain_mat;
		IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
	};

	static const unsigned m_node, m_dof, m_size;

	const double thickness;

	const bool reduced_scheme;

	vector<IntegrationPoint> int_pt;

	mat hourglassing;

	static void stack_stiffness(mat&, const mat&, const mat&, double);
public:
	CP4(unsigned,     // tag
	    uvec&&,       // node tag
	    unsigned,     // material tag
	    double = 1.,  // thickness
	    bool = false, // reduced integration
	    bool = false  // nonlinear geometry switch
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
	mat GetData(OutputType) override;
	void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
