﻿/*******************************************************************************
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
 * @class DKT4
 * @brief A DKT4 plate class.
 *
 * @author tlc
 * @date 27/03/2019
 * @version 0.1.0
 * @file DKT4.h
 * @addtogroup Plate
 * @ingroup Element
 * @{
 */

#ifndef DKT4_H
#define DKT4_H

#include <Element/MaterialElement.h>

class DKT4 final : public MaterialElement {
	static const unsigned p_node, p_dof, p_size;

	struct IntegrationPoint {
		struct SectionIntegrationPoint {
			const double eccentricity, factor;
			unique_ptr<Material> p_material;
			SectionIntegrationPoint(double, double, unique_ptr<Material>&&);
			SectionIntegrationPoint(const SectionIntegrationPoint&);
			SectionIntegrationPoint(SectionIntegrationPoint&&) noexcept = default;
			SectionIntegrationPoint& operator=(const SectionIntegrationPoint&) = delete;
			SectionIntegrationPoint& operator=(SectionIntegrationPoint&&) noexcept = delete;
			~SectionIntegrationPoint() = default;
		};

		vec coor;
		mat strain_mat;
		vector<SectionIntegrationPoint> sec_int_pt;
		explicit IntegrationPoint(vec&&);
	};

	const double thickness;
	const unsigned num_section_ip;

	vector<IntegrationPoint> int_pt;

	static field<mat> form_transform(const mat&);
public:
	DKT4(unsigned,    // element tag
	     uvec&&,      // node tag
	     unsigned,    // material tag
	     double,      // thickness
	     unsigned = 3 // integration points along thickness
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	int update_status() override;
	int clear_status() override;
	int commit_status() override;
	int reset_status() override;

#ifdef SUANPAN_VTK
	void Setup() override;
	void GetDisplacement(vtkSmartPointer<vtkDoubleArray>&) override;
	void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
