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

#ifndef VTKPARSER_H
#define VTKPARSER_H

#include <suanPan>

class DomainBase;

#ifdef SUANPAN_VTK

class vtkUnstructuredGrid;

int vtk_parser(const shared_ptr<DomainBase>&, istringstream&);

struct vtkInfo {
	OutputType type = OutputType::U;
	double scale = 1.;
	bool on_deformed = true;
	unsigned font_size = 8;
	unsigned canvas_size[2] = {500, 500};
	bool save_file = false;
	string file_name;
	string title_name;
	bool colorbar = true;
};

vtkInfo vtk_process(istringstream&);

void vtk_setup(const vtkSmartPointer<vtkUnstructuredGrid>&, const vtkInfo&);

void vtk_save(const vtkSmartPointer<vtkUnstructuredGrid>&, const vtkInfo&);

void vtk_plot_displacement(const shared_ptr<DomainBase>&, vtkInfo);

void vtk_plot_stress(const shared_ptr<DomainBase>&, vtkInfo);

#else

inline int vtk_parser(const shared_ptr<DomainBase>&, istringstream&) { return 0; }

#endif

#endif

//! @}
