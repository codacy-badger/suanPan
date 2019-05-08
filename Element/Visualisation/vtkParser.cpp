////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2019 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "vtkParser.h"

#ifdef SUANPAN_VTK

#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)  // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
VTK_MODULE_INIT(vtkInteractionStyle)  // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
VTK_MODULE_INIT(vtkRenderingFreeType) // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)

#include <vtkActor.h>
#include <vtkColorTransferFunction.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>

vtkInfo vtk_process(istringstream& command) {
	vtkInfo config;

	string keyword;

	while(!command.eof() && get_input(command, keyword))
		if(is_equal(keyword, "scale") && !get_input(command, config.scale)) config.scale = 1.;
		else if(is_equal(keyword, "deformed")) config.on_deformed = true;
		else if(is_equal(keyword, "undeformed")) config.on_deformed = false;
		else if(is_equal(keyword, "type") && get_input(command, keyword)) config.type = to_list(keyword.c_str());
		else if(is_equal(keyword, "fontsize") && !get_input(command, config.font_size)) config.font_size = 8;
		else if(is_equal(keyword, "save") && get_input(command, config.file_name)) config.save_file = true;
		else if(is_equal(keyword, "nobar")) config.colorbar = false;
		else if(is_equal(keyword, "size")) {
			if(!get_input(command, config.canvas_size[0])) config.canvas_size[0] = 500;
			if(!get_input(command, config.canvas_size[1])) config.canvas_size[1] = 500;
		}

	return config;
}

void vtk_setup(const vtkSmartPointer<vtkUnstructuredGrid>& grid, const vtkInfo& config) {
	auto color = vtkSmartPointer<vtkNamedColors>::New();
	auto table = vtkSmartPointer<vtkLookupTable>::New();
	auto func = vtkSmartPointer<vtkColorTransferFunction>::New();
	auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	auto bar = vtkSmartPointer<vtkScalarBarActor>::New();
	auto actor = vtkSmartPointer<vtkActor>::New();
	auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	auto renderer = vtkSmartPointer<vtkRenderer>::New();
	auto window = vtkSmartPointer<vtkRenderWindow>::New();

	func->SetColorSpaceToDiverging();
	func->AddRGBPoint(0., .230, .299, .754);
	func->AddRGBPoint(1., .706, .016, .150);
	table->SetNumberOfTableValues(256);
	table->Build();
	double rgb[4] = {0., 0., 0., 1.};
	for(auto I = 0; I < 256; ++I) {
		func->GetColor(double(I) / 256., rgb);
		table->SetTableValue(I, rgb);
	}

	mapper->SetInputDataObject(grid);
	mapper->SetLookupTable(table);
	mapper->SetScalarRange(grid->GetPointData()->GetScalars()->GetRange());

	bar->SetLookupTable(mapper->GetLookupTable());

	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(color->GetColor3d("DodgerBlue").GetData());
	actor->GetProperty()->EdgeVisibilityOn();
	actor->GetProperty()->SetOpacity(1.);

	renderer->AddActor(actor);
	if(config.colorbar) renderer->AddActor2D(bar);
	renderer->SetBackground(color->GetColor3d("Grey").GetData());
	renderer->ResetCamera();
	renderer->ResetCameraClippingRange();

	window->AddRenderer(renderer);
	window->SetSize(config.canvas_size[0], config.canvas_size[1]);
	window->Render();
	window->SetWindowName(config.title_name.c_str());

	interactor->SetRenderWindow(window);
	interactor->Start();
}

void vtk_save(const vtkSmartPointer<vtkUnstructuredGrid>& grid, const vtkInfo& config) {
	auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetInputData(grid);
	writer->SetFileName(config.file_name.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();
	suanpan_info("plot is written to %s.\n", config.file_name.c_str());
}

int vtk_parser(const shared_ptr<DomainBase>& domain, istringstream& command) {
	auto plot_info = vtk_process(command);
	if(!plot_info.on_deformed) plot_info.scale = 0.;

	if(plot_info.type == OutputType::U1 || plot_info.type == OutputType::U2 || plot_info.type == OutputType::U3 || plot_info.type == OutputType::UR1 || plot_info.type == OutputType::UR2 || plot_info.type == OutputType::UR3) domain->insert(make_shared<std::thread>(&vtk_plot_displacement, std::cref(domain), plot_info));
	else if(plot_info.type == OutputType::S11 || plot_info.type == OutputType::S22 || plot_info.type == OutputType::S33 || plot_info.type == OutputType::S12 || plot_info.type == OutputType::S23 || plot_info.type == OutputType::S13) domain->insert(make_shared<std::thread>(&vtk_plot_stress, std::cref(domain), plot_info));

	std::this_thread::sleep_for(std::chrono::milliseconds(500));

	return SUANPAN_SUCCESS;
}

void vtk_plot_displacement(const shared_ptr<DomainBase>& domain, vtkInfo config) {
	auto& t_node_pool = domain->get_node_pool();
	auto& t_element_pool = domain->get_element_pool();

	auto idx = 0;
	if(config.type == OutputType::U1) idx = 0;
	else if(config.type == OutputType::U2) idx = 1;
	else if(config.type == OutputType::U3) idx = 2;
	else if(config.type == OutputType::UR1) idx = 3;
	else if(config.type == OutputType::UR2) idx = 4;
	else if(config.type == OutputType::UR3) idx = 5;

	string title = "displacement componenet ";
	title += to_char(config.type);
	title += " plot";
	config.title_name = title;

	auto max_node = unsigned(t_node_pool.size());
	for(const auto& I : t_node_pool) if(I->get_tag() > max_node) max_node = I->get_tag();

	auto data = vtkSmartPointer<vtkDoubleArray>::New();
	auto node = vtkSmartPointer<vtkPoints>::New();

	data->SetNumberOfComponents(6);
	data->SetNumberOfTuples(++max_node);
	node->SetNumberOfPoints(max_node);

	for(unsigned I = 0; I < max_node; ++I) {
		node->SetPoint(I, 0., 0., 0.);
		data->SetTuple6(I, 0., 0., 0., 0., 0., 0.);
	}

	data->SetName("U");
	data->SetComponentName(0, "U1");
	data->SetComponentName(1, "U2");
	data->SetComponentName(2, "U3");
	data->SetComponentName(3, "UR1");
	data->SetComponentName(4, "UR2");
	data->SetComponentName(5, "UR3");

	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->Allocate(t_element_pool.size());
	suanpan_for_each(t_element_pool.cbegin(), t_element_pool.cend(), [&](const shared_ptr<Element>& t_element) {
			t_element->SetDeformation(node, config.scale);
			t_element->GetDisplacement(data);
			grid->InsertNextCell(t_element->GetCell()->GetCellType(), t_element->GetCell()->GetPointIds());
		});

	grid->SetPoints(node);

	if(config.save_file) {
		grid->GetPointData()->SetScalars(data);
		grid->GetPointData()->SetActiveScalars("U");
		vtk_save(grid, config);
	} else {
		auto sub_data = vtkSmartPointer<vtkDoubleArray>::New();

		sub_data->SetNumberOfTuples(data->GetNumberOfTuples());
		sub_data->CopyComponent(0, data, idx);

		grid->GetPointData()->SetScalars(sub_data);

		vtk_setup(grid, config);
	}
}

void vtk_plot_stress(const shared_ptr<DomainBase>& domain, vtkInfo config) {
	auto& t_node_pool = domain->get_node_pool();
	auto& t_element_pool = domain->get_element_pool();

	auto idx = 0;
	if(config.type == OutputType::S11) idx = 0;
	else if(config.type == OutputType::S22) idx = 1;
	else if(config.type == OutputType::S33) idx = 2;
	else if(config.type == OutputType::S12) idx = 3;
	else if(config.type == OutputType::S23) idx = 4;
	else if(config.type == OutputType::S13) idx = 5;

	string title = "stress componenet ";
	title += to_char(config.type);
	title += " plot";
	config.title_name = title;

	auto max_node = unsigned(t_node_pool.size());
	for(const auto& I : t_node_pool) if(I->get_tag() > max_node) max_node = I->get_tag();

	auto data = vtkSmartPointer<vtkDoubleArray>::New();
	auto node = vtkSmartPointer<vtkPoints>::New();

	data->SetNumberOfComponents(6);
	data->SetNumberOfTuples(++max_node);
	node->SetNumberOfPoints(max_node);

	for(unsigned I = 0; I < max_node; ++I) {
		node->SetPoint(I, 0., 0., 0.);
		data->SetTuple6(I, 0., 0., 0., 0., 0., 0.);
	}

	data->SetName("S");
	data->SetComponentName(0, "S11");
	data->SetComponentName(1, "S22");
	data->SetComponentName(2, "S33");
	data->SetComponentName(3, "S12");
	data->SetComponentName(4, "S23");
	data->SetComponentName(5, "S13");

	mat stress(6, max_node, fill::zeros);
	vec counter(max_node, fill::zeros);

	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->Allocate(t_element_pool.size());
	std::for_each(t_element_pool.cbegin(), t_element_pool.cend(), [&](const shared_ptr<Element>& t_element) {
		t_element->SetDeformation(node, config.scale);
		auto& t_encoding = t_element->get_node_encoding();
		counter(t_encoding) += 1.;
		stress.cols(t_encoding) += t_element->GetData(OutputType::S);
		grid->InsertNextCell(t_element->GetCell()->GetCellType(), t_element->GetCell()->GetPointIds());
	});

	for(unsigned I = 0; I < max_node; ++I)
		if(0. != counter(I)) {
			stress.col(I) /= counter(I);
			data->SetTuple(I, stress.colptr(I));
		}

	grid->SetPoints(node);

	if(config.save_file) {
		grid->GetPointData()->SetScalars(data);
		grid->GetPointData()->SetActiveScalars("S");
		vtk_save(grid, config);
	} else {
		auto sub_data = vtkSmartPointer<vtkDoubleArray>::New();

		sub_data->SetNumberOfTuples(data->GetNumberOfTuples());
		sub_data->CopyComponent(0, data, idx);

		grid->GetPointData()->SetScalars(sub_data);

		vtk_setup(grid, config);
	}
}

#endif
