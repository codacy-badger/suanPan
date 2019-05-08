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

#include <suanPan>

using std::ifstream;
using std::string;
using std::vector;

int process_command(const shared_ptr<Bead>& model, istringstream& command) {
	if(model == nullptr) return SUANPAN_SUCCESS;

	string command_id;
	if(!get_input(command, command_id)) return SUANPAN_SUCCESS;

	if(is_equal(command_id, "exit")) return SUANPAN_EXIT;
	if(is_equal(command_id, "quit")) return SUANPAN_EXIT;

	if(is_equal(command_id, "file")) return process_file(model, command);
	if(is_equal(command_id, "load")) return process_file(model, command);

	if(is_equal(command_id, "domain")) return create_new_domain(model, command);

	if(is_equal(command_id, "enable")) return enable_object(model, command);
	if(is_equal(command_id, "disable")) return disable_object(model, command);
	if(is_equal(command_id, "mute")) return disable_object(model, command);
	if(is_equal(command_id, "erase")) return erase_object(model, command);
	if(is_equal(command_id, "delete")) return erase_object(model, command);
	if(is_equal(command_id, "remove")) return erase_object(model, command);
	if(is_equal(command_id, "save")) return save_object(model, command);
	if(is_equal(command_id, "list")) return list_object(model, command);

	const auto& domain = get_current_domain(model);

	if(is_equal(command_id, "acceleration")) return create_new_acceleration(domain, command);
	if(is_equal(command_id, "amplitude")) return create_new_amplitude(domain, command);
	if(is_equal(command_id, "cload")) return create_new_cload(domain, command);
	if(is_equal(command_id, "converger")) return create_new_converger(domain, command);
	if(is_equal(command_id, "criterion")) return create_new_criterion(domain, command);
	if(is_equal(command_id, "disp")) return create_new_displacement(domain, command);
	if(is_equal(command_id, "displacement")) return create_new_displacement(domain, command);
	if(is_equal(command_id, "dispload")) return create_new_displacement(domain, command);
	if(is_equal(command_id, "element")) return create_new_element(domain, command);
	if(is_equal(command_id, "fix")) return create_new_bc(domain, command);
	if(is_equal(command_id, "fix2")) return create_new_bc2(domain, command);
	if(is_equal(command_id, "import")) return create_new_external_module(domain, command);
	if(is_equal(command_id, "initial")) return create_new_initial(domain, command);
	if(is_equal(command_id, "integrator")) return create_new_integrator(domain, command);
	if(is_equal(command_id, "material")) return create_new_material(domain, command);
	if(is_equal(command_id, "modifier")) return create_new_modifier(domain, command);
	if(is_equal(command_id, "mass")) return create_new_mass(domain, command);
	if(is_equal(command_id, "node")) return create_new_node(domain, command);
	if(is_equal(command_id, "orientation")) return create_new_orientation(domain, command);
	if(is_equal(command_id, "recorder")) return create_new_recorder(domain, command);
	if(is_equal(command_id, "plainrecorder")) return create_new_plainrecorder(domain, command);
	if(is_equal(command_id, "hdf5recorder")) return create_new_hdf5recorder(domain, command);
	if(is_equal(command_id, "section")) return create_new_section(domain, command);
	if(is_equal(command_id, "solver")) return create_new_solver(domain, command);
	if(is_equal(command_id, "step")) return create_new_step(domain, command);

	if(is_equal(command_id, "set")) return set_property(domain, command);

	if(is_equal(command_id, "materialtest1d")) return test_material1d(domain, command);
	if(is_equal(command_id, "materialtest2d")) return test_material2d(domain, command);
	if(is_equal(command_id, "materialtest3d")) return test_material3d(domain, command);

	if(is_equal(command_id, "plot")) return vtk_parser(domain, command);

	if(is_equal(command_id, "peek")) return print_info(domain, command);

	if(is_equal(command_id, "analyze")) return model->analyze();
	if(is_equal(command_id, "analyse")) return model->analyze();

	if(is_equal(command_id, "clear")) {
		domain->clear_status();
		return SUANPAN_SUCCESS;
	}

	if(is_equal(command_id, "reset")) {
		domain->reset_status();
		return SUANPAN_SUCCESS;
	}

	if(is_equal(command_id, "summary")) {
		domain->summary();
		return SUANPAN_SUCCESS;
	}

	if(is_equal(command_id, "version")) print_version();
	else suanpan_info("command not found.\n");

	return SUANPAN_SUCCESS;
}

int process_file(const shared_ptr<Bead>& model, const char* file_name) {
	ifstream input_file(file_name);

	if(!input_file.is_open()) {
		string new_name = file_name;
		new_name += ".supan";
		input_file.open(new_name);
		if(!input_file.is_open()) {
			new_name = file_name;
			std::transform(new_name.begin(), new_name.end(), new_name.begin(), suanpan::to_upper);
			new_name += ".supan";
			input_file.open(new_name);
			if(!input_file.is_open()) {
				new_name = file_name;
				std::transform(new_name.begin(), new_name.end(), new_name.begin(), suanpan::to_lower);
				new_name += ".supan";
				input_file.open(new_name);
				if(!input_file.is_open()) {
					suanpan_error("process_file() cannot open the input file.\n");
					return SUANPAN_EXIT;
				}
			}
		}
	}

	string all_line, command_line;
	while(!getline(input_file, command_line).fail())
		if(!command_line.empty() && command_line[0] != '#' && command_line[0] != '!') {
			while(*command_line.crbegin() == ' ' || *command_line.crbegin() == '\t') command_line.pop_back();
			if(*command_line.crbegin() == '\\') {
				command_line.pop_back();
				all_line.append(command_line);
			} else {
				all_line.append(command_line);
				const auto if_comment = all_line.find('!');
				if(string::npos != if_comment) all_line.erase(if_comment);
				istringstream tmp_str(all_line);
				if(process_command(model, tmp_str) == SUANPAN_EXIT) return SUANPAN_EXIT;
				all_line.clear();
			}
		}

	return SUANPAN_SUCCESS;
}

int process_file(const shared_ptr<Bead>& model, istringstream& command) {
	string file_name;
	if(!get_input(command, file_name)) {
		suanpan_info("process_file() needs a file name.\n");
		return SUANPAN_SUCCESS;
	}

	return process_file(model, file_name.c_str());
}

int create_new_domain(const shared_ptr<Bead>& model, istringstream& command) {
	unsigned domain_id;
	if(!get_input(command, domain_id)) {
		suanpan_info("create_new_domain() requires a tag.\n");
		return SUANPAN_SUCCESS;
	}

	model->set_current_domain_tag(domain_id);

	auto& tmp_domain = get_domain(model, domain_id);

	if(tmp_domain == nullptr) {
		tmp_domain = make_shared<Domain>(domain_id);
		if(tmp_domain != nullptr) suanpan_info("create_new_domain() successfully creates Domain %u.\n", domain_id);
	} else suanpan_info("create_new_domain() switches to Domain %u.\n", domain_id);

	return SUANPAN_SUCCESS;
}

int disable_object(const shared_ptr<Bead>& model, istringstream& command) {
	const auto& domain = get_current_domain(model);
	if(domain == nullptr) {
		suanpan_info("disable_object() needs a valid domain.\n");
		return SUANPAN_SUCCESS;
	}

	string object_type;
	if(!get_input(command, object_type)) {
		suanpan_info("disable_object() needs object type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(is_equal(object_type, "domain")) while(get_input(command, tag)) model->disable_domain(tag);
	else if(is_equal(object_type, "step")) while(get_input(command, tag)) domain->disable_step(tag);
	else if(is_equal(object_type, "converger")) while(get_input(command, tag)) domain->disable_converger(tag);
	else if(is_equal(object_type, "bc")) while(get_input(command, tag)) domain->disable_constraint(tag);
	else if(is_equal(object_type, "constraint")) while(get_input(command, tag)) domain->disable_constraint(tag);
	else if(is_equal(object_type, "element")) while(get_input(command, tag)) domain->disable_element(tag);
	else if(is_equal(object_type, "load")) while(get_input(command, tag)) domain->disable_load(tag);
	else if(is_equal(object_type, "material")) while(get_input(command, tag)) domain->disable_material(tag);
	else if(is_equal(object_type, "node")) while(get_input(command, tag)) domain->disable_node(tag);
	else if(is_equal(object_type, "recorder")) while(get_input(command, tag)) domain->disable_recorder(tag);

	return SUANPAN_SUCCESS;
}

int enable_object(const shared_ptr<Bead>& model, istringstream& command) {
	const auto& domain = get_current_domain(model);
	if(domain == nullptr) {
		suanpan_info("enable_object() needs a valid domain.\n");
		return SUANPAN_SUCCESS;
	}

	string object_type;
	if(!get_input(command, object_type)) {
		suanpan_info("enable_object() needs object type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(is_equal(object_type, "domain")) while(get_input(command, tag)) model->enable_domain(tag);
	else if(is_equal(object_type, "step")) while(get_input(command, tag)) domain->enable_step(tag);
	else if(is_equal(object_type, "converger")) while(get_input(command, tag)) domain->enable_converger(tag);
	else if(is_equal(object_type, "bc")) while(get_input(command, tag)) domain->enable_constraint(tag);
	else if(is_equal(object_type, "constraint")) while(get_input(command, tag)) domain->enable_constraint(tag);
	else if(is_equal(object_type, "element")) while(get_input(command, tag)) domain->enable_element(tag);
	else if(is_equal(object_type, "load")) while(get_input(command, tag)) domain->enable_load(tag);
	else if(is_equal(object_type, "material")) while(get_input(command, tag)) domain->enable_material(tag);
	else if(is_equal(object_type, "node")) while(get_input(command, tag)) domain->enable_node(tag);
	else if(is_equal(object_type, "recorder")) while(get_input(command, tag)) domain->enable_recorder(tag);

	return SUANPAN_SUCCESS;
}

int erase_object(const shared_ptr<Bead>& model, istringstream& command) {
	const auto& domain = get_current_domain(model);
	if(domain == nullptr) {
		suanpan_info("erase_object() needs a valid domain.\n");
		return SUANPAN_SUCCESS;
	}

	string object_type;
	if(!get_input(command, object_type)) {
		suanpan_info("erase_object() needs object type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(is_equal(object_type, "domain")) while(get_input(command, tag)) model->erase_domain(tag);
	else if(is_equal(object_type, "step")) while(get_input(command, tag)) domain->erase_step(tag);
	else if(is_equal(object_type, "converger")) while(get_input(command, tag)) domain->erase_converger(tag);
	else if(is_equal(object_type, "bc")) while(get_input(command, tag)) domain->erase_constraint(tag);
	else if(is_equal(object_type, "constraint")) while(get_input(command, tag)) domain->erase_constraint(tag);
	else if(is_equal(object_type, "element")) while(get_input(command, tag)) domain->erase_element(tag);
	else if(is_equal(object_type, "load")) while(get_input(command, tag)) domain->erase_load(tag);
	else if(is_equal(object_type, "material")) while(get_input(command, tag)) domain->erase_material(tag);
	else if(is_equal(object_type, "node")) while(get_input(command, tag)) domain->erase_node(tag);
	else if(is_equal(object_type, "recorder")) while(get_input(command, tag)) domain->erase_recorder(tag);

	return SUANPAN_SUCCESS;
}

int save_object(const shared_ptr<Bead>& model, istringstream& command) {
	const auto& domain = get_current_domain(model);
	if(domain == nullptr) {
		suanpan_info("erase_object() needs a valid domain.\n");
		return SUANPAN_SUCCESS;
	}

	string object_id;
	if(!get_input(command, object_id)) {
		suanpan_info("save_object() needs a valid object type.\n");
		return SUANPAN_SUCCESS;
	}

	if(is_equal(object_id, "Recorder")) {
		unsigned tag;
		while(get_input(command, tag)) if(domain->find_recorder(tag)) domain->get_recorder(tag)->save();
	} else if(is_equal(object_id, "Stiffness")) {
		string name;
		if(!command.eof() && !get_input(command, name)) name = "K";
		domain->get_factory()->get_stiffness()->save(name.c_str());
	} else if(is_equal(object_id, "Mass")) {
		string name;
		if(!command.eof() && !get_input(command, name)) name = "M";
		domain->get_factory()->get_mass()->save(name.c_str());
	} else if(is_equal(object_id, "Damping")) {
		string name;
		if(!command.eof() && !get_input(command, name)) name = "C";
		domain->get_factory()->get_damping()->save(name.c_str());
	}

	return SUANPAN_SUCCESS;
}

int list_object(const shared_ptr<Bead>& model, istringstream& command) {
	const auto& domain = get_current_domain(model);
	if(domain == nullptr) {
		suanpan_info("list_object() needs a valid domain.\n");
		return SUANPAN_SUCCESS;
	}

	string object_type;
	if(!get_input(command, object_type)) {
		suanpan_info("list_object() needs object type.\n");
		return SUANPAN_SUCCESS;
	}

	vector<unsigned> list;
	if(is_equal(object_type, "converger")) for(const auto& I : domain->get_converger_pool()) list.emplace_back(I->get_tag());
	else if(is_equal(object_type, "constraint")) for(const auto& I : domain->get_constraint_pool()) list.emplace_back(I->get_tag());
	else if(is_equal(object_type, "element")) for(const auto& I : domain->get_element_pool()) list.emplace_back(I->get_tag());
	else if(is_equal(object_type, "load")) for(const auto& I : domain->get_load_pool()) list.emplace_back(I->get_tag());
	else if(is_equal(object_type, "material")) for(const auto& I : domain->get_material_pool()) list.emplace_back(I->get_tag());
	else if(is_equal(object_type, "node")) for(const auto& I : domain->get_node_pool()) list.emplace_back(I->get_tag());
	else if(is_equal(object_type, "recorder")) for(const auto& I : domain->get_recorder_pool()) list.emplace_back(I->get_tag());

	suanpan_info("This domain has the following %ss:", object_type.c_str());
	for(const auto& I : list) suanpan_info("\t%u", I);
	suanpan_info(".\n");

	return SUANPAN_SUCCESS;
}

int create_new_acceleration(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned load_id;
	if(!get_input(command, load_id)) {
		suanpan_info("create_new_acceleration() needs a tag.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned amplitude_id;
	if(!get_input(command, amplitude_id)) {
		suanpan_info("create_new_acceleration() needs a valid amplitude tag.\n");
		return SUANPAN_SUCCESS;
	}

	double magnitude;
	if(!get_input(command, magnitude)) {
		suanpan_info("create_new_acceleration() needs load magnitude.\n");
		return SUANPAN_SUCCESS;
	}

	uword dof_id;
	vector<uword> dof_pool;
	while(get_input(command, dof_id)) dof_pool.emplace_back(dof_id);

	const auto& step_tag = domain->get_current_step_tag();

	if(!domain->insert(make_shared<Acceleration>(load_id, step_tag, magnitude, uvec(dof_pool), amplitude_id))) suanpan_error("create_new_acceleration() fails to create new load.\n");

	return SUANPAN_SUCCESS;
}

int create_new_amplitude(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string amplitude_type;
	if(!get_input(command, amplitude_type)) {
		suanpan_info("create_new_amplitude() needs a valid amplitude type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_info("create_new_amplitude() needs a valid amplitude type.\n");
		return SUANPAN_SUCCESS;
	}

	const auto step_tag = domain->get_current_step_tag();

	if(is_equal(amplitude_type, "Constant")) domain->insert(make_shared<Constant>(tag, step_tag));
	else if(is_equal(amplitude_type, "Tabular")) {
		string file_name;
		if(!get_input(command, file_name)) {
			suanpan_info("create_new_amplitude() needs a valid file.\n");
			return SUANPAN_SUCCESS;
		}
		domain->insert(make_shared<Tabular>(tag, file_name.c_str(), step_tag));
	} else if(is_equal(amplitude_type, "Decay")) {
		double A;
		if(!get_input(command, A)) {
			suanpan_info("create_new_amplitude() needs a A.\n");
			return SUANPAN_SUCCESS;
		}
		double TD;
		if(!get_input(command, TD)) {
			suanpan_info("create_new_amplitude() needs a TD.\n");
			return SUANPAN_SUCCESS;
		}
		domain->insert(make_shared<Decay>(tag, A, TD, step_tag));
	} else if(is_equal(amplitude_type, "Linear")) {
		double A;
		if(!get_input(command, A)) {
			suanpan_info("create_new_amplitude() needs a slope.\n");
			return SUANPAN_SUCCESS;
		}
		domain->insert(make_shared<Linear>(tag, A, step_tag));
	} else if(is_equal(amplitude_type, "Combine")) {
		vector<uword> tag_pool;
		uword t_tag;
		while(get_input(command, t_tag)) tag_pool.emplace_back(t_tag);
		domain->insert(make_shared<Combine>(tag, uvec(tag_pool), step_tag));
	} else if(is_equal(amplitude_type, "Modulated") || is_equal(amplitude_type, "Sine") || is_equal(amplitude_type, "Cosine")) {
		double W;
		if(!get_input(command, W)) {
			suanpan_info("create_new_amplitude() needs a period/amplitude.\n");
			return SUANPAN_SUCCESS;
		}

		double amp;
		vector<double> A;
		while(get_input(command, amp)) A.emplace_back(amp);

		if(is_equal(amplitude_type, "Modulated")) domain->insert(make_shared<Modulated>(tag, W, std::move(A), step_tag));
		else if(is_equal(amplitude_type, "Sine")) domain->insert(make_shared<Sine>(tag, W, std::move(A), step_tag));
		else if(is_equal(amplitude_type, "Cosine")) domain->insert(make_shared<Cosine>(tag, W, std::move(A), step_tag));
	}

	return SUANPAN_SUCCESS;
}

int create_new_bc(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned bc_id;
	if(!get_input(command, bc_id)) {
		suanpan_info("create_new_bc() needs BC type.\n");
		return SUANPAN_SUCCESS;
	}

	string dof_id;
	if(!get_input(command, dof_id)) {
		suanpan_info("create_new_bc() needs valid DoFs.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned node;
	vector<uword> node_tag;
	while(get_input(command, node)) node_tag.push_back(node);

	const auto bc_type = suanpan::to_lower(dof_id[0]);

	const auto& step_tag = domain->get_current_step_tag();

	if(is_equal(bc_type, 'p')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), "PINNED"));
	else if(is_equal(bc_type, 'e')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), "ENCASTRE"));
	else if(is_equal(bc_type, 'x')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), "XSYMM"));
	else if(is_equal(bc_type, 'y')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), "YSYMM"));
	else if(is_equal(bc_type, 'z')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), "ZSYMM"));
	else if(is_equal(bc_type, '1')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 1));
	else if(is_equal(bc_type, '2')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 2));
	else if(is_equal(bc_type, '3')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 3));
	else if(is_equal(bc_type, '4')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 4));
	else if(is_equal(bc_type, '5')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 5));
	else if(is_equal(bc_type, '6')) domain->insert(make_shared<BC>(bc_id, step_tag, uvec(node_tag), 6));

	return SUANPAN_SUCCESS;
}

int create_new_bc2(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned bc_id;
	if(!get_input(command, bc_id)) {
		suanpan_info("create_new_bc2() needs BC type.\n");
		return SUANPAN_SUCCESS;
	}

	string dof_id;
	if(!get_input(command, dof_id)) {
		suanpan_info("create_new_bc2() needs valid DoFs.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned node;
	vector<uword> node_tag;
	while(get_input(command, node)) node_tag.push_back(node);

	const auto bc_type = suanpan::to_lower(dof_id[0]);

	const auto& step_tag = domain->get_current_step_tag();

	if(is_equal(bc_type, 'p')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), "PINNED"));
	else if(is_equal(bc_type, 'e')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), "ENCASTRE"));
	else if(is_equal(bc_type, 'x')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), "XSYMM"));
	else if(is_equal(bc_type, 'y')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), "YSYMM"));
	else if(is_equal(bc_type, 'z')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), "ZSYMM"));
	else if(is_equal(bc_type, '1')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), 1));
	else if(is_equal(bc_type, '2')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), 2));
	else if(is_equal(bc_type, '3')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), 3));
	else if(is_equal(bc_type, '4')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), 4));
	else if(is_equal(bc_type, '5')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), 5));
	else if(is_equal(bc_type, '6')) domain->insert(make_shared<BC2>(bc_id, step_tag, uvec(node_tag), 6));

	return SUANPAN_SUCCESS;
}

int create_new_cload(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned load_id;
	if(!get_input(command, load_id)) {
		suanpan_info("create_new_cload() needs a tag.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned amplitude_id;
	if(!get_input(command, amplitude_id)) {
		suanpan_info("create_new_cload() needs a valid amplitude tag.\n");
		return SUANPAN_SUCCESS;
	}

	double magnitude;
	if(!get_input(command, magnitude)) {
		suanpan_info("create_new_cload() needs load magnitude.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned dof_id;
	if(!get_input(command, dof_id)) {
		suanpan_info("create_new_cload() needs a valid DoF.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned node;
	vector<uword> node_tag;
	while(get_input(command, node)) node_tag.push_back(node);

	if(!domain->insert(make_shared<NodalLoad>(load_id, domain->get_current_step_tag(), magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_cload() fails to create new load.\n");

	return SUANPAN_SUCCESS;
}

int create_new_converger(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string converger_id;
	if(!get_input(command, converger_id)) {
		suanpan_info("create_new_converger() requires converger type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_info("create_new_converger() requires a tag.\n");
		return SUANPAN_SUCCESS;
	}

	auto tolerance = 1E-6;
	if(!command.eof() && !get_input(command, tolerance)) {
		suanpan_info("create_new_converger() reads wrong tolerance.\n");
		return SUANPAN_SUCCESS;
	}

	auto max_iteration = 10;
	if(!command.eof() && !get_input(command, max_iteration)) {
		suanpan_info("create_new_converger() reads wrong max iteration.\n");
		return SUANPAN_SUCCESS;
	}

	string print_flag = "false";
	if(!command.eof() && !get_input(command, print_flag)) {
		suanpan_info("create_new_converger() reads wrong print flag.\n");
		return SUANPAN_SUCCESS;
	}

	auto code = 0;
	if(is_equal(converger_id, "AbsResidual") && domain->insert(make_shared<AbsResidual>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
	else if(is_equal(converger_id, "RelResidual") && domain->insert(make_shared<RelResidual>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
	else if(is_equal(converger_id, "AbsIncreDisp") && domain->insert(make_shared<AbsIncreDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
	else if(is_equal(converger_id, "RelIncreDisp") && domain->insert(make_shared<RelIncreDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
	else if(is_equal(converger_id, "AbsDisp") && domain->insert(make_shared<AbsDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
	else if(is_equal(converger_id, "RelDisp") && domain->insert(make_shared<RelDisp>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
	else if(is_equal(converger_id, "AbsError") && domain->insert(make_shared<AbsError>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
	else if(is_equal(converger_id, "RelError") && domain->insert(make_shared<RelError>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
	else if(is_equal(converger_id, "AbsIncreEnergy") && domain->insert(make_shared<AbsIncreEnergy>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
	else if(is_equal(converger_id, "RelIncreEnergy") && domain->insert(make_shared<RelIncreEnergy>(tag, tolerance, max_iteration, is_true(print_flag)))) code = 1;
	else if(is_equal(converger_id, "FixedNumber") && domain->insert(make_shared<FixedNumber>(tag, max_iteration, is_true(print_flag)))) code = 1;
	else suanpan_info("create_new_converger() cannot identify the converger type.\n");

	if(code == 1) {
		if(domain->get_current_step_tag() != 0) domain->get_current_step()->set_converger_tag(tag);
		domain->set_current_converger_tag(tag);
	} else suanpan_info("create_new_converger() fails to create the new converger.\n");

	return SUANPAN_SUCCESS;
}

int create_new_criterion(const shared_ptr<DomainBase>& domain, istringstream& command) {
	const auto& step_tag = domain->get_current_step_tag();
	if(step_tag == 0) {
		suanpan_info("create_new_criterion() needs a valid step.\n");
		return SUANPAN_SUCCESS;
	}

	string criterion_type;
	if(!get_input(command, criterion_type)) {
		suanpan_info("create_new_criterion() need a criterion type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_info("create_new_criterion() requires a tag.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned node;
	if(!get_input(command, node)) {
		suanpan_info("create_new_criterion() requires a node.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned dof;
	if(!get_input(command, dof)) {
		suanpan_info("create_new_criterion() requires a dof.\n");
		return SUANPAN_SUCCESS;
	}

	double limit;
	if(!get_input(command, limit)) {
		suanpan_info("create_new_criterion() requires a limit.\n");
		return SUANPAN_SUCCESS;
	}

	if(is_equal(criterion_type, "MaxDisplacement")) domain->insert(make_shared<MaxDisplacement>(tag, step_tag, node, dof, limit));
	else if(is_equal(criterion_type, "MinDisplacement")) domain->insert(make_shared<MinDisplacement>(tag, step_tag, node, dof, limit));

	return SUANPAN_SUCCESS;
}

int create_new_displacement(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned load_id;
	if(!get_input(command, load_id)) {
		suanpan_info("create_new_displacement() needs a tag.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned amplitude_id;
	if(!get_input(command, amplitude_id)) {
		suanpan_info("create_new_displacement() needs a valid amplitude tag.\n");
		return SUANPAN_SUCCESS;
	}

	double magnitude;
	if(!get_input(command, magnitude)) {
		suanpan_info("create_new_displacement() needs load magnitude.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned dof_id;
	if(!get_input(command, dof_id)) {
		suanpan_info("create_new_displacement() needs a valid DoF.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned node;
	vector<uword> node_tag;
	while(get_input(command, node)) node_tag.push_back(node);

	const auto& step_tag = domain->get_current_step_tag();

	if(!domain->insert(make_shared<NodalDisplacement>(load_id, step_tag, magnitude, uvec(node_tag), dof_id, amplitude_id))) suanpan_error("create_new_displacement() fails to create new load.\n");

	return SUANPAN_SUCCESS;
}

int create_new_external_module(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string library_name;

	if(!get_input(command, library_name)) {
		suanpan_info("create_new_external_module() needs module name.\n");
		return SUANPAN_SUCCESS;
	}

	auto code = 0;
	for(const auto& I : domain->get_external_module_pool())
		if(is_equal(I->library_name, library_name) || I->locate_module(library_name)) {
			code = 1;
			break;
		}

	if(code == 0) domain->insert(make_shared<ExternalModule>(library_name));

	return SUANPAN_SUCCESS;
}

int create_new_initial(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string variable_type;
	if(!get_input(command, variable_type)) {
		suanpan_error("create_new_initial() needs a valid variable type.\n");
		return SUANPAN_SUCCESS;
	}

	double magnitude;
	if(!get_input(command, magnitude)) {
		suanpan_error("create_new_initial() needs a valid magnitude.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned dof_tag;
	if(!get_input(command, dof_tag)) {
		suanpan_error("create_new_initial() needs a valid dof tag.\n");
		return SUANPAN_SUCCESS;
	}

	if(is_equal("displacement", variable_type))
		while(!command.eof()) {
			unsigned node_tag;
			if(get_input(command, node_tag) && domain->find_node(node_tag)) {
				auto t_variable = domain->get_node(node_tag)->get_current_displacement();
				if(t_variable.n_elem < dof_tag) t_variable.resize(dof_tag);
				t_variable(dof_tag - 1) = magnitude;
				domain->get_node(node_tag)->set_current_displacement(t_variable);
			}
		}
	else if(is_equal("velocity", variable_type))
		while(!command.eof()) {
			unsigned node_tag;
			if(get_input(command, node_tag) && domain->find_node(node_tag)) {
				auto t_variable = domain->get_node(node_tag)->get_current_velocity();
				if(t_variable.n_elem < dof_tag) t_variable.resize(dof_tag);
				t_variable(dof_tag - 1) = magnitude;
				domain->get_node(node_tag)->set_current_velocity(t_variable);
			}
		}
	else if(is_equal("acceleration", variable_type))
		while(!command.eof()) {
			unsigned node_tag;
			if(get_input(command, node_tag) && domain->find_node(node_tag)) {
				auto t_variable = domain->get_node(node_tag)->get_current_acceleration();
				if(t_variable.n_elem < dof_tag) t_variable.resize(dof_tag);
				t_variable(dof_tag - 1) = magnitude;
				domain->get_node(node_tag)->set_current_acceleration(t_variable);
			}
		}

	return SUANPAN_SUCCESS;
}

int create_new_integrator(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string integrator_type;
	if(!get_input(command, integrator_type)) {
		suanpan_error("create_new_integrator() needs a valid integrator type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("create_new_integrator() needs a valid tag.\n");
		return SUANPAN_SUCCESS;
	}

	auto code = 0;
	if(is_equal(integrator_type, "Newmark")) {
		auto alpha = .25, beta = .5;
		if(!command.eof()) {
			if(!get_input(command, alpha)) {
				suanpan_error("create_new_integrator() needs a valid alpha.\n");
				return SUANPAN_SUCCESS;
			}
			if(!get_input(command, beta)) {
				suanpan_error("create_new_integrator() needs a valid beta.\n");
				return SUANPAN_SUCCESS;
			}
		}
		if(domain->insert(make_shared<Newmark>(tag, alpha, beta))) code = 1;
	} else if(is_equal(integrator_type, "RayleighNewmark")) {
		double damping_alpha, damping_beta;
		if(!get_input(command, damping_alpha)) {
			suanpan_error("create_new_integrator() needs a valid alpha for Rayleigh damping.\n");
			return SUANPAN_SUCCESS;
		}
		if(!get_input(command, damping_beta)) {
			suanpan_error("create_new_integrator() needs a valid beta for Rayleigh damping.\n");
			return SUANPAN_SUCCESS;
		}

		auto alpha = .25, beta = .5;
		if(!command.eof()) {
			if(!get_input(command, alpha)) {
				suanpan_error("create_new_integrator() needs a valid alpha.\n");
				return SUANPAN_SUCCESS;
			}
			if(!get_input(command, beta)) {
				suanpan_error("create_new_integrator() needs a valid beta.\n");
				return SUANPAN_SUCCESS;
			}
		}
		if(domain->insert(make_shared<RayleighNewmark>(tag, damping_alpha, damping_beta, alpha, beta))) code = 1;
	} else if(is_equal(integrator_type, "GeneralizedAlpha")) {
		auto alpha_m = .0, alpha_f = .0;
		if(!command.eof()) {
			if(!get_input(command, alpha_m)) {
				suanpan_error("create_new_integrator() needs a valid alpha_m.\n");
				return SUANPAN_SUCCESS;
			}
			if(!get_input(command, alpha_f)) {
				suanpan_error("create_new_integrator() needs a valid alpha_f.\n");
				return SUANPAN_SUCCESS;
			}
		}
		if(domain->insert(make_shared<GeneralizedAlpha>(tag, alpha_m, alpha_f))) code = 1;
	} else if(is_equal(integrator_type, "CentralDifference")) if(domain->insert(make_shared<CentralDifference>(tag))) code = 1;

	if(code == 1) {
		if(domain->get_current_step_tag() != 0) domain->get_current_step()->set_integrator_tag(tag);
		domain->set_current_integrator_tag(tag);
	} else suanpan_info("create_new_integrator() fails to create the new integrator.\n");

	return SUANPAN_SUCCESS;
}

int create_new_mass(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_debug("create_new_mass() needs a valid tag.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned node;
	if(!get_input(command, node)) {
		suanpan_debug("create_new_mass() needs one valid node.\n");
		return SUANPAN_SUCCESS;
	}

	double magnitude;
	if(!get_input(command, magnitude)) {
		suanpan_debug("create_new_mass() needs a valid magnitude.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned dof;
	vector<uword> dof_tag;
	while(get_input(command, dof)) dof_tag.push_back(dof);

	domain->insert(make_shared<Mass>(tag, node, magnitude, uvec(dof_tag)));

	return SUANPAN_SUCCESS;
}

int create_new_modifier(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string modifier_type;
	if(!get_input(command, modifier_type)) {
		suanpan_error("create_new_modifier() needs a valid modifier type.\n");
		return SUANPAN_SUCCESS;
	}

	unique_ptr<Modifier> new_modifier = nullptr;

	if(is_equal(modifier_type, "LumpedSimple")) {
		unsigned tag;
		if(!get_input(command, tag)) {
			suanpan_info("create_new_modifier() needs a valid tag.\n");
			return SUANPAN_SUCCESS;
		}

		vector<uword> element_tag;
		unsigned e_tag;
		while(!command.eof()) if(get_input(command, e_tag)) element_tag.emplace_back(e_tag);

		new_modifier = make_unique<LumpedSimple>(tag, uvec(element_tag));
	} else if(is_equal(modifier_type, "LumpedScale")) {
		unsigned tag;
		if(!get_input(command, tag)) {
			suanpan_info("create_new_modifier() needs a valid tag.\n");
			return SUANPAN_SUCCESS;
		}

		vector<uword> element_tag;
		unsigned e_tag;
		while(!command.eof()) if(get_input(command, e_tag)) element_tag.emplace_back(e_tag);

		new_modifier = make_unique<LumpedScale>(tag, uvec(element_tag));
	} else if(is_equal(modifier_type, "Rayleigh")) {
		unsigned tag;
		if(!get_input(command, tag)) {
			suanpan_info("create_new_modifier() needs a valid tag.\n");
			return SUANPAN_SUCCESS;
		}

		double a, b, c, d;
		if(!get_input(command, a)) {
			suanpan_info("create_new_modifier() needs four valid numbers.\n");
			return SUANPAN_SUCCESS;
		}
		if(!get_input(command, b)) {
			suanpan_info("create_new_modifier() needs four valid numbers.\n");
			return SUANPAN_SUCCESS;
		}
		if(!get_input(command, c)) {
			suanpan_info("create_new_modifier() needs four valid numbers.\n");
			return SUANPAN_SUCCESS;
		}
		if(!get_input(command, d)) {
			suanpan_info("create_new_modifier() needs four valid numbers.\n");
			return SUANPAN_SUCCESS;
		}
		vector<uword> element_tag;
		unsigned e_tag;
		while(!command.eof()) if(get_input(command, e_tag)) element_tag.emplace_back(e_tag);

		new_modifier = make_unique<Rayleigh>(tag, a, b, c, d, uvec(element_tag));
	} else {
		// check if the library is already loaded
		auto code = false;
		for(const auto& I : domain->get_external_module_pool())
			if(is_equal(I->library_name, modifier_type) || I->locate_module(modifier_type)) {
				code = true;
				break;
			}

		// not loaded then try load it
		if(!code && domain->insert(make_shared<ExternalModule>(modifier_type))) code = true;

		// if loaded find corresponding function
		if(code)
			for(const auto& I : domain->get_external_module_pool()) {
				if(I->locate_module(modifier_type)) I->new_object(new_modifier, command);
				if(new_modifier != nullptr) break;
			}
	}

	if(new_modifier == nullptr || !domain->insert(std::move(new_modifier))) suanpan_error("create_new_modifier() fails to create new modifier.\n");

	return SUANPAN_SUCCESS;
}

int create_new_node(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned node_id;
	if(!get_input(command, node_id)) {
		suanpan_info("create_new_node() needs a tag.\n");
		return SUANPAN_SUCCESS;
	}

	vector<double> coor;
	double X;
	while(get_input(command, X)) coor.push_back(X);

	if(!domain->insert(make_shared<Node>(node_id, vec(coor)))) suanpan_debug("create_new_node() fails to insert Node %u.\n", node_id);

	return SUANPAN_SUCCESS;
}

int create_new_orientation(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string file_type;
	if(!get_input(command, file_type)) {
		suanpan_info("create_new_orientation() needs a valid type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_info("create_new_orientation() needs a valid tag.\n");
		return SUANPAN_SUCCESS;
	}

	vec xyz(3);
	for(auto I = 0; I < 3; ++I) get_input(command, xyz(I));

	if(is_equal(file_type, "B3DL")) domain->insert(make_shared<B3DL>(tag, std::move(xyz)));

	return SUANPAN_SUCCESS;
}

int create_new_recorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_info("create_new_recorder() needs a valid tag.\n");
		return SUANPAN_SUCCESS;
	}

	string file_type;
	if(!get_input(command, file_type)) {
		suanpan_info("create_new_recorder() needs a valid object type.\n");
		return SUANPAN_SUCCESS;
	}

	string object_type;
	if(!get_input(command, object_type)) {
		suanpan_info("create_new_recorder() needs a valid object type.\n");
		return SUANPAN_SUCCESS;
	}

	if(is_equal(object_type, "Eigen")) {
		if(!domain->insert(make_shared<EigenRecorder>(tag, is_equal(file_type[0], 'h')))) suanpan_info("create_new_recorder() fails to create a new eigen recorder.\n");
		return SUANPAN_SUCCESS;
	}

	string variable_type;
	if(!get_input(command, variable_type)) {
		suanpan_info("create_new_recorder() needs a valid recorder type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned inverval = 1;

	while(true) {
		if(command.peek() == EOF) return SUANPAN_SUCCESS;
		if(command.peek() == '\t' || command.peek() == ' ') command.ignore();
		else break;
	}

	if(is_equal(command.peek(), 'e') || is_equal(command.peek(), 'i')) {
		string tmp_string;
		get_input(command, tmp_string);
		if(!get_input(command, inverval)) return SUANPAN_SUCCESS;
	}

	if(is_equal(object_type, "Frame")) {
		if(!domain->insert(make_shared<FrameRecorder>(tag, to_list(variable_type.c_str()), inverval))) suanpan_info("create_new_recorder() fails to create a new eigen recorder.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned s_object_tag;
	vector<uword> object_tag;
	while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

	if(is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), inverval, true, is_equal(file_type[0], 'h')))) suanpan_info("create_new_recorder() fails to create a new node recorder.\n");
	else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), inverval, true, is_equal(file_type[0], 'h')))) suanpan_info("create_new_recorder() fails to create a new element recorder.\n");

	return SUANPAN_SUCCESS;
}

int create_new_plainrecorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_info("create_new_plainrecorder() needs a valid tag.\n");
		return SUANPAN_SUCCESS;
	}

	string object_type;
	if(!get_input(command, object_type)) {
		suanpan_info("create_new_plainrecorder() needs a valid object type.\n");
		return SUANPAN_SUCCESS;
	}

	if(is_equal(object_type, "Eigen")) {
		if(!domain->insert(make_shared<EigenRecorder>(tag, false))) suanpan_info("create_new_plainrecorder() fails to create a new eigen recorder.\n");
		return SUANPAN_SUCCESS;
	}

	string variable_type;
	if(!get_input(command, variable_type)) {
		suanpan_info("create_new_plainrecorder() needs a valid recorder type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned inverval = 1;

	while(true) {
		if(command.peek() == EOF) break;
		if(command.peek() != '\t' && command.peek() != ' ') break;
		command.ignore();
	}

	if(is_equal(command.peek(), 'e') || is_equal(command.peek(), 'i')) {
		string tmp_string;
		get_input(command, tmp_string);
		if(!get_input(command, inverval)) return SUANPAN_SUCCESS;
	}

	unsigned s_object_tag;
	vector<uword> object_tag;
	while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

	if(is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), inverval, true, false))) suanpan_info("create_new_plainrecorder() fails to create a new node recorder.\n");
	else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), inverval, true, false))) suanpan_info("create_new_plainrecorder() fails to create a new element recorder.\n");

	return SUANPAN_SUCCESS;
}

int create_new_hdf5recorder(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_info("create_new_hdf5recorder() needs a valid tag.\n");
		return SUANPAN_SUCCESS;
	}

	string object_type;
	if(!get_input(command, object_type)) {
		suanpan_info("create_new_hdf5recorder() needs a valid object type.\n");
		return SUANPAN_SUCCESS;
	}

	if(is_equal(object_type, "Eigen")) {
		if(!domain->insert(make_shared<EigenRecorder>(tag, true))) suanpan_info("create_new_hdf5recorder() fails to create a new eigen recorder.\n");
		return SUANPAN_SUCCESS;
	}

	string variable_type;
	if(!get_input(command, variable_type)) {
		suanpan_info("create_new_hdf5recorder() needs a valid recorder type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned inverval = 1;

	while(true) {
		if(command.peek() == EOF) break;
		if(command.peek() != '\t' && command.peek() != ' ') break;
		command.ignore();
	}

	if(is_equal(command.peek(), 'e') || is_equal(command.peek(), 'i')) {
		string tmp_string;
		get_input(command, tmp_string);
		if(!get_input(command, inverval)) return SUANPAN_SUCCESS;
	}

	if(is_equal(object_type, "Frame")) {
		if(!domain->insert(make_shared<FrameRecorder>(tag, to_list(variable_type.c_str()), inverval))) suanpan_info("create_new_recorder() fails to create a new eigen recorder.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned s_object_tag;
	vector<uword> object_tag;
	while(!command.eof() && get_input(command, s_object_tag)) object_tag.emplace_back(s_object_tag);

	if(is_equal(object_type, "Node") && !domain->insert(make_shared<NodeRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), inverval, true, true))) suanpan_info("create_new_hdf5recorder() fails to create a new node recorder.\n");
	else if(is_equal(object_type, "Element") && !domain->insert(make_shared<ElementRecorder>(tag, uvec(object_tag), to_list(variable_type.c_str()), inverval, true, true))) suanpan_info("create_new_hdf5recorder() fails to create a new element recorder.\n");

	return SUANPAN_SUCCESS;
}

int create_new_solver(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string solver_type;
	if(!get_input(command, solver_type)) {
		suanpan_info("create_new_solver() requires solver type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_info("create_new_solver() requires a tag.\n");
		return SUANPAN_SUCCESS;
	}

	auto code = 0;
	if(is_equal(solver_type, "Newton")) { if(domain->insert(make_shared<Newton>(tag))) code = 1; } else if(is_equal(solver_type, "modifiedNewton")) { if(domain->insert(make_shared<Newton>(tag, true))) code = 1; } else if(is_equal(solver_type, "mNewton")) { if(domain->insert(make_shared<Newton>(tag, true))) code = 1; } else if(is_equal(solver_type, "BFGS")) { if(domain->insert(make_shared<BFGS>(tag))) code = 1; } else if(is_equal(solver_type, "LBFGS")) {
		auto max_history = 20;
		if(!command.eof() && !get_input(command, max_history)) {
			suanpan_error("create_new_solver() requires a valid maximum step for LBFGS algorithm.\n");
			return SUANPAN_SUCCESS;
		}

		if(domain->insert(make_shared<BFGS>(tag, max_history))) code = 1;
	} else if(is_equal(solver_type, "Ramm")) {
		auto arc_length = .1;
		string fixed_arc_length = "False";

		if(!command.eof() && !get_input(command, arc_length)) {
			suanpan_error("create_new_solver() requires a valid arc length.\n");
			return SUANPAN_SUCCESS;
		}
		if(!command.eof() && !get_input(command, fixed_arc_length)) {
			suanpan_error("create_new_solver() requires a valid arc length switch.\n");
			return SUANPAN_SUCCESS;
		}

		if(domain->insert(make_shared<Ramm>(tag, arc_length, is_true(fixed_arc_length)))) code = 1;
	} else if(is_equal(solver_type, "DisplacementControl")) { if(domain->insert(make_shared<MPDC>(tag))) code = 1; } else if(is_equal(solver_type, "MPDC")) { if(domain->insert(make_shared<MPDC>(tag))) code = 1; } else suanpan_error("create_new_solver() cannot identify solver type.\n");

	if(code == 1) {
		if(domain->get_current_step_tag() != 0) domain->get_current_step()->set_solver_tag(tag);
		domain->set_current_solver_tag(tag);
	} else suanpan_error("create_new_solver() cannot create the new solver.\n");

	return SUANPAN_SUCCESS;
}

int create_new_step(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string step_type;
	if(!get_input(command, step_type)) {
		suanpan_info("create_new_step() requires step type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_info("create_new_step() requires a tag.\n");
		return SUANPAN_SUCCESS;
	}

	if(is_equal(step_type, "Frequency")) {
		auto eigen_number = 1;
		if(!command.eof() && !get_input(command, eigen_number)) {
			suanpan_info("create_new_step() reads a wrong number of eigenvalues.\n");
			return SUANPAN_SUCCESS;
		}
		if(domain->insert(make_shared<Frequency>(tag, eigen_number))) domain->set_current_step_tag(tag);
		else suanpan_error("create_new_step() cannot create the new step.\n");
	} else if(is_equal(step_type, "Buckling") || is_equal(step_type, "Buckle")) {
		if(domain->insert(make_shared<Buckle>(tag))) domain->set_current_step_tag(tag);
		else suanpan_error("create_new_step() cannot create the new step.\n");
	} else if(is_equal(step_type, "Static")) {
		auto time = 1.;
		if(!command.eof() && !get_input(command, time)) {
			suanpan_info("create_new_step() reads a wrong time period.\n");
			return SUANPAN_SUCCESS;
		}
		if(domain->insert(make_shared<Static>(tag, time))) domain->set_current_step_tag(tag);
		else suanpan_error("create_new_step() cannot create the new step.\n");
	} else if(is_equal(step_type, "Dynamic")) {
		auto time = 1.;
		if(!command.eof() && !get_input(command, time)) {
			suanpan_info("create_new_step() reads a wrong time period.\n");
			return SUANPAN_SUCCESS;
		}
		if(domain->insert(make_shared<Dynamic>(tag, time))) domain->set_current_step_tag(tag);
		else suanpan_error("create_new_step() cannot create the new step.\n");
	} else if(is_equal(step_type, "ArcLength")) {
		unsigned node;
		if(!get_input(command, node)) {
			suanpan_info("create_new_step() requires a node.\n");
			return SUANPAN_SUCCESS;
		}

		unsigned dof;
		if(!get_input(command, dof)) {
			suanpan_info("create_new_step() requires a dof.\n");
			return SUANPAN_SUCCESS;
		}

		double magnitude;
		if(!get_input(command, magnitude)) {
			suanpan_info("create_new_step() requires a magnitude.\n");
			return SUANPAN_SUCCESS;
		}

		if(domain->insert(make_shared<ArcLength>(tag, node, dof, magnitude))) domain->set_current_step_tag(tag);
		else suanpan_error("create_new_step() cannot create the new step.\n");
	} else suanpan_info("create_new_step() cannot identify step type.\n");

	return SUANPAN_SUCCESS;
}

int test_material1d(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned material_tag;
	if(!get_input(command, material_tag)) {
		suanpan_error("test_material1d() needs a valid material tag.\n");
		return SUANPAN_SUCCESS;
	}

	double incre;
	if(!get_input(command, incre)) {
		suanpan_error("test_material1d() needs a valid step size.\n");
		return SUANPAN_SUCCESS;
	}

	vector<unsigned> load_step;
	int step;
	while(get_input(command, step)) load_step.push_back(step < 0 ? -step : step);

	if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

	auto& material_proto = domain->get_material(material_tag);

	const auto result = material_tester(material_proto->get_copy(), load_step, {incre});

#ifdef SUANPAN_HDF5
	if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

	if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

	std::ofstream gnuplot("RESULT.plt");

	if(gnuplot.is_open()) {
		gnuplot << "reset\n";
		gnuplot << "set term tikz size 14cm,10cm\n";
		gnuplot << "set output \"RESULT.tex\"\n";
		gnuplot << "unset key\n";
		gnuplot << "set xrange [*:*]\n";
		gnuplot << "set yrange [*:*]\n";
		gnuplot << "set xlabel \"input\"\n";
		gnuplot << "set ylabel \"output\"\n";
		gnuplot << "set grid\n";
		gnuplot << "plot \"RESULT.txt\" u 1:2 w l lw 2\n";
		gnuplot << "set output\n";

		gnuplot.close();
	}

	return SUANPAN_SUCCESS;
}

int test_material2d(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned material_tag;
	if(!get_input(command, material_tag)) {
		suanpan_error("test_material2d() needs a valid material tag.\n");
		return SUANPAN_SUCCESS;
	}

	vec incre(3);
	for(auto I = 0; I < 3; ++I) {
		if(!get_input(command, incre(I))) {
			suanpan_error("test_material2d() needs a valid step size.\n");
			return SUANPAN_SUCCESS;
		}
	}

	vector<unsigned> load_step;
	int step;
	while(get_input(command, step)) load_step.push_back(step < 0 ? -step : step);

	if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

	auto& material_proto = domain->get_material(material_tag);

	const auto result = material_tester(material_proto->get_copy(), load_step, incre);

#ifdef SUANPAN_HDF5
	if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

	if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

	return SUANPAN_SUCCESS;
}

int test_material3d(const shared_ptr<DomainBase>& domain, istringstream& command) {
	unsigned material_tag;
	if(!get_input(command, material_tag)) {
		suanpan_error("test_material3d() needs a valid material tag.\n");
		return SUANPAN_SUCCESS;
	}

	vec incre(6);
	for(auto I = 0; I < 6; ++I)
		if(!get_input(command, incre(I))) {
			suanpan_error("test_material3d() needs a valid step size.\n");
			return SUANPAN_SUCCESS;
		}

	vector<unsigned> load_step;
	int step;
	while(get_input(command, step)) load_step.push_back(step < 0 ? -step : step);

	if(!domain->find_material(material_tag)) return SUANPAN_SUCCESS;

	auto& material_proto = domain->get_material(material_tag);

	const auto result = material_tester(material_proto->get_copy(), load_step, incre);

#ifdef SUANPAN_HDF5
	if(!result.save("RESULT.h5", hdf5_binary_trans)) suanpan_error("fail to save file.\n");
#endif

	if(!result.save("RESULT.txt", raw_ascii)) suanpan_error("fail to save file.\n");

	return SUANPAN_SUCCESS;
}

int set_property(const shared_ptr<DomainBase>& domain, istringstream& command) {
	if(domain->get_current_step_tag() == 0) return SUANPAN_SUCCESS;

	const auto& tmp_step = domain->get_current_step();

	string property_id;
	if(!get_input(command, property_id)) {
		suanpan_info("set_property() need a property type.\n");
		return SUANPAN_SUCCESS;
	}

	if(is_equal(property_id, "fixed_step_size")) {
		string value;
		get_input(command, value) ? tmp_step->set_fixed_step_size(is_true(value)) : suanpan_info("set_property() need a valid value.\n");
	} else if(is_equal(property_id, "symm_mat")) {
		string value;
		get_input(command, value) ? tmp_step->set_symm(is_true(value)) : suanpan_info("set_property() need a valid value.\n");
	} else if(is_equal(property_id, "band_mat")) {
		string value;
		get_input(command, value) ? tmp_step->set_band(is_true(value)) : suanpan_info("set_property() need a valid value.\n");
	} else if(is_equal(property_id, "sparse_mat")) {
		string value;
		get_input(command, value) ? tmp_step->set_sparse(is_true(value)) : suanpan_info("set_property() need a valid value.\n");
	} else if(is_equal(property_id, "ini_step_size")) {
		double step_time;
		get_input(command, step_time) ? tmp_step->set_ini_step_size(step_time) : suanpan_info("set_property() need a valid value.\n");
	} else if(is_equal(property_id, "min_step_size")) {
		double step_time;
		get_input(command, step_time) ? tmp_step->set_min_step_size(step_time) : suanpan_info("set_property() need a valid value.\n");
	} else if(is_equal(property_id, "max_step_size")) {
		double step_time;
		get_input(command, step_time) ? tmp_step->set_max_step_size(step_time) : suanpan_info("set_property() need a valid value.\n");
	} else if(is_equal(property_id, "max_iteration")) {
		unsigned max_number;
		get_input(command, max_number) ? tmp_step->set_max_substep(max_number) : suanpan_info("set_property() need a valid value.\n");
	} else if(is_equal(property_id, "eigen_number")) {
		unsigned eigen_number;
		if(get_input(command, eigen_number)) {
			const auto eigen_step = std::dynamic_pointer_cast<Frequency>(tmp_step);
			if(eigen_step == nullptr) suanpan_info("set_property() cannot set eigen number for noneigen step.\n");
			else eigen_step->set_eigen_number(eigen_number);
		} else suanpan_info("set_property() need a valid eigen number.\n");
	}

	return SUANPAN_SUCCESS;
}

int print_info(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string object_type;
	if(!get_input(command, object_type)) {
		suanpan_info("print_info() needs object type.\n");
		return SUANPAN_SUCCESS;
	}

	unsigned tag;
	if(is_equal(object_type, "node"))
		while(get_input(command, tag) && domain->find_node(tag)) {
			get_node(domain, tag)->print();
			suanpan_info("\n");
		}
	else if(is_equal(object_type, "element"))
		while(get_input(command, tag) && domain->find_element(tag)) {
			get_element(domain, tag)->print();
			suanpan_info("\n");
		}
	else if(is_equal(object_type, "material"))
		while(get_input(command, tag) && domain->find_material(tag)) {
			get_material(domain, tag)->print();
			suanpan_info("\n");
		}
	else if(is_equal(object_type, "constraint"))
		while(get_input(command, tag) && domain->find_constraint(tag)) {
			get_constraint(domain, tag)->print();
			suanpan_info("\n");
		}
	else if(is_equal(object_type, "recorder"))
		while(get_input(command, tag) && domain->find_recorder(tag)) {
			get_recorder(domain, tag)->print();
			suanpan_info("\n");
		}
	else if(is_equal(object_type, "solver"))
		while(get_input(command, tag) && domain->find_solver(tag)) {
			get_solver(domain, tag)->print();
			suanpan_info("\n");
		}
	else if(is_equal(object_type, "integrator"))
		while(get_input(command, tag) && domain->find_integrator(tag)) {
			get_integrator(domain, tag)->print();
			suanpan_info("\n");
		}
	else if(is_equal(object_type, "eigenvalue")) {
		suanpan_info("Eigenvalues:\n");
		domain->get_factory()->get_eigenvalue().print();
		suanpan_info("\n");
	}

	return SUANPAN_SUCCESS;
}
