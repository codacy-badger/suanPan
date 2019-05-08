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

#include "MaterialTester.h"

mat material_tester(const shared_ptr<Material>& obj, const vector<unsigned>& idx, const vec& incre) {
	if(!obj->is_initialized()) {
		obj->Material::initialize();
		obj->initialize();
		obj->set_initialized(true);
	}

	if(obj->get_material_type() == MaterialType::D1 && incre.n_elem != 1) {
		suanpan_error("the tester cannot be appiled to the given material model.\n");
		return {};
	}
	if(obj->get_material_type() == MaterialType::D2 && incre.n_elem != 3) {
		suanpan_error("the tester cannot be appiled to the given material model.\n");
		return {};
	}
	if(obj->get_material_type() == MaterialType::D3 && incre.n_elem != 6) {
		suanpan_error("the tester cannot be appiled to the given material model.\n");
		return {};
	}

	unsigned total_size = 1;
	for(const auto& I : idx) total_size += I;

	mat response(total_size, 2 * incre.n_elem, fill::zeros);

	const span span_a(0, incre.n_elem - 1);
	const span span_b(incre.n_elem, 2 * incre.n_elem - 1);

	auto current_pos = 0;

	response(current_pos, span_a) = obj->get_current_strain().t();
	response(current_pos++, span_b) = obj->get_current_stress().t();

	auto incre_strain = incre;
	for(const auto& I : idx) {
		for(unsigned J = 0; J < I; ++J) {
			if(obj->update_incre_status(incre_strain) != SUANPAN_SUCCESS) break;
			obj->commit_status();
			response(current_pos, span_a) = obj->get_current_strain().t();
			response(current_pos++, span_b) = obj->get_current_stress().t();
		}
		incre_strain = -incre_strain;
	}

	obj->print();
	obj->reset_status();
	obj->clear_status();

	return response;
}
