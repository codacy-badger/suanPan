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

#include "Laminated.h"
#include <Domain/DomainBase.h>

Laminated::Laminated(const unsigned T, uvec&& MT)
	: Material2D(T)
	, mat_tag(std::forward<uvec>(MT)) {}

Laminated::Laminated(const Laminated& old_obj)
	: Material2D(old_obj)
	, mat_tag(old_obj.mat_tag) {
	mat_pool.clear(), mat_pool.reserve(old_obj.mat_pool.size());
	for(const auto& I : old_obj.mat_pool) if(I != nullptr) mat_pool.emplace_back(I->get_copy());
}

void Laminated::initialize(const shared_ptr<DomainBase>& D) {
	if(D == nullptr) return;

	initial_stiffness.zeros(3, 3);
	mat_pool.clear(), mat_pool.reserve(mat_tag.n_elem);
	for(const auto I : mat_tag) {
		if(!D->find_material(unsigned(I)) || D->get_material(unsigned(I))->get_material_type() != MaterialType::D2) {
			D->disable_material(get_tag());
			return;
		}
		mat_pool.emplace_back(D->get_material(unsigned(I))->get_copy());
		mat_pool.back()->Material::initialize(D);
		mat_pool.back()->initialize(D);
		access::rw(density) += mat_pool.back()->get_parameter(ParameterType::DENSITY);
		initial_stiffness += mat_pool.back()->get_initial_stiffness();
	}

	trial_stiffness = current_stiffness = initial_stiffness;
}

unique_ptr<Material> Laminated::get_copy() { return make_unique<Laminated>(*this); }

int Laminated::update_trial_status(const vec& t_strain) {
	trial_strain = t_strain;

	trial_stress.zeros(3);
	trial_stiffness.zeros(3, 3);
	for(const auto& I : mat_pool) {
		if(I->update_trial_status(trial_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		trial_stress += I->get_trial_stress();
		trial_stiffness += I->get_trial_stiffness();
	}

	return SUANPAN_SUCCESS;
}

int Laminated::clear_status() {
	current_strain.zeros();
	trial_strain.zeros();
	current_stress.zeros();
	trial_stress.zeros();
	trial_stiffness = current_stiffness = initial_stiffness;
	auto code = 0;
	for(const auto& I : mat_pool) code += I->clear_status();
	return code;
}

int Laminated::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_stiffness = trial_stiffness;
	auto code = 0;
	for(const auto& I : mat_pool) code += I->commit_status();
	return code;
}

int Laminated::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_stiffness = current_stiffness;
	auto code = 0;
	for(const auto& I : mat_pool) code += I->reset_status();
	return code;
}

vector<vec> Laminated::record(const OutputType P) {
	vector<vec> data;

	for(const auto& I : mat_pool) for(const auto& J : I->record(P)) data.emplace_back(J);

	return data;
}

void Laminated::print() {
	suanpan_info("A multilayer wrapper for planar problems.\n");
	unsigned t_tag = 0;
	for(const auto& I : mat_pool) {
		suanpan_info("component %u: ", ++t_tag);
		I->print();
	}
}
