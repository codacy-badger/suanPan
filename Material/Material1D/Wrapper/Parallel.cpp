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

#include "Parallel.h"
#include <Domain/DomainBase.h>

Parallel::Parallel(const unsigned T, uvec&& MT)
	: Material1D(T, 0.)
	, mat_tag(std::forward<uvec>(MT)) {}

Parallel::Parallel(const Parallel& old_obj)
	: Material1D(old_obj)
	, mat_tag(old_obj.mat_tag) { for(const auto& I : old_obj.mat_pool) if(I != nullptr) mat_pool.emplace_back(I->get_copy()); }

void Parallel::initialize(const shared_ptr<DomainBase>& D) {
	if(D == nullptr) {
		D->disable_material(get_tag());
		return;
	}

	mat_pool.clear(), mat_pool.reserve(mat_tag.n_elem);
	for(unsigned I = 0; I < mat_tag.n_elem; ++I) {
		const auto t_tag = unsigned(mat_tag(I));
		if(!D->find_material(t_tag) || D->get_material(t_tag)->get_material_type() != MaterialType::D1) {
			D->disable_material(get_tag());
			return;
		}
		mat_pool.emplace_back(get_material(D, t_tag)->get_copy());
		access::rw(density) += mat_pool.back()->get_parameter(ParameterType::DENSITY);
	}

	initial_stiffness.zeros(1);
	for(const auto& I : mat_pool) {
		I->Material::initialize(D);
		I->initialize(D);
		initial_stiffness += I->get_initial_stiffness();
	}
	trial_stiffness = current_stiffness = initial_stiffness;
}

unique_ptr<Material> Parallel::get_copy() { return make_unique<Parallel>(*this); }

int Parallel::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

	trial_stress.zeros();
	trial_stiffness.zeros();
	for(const auto& I : mat_pool) {
		if(I->update_trial_status(trial_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		trial_stress += I->get_trial_stress();
		trial_stiffness += I->get_trial_stiffness();
	}

	return SUANPAN_SUCCESS;
}

int Parallel::clear_status() {
	auto code = 0;
	for(const auto& I : mat_pool) code += I->clear_status();
	return code;
}

int Parallel::commit_status() {
	auto code = 0;
	for(const auto& I : mat_pool) code += I->commit_status();
	return code;
}

int Parallel::reset_status() {
	auto code = 0;
	for(const auto& I : mat_pool) code += I->reset_status();
	return code;
}

vector<vec> Parallel::record(const OutputType P) {
	vector<vec> data;

	for(const auto& I : mat_pool) for(const auto& J : I->record(P)) data.emplace_back(J);

	return data;
}

void Parallel::print() {
	suanpan_info("A Parallel container that holds following models: ");
	mat_tag.t().print();
	for(const auto& I : mat_pool) I->print();
}
