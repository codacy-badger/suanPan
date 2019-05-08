﻿////////////////////////////////////////////////////////////////////////////////
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

#include "Sequential.h"
#include <Domain/DomainBase.h>

const unsigned Sequential::max_iteration = 20;

Sequential::Sequential(const unsigned T, uvec&& MT)
	: Material1D(T, 0.)
	, mat_size(MT.n_elem - 1)
	, mat_tag(std::forward<uvec>(MT)) {}

Sequential::Sequential(const Sequential& old_obj)
	: Material1D(old_obj)
	, mat_size(old_obj.mat_size)
	, mat_tag(old_obj.mat_tag)
	, jacobian(old_obj.jacobian) { for(const auto& I : old_obj.mat_pool) if(I != nullptr) mat_pool.emplace_back(I->get_copy()); }

void Sequential::initialize(const shared_ptr<DomainBase>& D) {
	if(D == nullptr) return;

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

	jacobian.zeros(mat_tag.n_elem, mat_tag.n_elem);
	jacobian.row(0).fill(1.);

	initial_stiffness.zeros(1);
	for(const auto& I : mat_pool) {
		I->Material::initialize(D);
		I->initialize(D);
		initial_stiffness += 1. / I->get_initial_stiffness();
	}
	trial_stiffness = current_stiffness = initial_stiffness = 1. / initial_stiffness;
}

unique_ptr<Material> Sequential::get_copy() { return make_unique<Sequential>(*this); }

int Sequential::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

	unsigned counter = 0;
	while(++counter < max_iteration) {
		vec residual(mat_tag.n_elem, fill::zeros);
		residual(0) = trial_strain(0) - mat_pool.front()->get_trial_strain().at(0) - mat_pool.back()->get_trial_strain().at(0);
		residual(1) = mat_pool.front()->get_trial_stress().at(0);
		residual(mat_size) -= mat_pool.back()->get_trial_stress().at(0);

		jacobian(1, 0) = -mat_pool.front()->get_trial_stiffness().at(0);
		jacobian(mat_size, mat_size) = mat_pool.back()->get_trial_stiffness().at(0);

		for(uword I = 1; I < mat_size; ++I) {
			residual(0) -= mat_pool[I]->get_trial_strain().at(0);
			residual(I) -= mat_pool[I]->get_trial_stress().at(0);
			residual(I + 1) += mat_pool[I]->get_trial_stress().at(0);
			jacobian(I + 1, I) = -(jacobian(I, I) = mat_pool[I]->get_trial_stiffness().at(0));
		}

		const vec i_strain = solve(jacobian, residual);

		for(size_t I = 0; I < mat_pool.size(); ++I) mat_pool[I]->update_trial_status(mat_pool[I]->get_trial_strain() + i_strain[I]);

		const auto error = norm(i_strain);
		suanpan_debug("Sequential material local iteration error: %.4E.\n", error);
		if(error <= tolerance) break;
	}

	if(max_iteration == counter) {
		suanpan_error("Sequential material cannot converge within %u iterations.\n", max_iteration);
		return SUANPAN_FAIL;
	}

	trial_stress = mat_pool.front()->get_trial_stress();

	auto t_stiff = 0.;
	for(const auto& I : mat_pool) t_stiff += 1. / I->get_trial_stiffness().at(0);

	trial_stiffness = 1. / t_stiff;

	return SUANPAN_SUCCESS;
}

int Sequential::clear_status() {
	auto code = 0;
	for(const auto& I : mat_pool) code += I->clear_status();
	return code;
}

int Sequential::commit_status() {
	auto code = 0;
	for(const auto& I : mat_pool) code += I->commit_status();
	return code;
}

int Sequential::reset_status() {
	auto code = 0;
	for(const auto& I : mat_pool) code += I->reset_status();
	return code;
}

vector<vec> Sequential::record(const OutputType P) {
	vector<vec> data;

	for(const auto& I : mat_pool) for(const auto& J : I->record(P)) data.emplace_back(J);

	return data;
}

void Sequential::print() {
	suanpan_info("A wrapper of several material models.\n");
	for(const auto& I : mat_pool) I->print();
}
