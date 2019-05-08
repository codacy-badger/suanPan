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

#include "Maxwell.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Recorder/OutputType.h>

const unsigned Maxwell::max_iteration = 20;

Maxwell::Maxwell(const unsigned T, const unsigned DT, const unsigned ST, const bool UM, const unsigned PC)
	: Material1D(T, 0.)
	, damper_tag(DT)
	, spring_tag(ST)
	, proceed(PC)
	, use_matrix(UM) { access::rw(tolerance) = 1E-11; }

Maxwell::Maxwell(const Maxwell& old_obj)
	: Material1D(old_obj)
	, incre_time(old_obj.incre_time)
	, damper_tag(old_obj.damper_tag)
	, spring_tag(old_obj.spring_tag)
	, proceed(old_obj.proceed)
	, use_matrix(old_obj.use_matrix)
	, damper(old_obj.damper == nullptr ? nullptr : old_obj.damper->get_copy())
	, spring(old_obj.spring == nullptr ? nullptr : old_obj.spring->get_copy()) {}

void Maxwell::initialize(const shared_ptr<DomainBase>& D) {
	if(D == nullptr) return;

	if(!D->find_material(damper_tag) || !D->find_material(spring_tag)) {
		D->disable_material(get_tag());
		return;
	}

	incre_time = &D->get_factory()->get_incre_time();

	damper = get_material(D, damper_tag)->get_copy();
	spring = get_material(D, spring_tag)->get_copy();

	damper->Material::initialize(D);
	damper->initialize(D);
	spring->Material::initialize(D);
	spring->initialize(D);

	trial_strain_rate = current_strain_rate = incre_strain_rate.zeros(1);

	const auto& K1 = spring->get_initial_stiffness().at(0);
	const auto& K2 = damper->get_initial_stiffness().at(0);
	const auto& K3 = damper->get_initial_damping().at(0);

	trial_damping = current_damping = initial_damping = K3;
	trial_stiffness = current_stiffness = initial_stiffness = K1 * K2 / (K1 + K2);
}

unique_ptr<Material> Maxwell::get_copy() { return make_unique<Maxwell>(*this); }

int Maxwell::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
	incre_strain = (trial_strain = t_strain) - current_strain;
	incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

	if(fabs(incre_strain(0) + fabs(incre_strain_rate(0))) <= datum::eps) return SUANPAN_SUCCESS;

	const auto& K1 = spring->get_trial_stiffness().at(0);
	const auto& K2 = damper->get_trial_stiffness().at(0);
	const auto& K3 = damper->get_trial_damping().at(0);
	const auto& F1 = spring->get_trial_stress().at(0);
	const auto& F2 = damper->get_trial_stress().at(0);

	// \beta\Delta{}t
	const auto factor_a = (incre_strain(0) - *incre_time * current_strain_rate(0)) / incre_strain_rate(0);

	const auto target = incre_strain(0) - *incre_time * damper->get_current_strain_rate().at(0);

	vec solution(3, fill::zeros);

	double error;
	auto ref_error = 1.;

	counter = 0;

	if(use_matrix) {
		mat jacobian(3, 3, fill::zeros);
		jacobian(0, 0) = jacobian(0, 1) = jacobian(1, 0) = 1.;
		jacobian(1, 2) = factor_a;

		while(++counter < max_iteration) {
			const vec residual{incre_strain(0) - solution(0) - solution(1), target - solution(0) - factor_a * solution(2), F1 - F2};
			jacobian(2, 0) = -K1;
			jacobian(2, 1) = K2;
			jacobian(2, 2) = K3;
			const vec incre = solve(jacobian, residual, solve_opts::equilibrate);
			if(1 == counter) ref_error = std::max(1., norm(residual));
			suanpan_extra_debug("Maxwell local iteration error: %.4E.\n", error = norm(residual) / ref_error);
			if(norm(incre) <= tolerance && error <= tolerance) break;
			solution += incre;
			spring->update_incre_status(solution(0));
			damper->update_incre_status(solution(1), solution(2));
		}
	} else
		while(++counter < max_iteration) {
			const auto residual_a = incre_strain(0) - solution(0) - solution(1);
			const auto residual_b = target - solution(0) - factor_a * solution(2);
			const auto residual_c = F1 - F2;
			const auto residual = residual_a * K2 + residual_b / factor_a * K3 - residual_c;
			const auto jacobian = K1 + K2 + K3 / factor_a;
			const auto incre = residual / jacobian;
			if(1 == counter) ref_error = std::max(1., fabs(residual));
			suanpan_extra_debug("Maxwell local iteration error: %.4E.\n", error = fabs(residual) / ref_error);
			if(fabs(incre) <= tolerance && error <= tolerance) break;
			solution(0) += incre;
			solution(1) += residual_a - incre;
			solution(2) += (residual_b - incre) / factor_a;
			spring->update_incre_status(solution(0));
			damper->update_incre_status(solution(1), solution(2));
		}

	if(max_iteration != counter) delay_counter = 0;
	else if(1 == proceed) reset_status();
	else if(0 == proceed || ++delay_counter == proceed) {
		suanpan_warning("Maxwell local iteration cannot converge within %u interations.\n", max_iteration);
		return SUANPAN_FAIL;
	}

	trial_stress = .5 * (spring->get_trial_stress() + damper->get_trial_stress());

	trial_damping = K3;
	trial_stiffness = K1 * K2 / (K1 + K2);

	return SUANPAN_SUCCESS;
}

int Maxwell::clear_status() {
	trial_strain = current_strain.zeros();
	trial_stress = current_stress.zeros();
	trial_strain_rate = current_strain_rate.zeros();
	trial_damping = current_damping = initial_damping;
	trial_stiffness = current_stiffness = initial_stiffness;
	return spring->clear_status() + damper->clear_status();
}

int Maxwell::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_strain_rate = trial_strain_rate;
	current_damping = trial_damping;
	current_stiffness = trial_stiffness;
	return spring->commit_status() + damper->commit_status();
}

int Maxwell::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_strain_rate = current_strain_rate;
	trial_damping = current_damping;
	trial_stiffness = current_stiffness;
	return spring->reset_status() + damper->reset_status();
}

vector<vec> Maxwell::record(const OutputType P) {
	vector<vec> data;

	if(OutputType::SD == P) data.emplace_back(current_stress);
	else if(OutputType::SS == P) data.emplace_back(current_stress);
	else if(OutputType::S == P) data.emplace_back(current_stress);
	else if(OutputType::ED == P) data.emplace_back(damper->get_current_strain());
	else if(OutputType::VD == P) data.emplace_back(damper->get_current_strain_rate());
	else if(OutputType::ES == P) data.emplace_back(spring->get_current_strain());
	else if(OutputType::VS == P) data.emplace_back(current_strain_rate - damper->get_current_strain_rate());
	else if(OutputType::E == P) data.emplace_back(current_strain);
	else if(OutputType::V == P) data.emplace_back(current_strain_rate);
	else if(OutputType::LITR == P) data.emplace_back(vec{double(counter)});

	return data;
}

void Maxwell::print() { suanpan_info("A Maxwell material model.\n"); }
