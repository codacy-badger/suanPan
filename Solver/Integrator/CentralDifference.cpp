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

#include "CentralDifference.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Toolbox/arpack_wrapper.h>
#include <future>

CentralDifference::CentralDifference(const unsigned T)
	: Integrator(T) {}

int CentralDifference::initialize() {
	const auto code = Integrator::initialize();

	if(code == 0) {
		const auto& W = get_domain().lock()->get_factory();

		vec eig_val;
		mat eig_vec;

		eig_solve(eig_val, eig_vec, W->get_stiffness(), W->get_mass(), std::min(unsigned(10), W->get_size() / 2));

		max_dt = datum::pi / sqrt(eig_val.min());
	}

	return code;
}

void CentralDifference::update_parameter() {
	const auto& W = get_domain().lock()->get_factory();

	if(DT != W->get_incre_time() || W->get_pre_displacement().is_empty()) {
		DT = W->get_incre_time();
		if(DT > max_dt) suanpan_error("update_status() requires a smaller time increment.\n");

		C0 = 1. / DT / DT;
		C1 = .5 / DT;
		C2 = 2. * C0;
		C3 = 1. / C2;

		W->set_pre_displacement(W->get_current_displacement() - DT * W->get_current_velocity() + C3 * W->get_current_acceleration());
	}
}

void CentralDifference::assemble_resistance() {
	update_parameter();

	const auto& D = get_domain().lock();
	const auto& W = D->get_factory();

	// done in update_trail_status()
	// D->assemble_resistance();

	D->assemble_mass();
	D->assemble_damping();

	auto t_vector_a(std::async([&]() { return vec{get_mass(W) * (C0 * (W->get_trial_displacement() + W->get_pre_displacement()) - C2 * W->get_current_displacement())}; }));
	auto t_vector_b(std::async([&]() { return vec{get_damping(W) * (C1 * (W->get_trial_displacement() - W->get_pre_displacement()))}; }));

	get_sushi(W) += t_vector_a.get() + t_vector_b.get();
}

void CentralDifference::assemble_matrix() {
	update_parameter();

	const auto& W = get_domain().lock()->get_factory();

	auto& t_stiff = get_stiffness(W);

	1. == C0 ? t_stiff += get_mass(W) : -1 == C0 ? t_stiff -= get_mass(W) : t_stiff += C0 * get_mass(W);
	1. == C1 ? t_stiff += get_damping(W) : -1 == C1 ? t_stiff -= get_damping(W) : t_stiff += C1 * get_damping(W);
}

void CentralDifference::commit_status() const {
	const auto& D = get_domain().lock();
	const auto& W = D->get_factory();

	W->update_current_velocity(C1 * (W->get_trial_displacement() - W->get_pre_displacement()));

	W->update_current_acceleration(C0 * (W->get_pre_displacement() - 2. * W->get_current_displacement() + W->get_trial_displacement()));

	W->commit_pre_displacement();

	D->commit_status();
}