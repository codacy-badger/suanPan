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

#include "NonlinearPerzyna.h"
#include <Toolbox/tensorToolbox.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

const unsigned NonlinearPerzyna::max_iteration = 20;
const double NonlinearPerzyna::root_three_two = sqrt(1.5);
const mat NonlinearPerzyna::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

NonlinearPerzyna::NonlinearPerzyna(const unsigned T, const double E, const double V, const double MU, const double EPS, const double R)
	: DataNonlinearPerzyna{E, V, MU, EPS, E / (2. + 2. * V), E / (1. + V), 3. * E / (2. + 2. * V), MU * EPS}
	, Material3D(T, R) {}

void NonlinearPerzyna::initialize(const shared_ptr<DomainBase>& D) {
	if(D == nullptr) return;

	incre_time = &D->get_factory()->get_incre_time();

	trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

	trial_history = current_history.zeros(1);
}

double NonlinearPerzyna::get_parameter(const ParameterType P) const {
	switch(P) {
	case ParameterType::E:
	case ParameterType::ELASTICMODULUS:
	case ParameterType::YOUNGSMODULUS:
		return elastic_modulus;
	case ParameterType::POISSONSRATIO:
		return poissons_ratio;
	case ParameterType::SHEARMODULUS:
	case ParameterType::G:
		return shear_modulus;
	case ParameterType::DENSITY:
		return density;
	default:
		return 0.;
	}
}

int NonlinearPerzyna::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	if(norm(incre_strain) <= tolerance) return SUANPAN_SUCCESS;

	trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

	trial_history = current_history;
	auto& plastic_strain = trial_history(0);

	const auto dev_stress = tensor::dev(trial_stress);
	const auto eqv_stress = root_three_two * tensor::stress::norm(dev_stress);

	auto residual = eqv_stress - std::max(0., compute_k(plastic_strain));

	if(residual < 0.) return SUANPAN_SUCCESS;

	auto gamma = 0., dk = 0., pow_term = 1., denom = *incre_time;
	unsigned counter = 0;
	while(++counter < max_iteration) {
		const auto gradient = (triple_shear + factor_a * (eqv_stress - triple_shear * gamma) / denom) * pow_term - (dk = compute_dk(plastic_strain));
		const auto incre_gamma = residual / gradient;
		suanpan_extra_debug("NonlinearPerzyna local iterative loop error: %.5E.\n", incre_gamma);
		if(fabs(incre_gamma) <= tolerance) break;
		plastic_strain = current_history(0) + (gamma += incre_gamma);
		denom = *incre_time + mu * gamma;
		pow_term = pow(*incre_time / denom, epsilon);
		residual = (eqv_stress - triple_shear * gamma) * pow_term - std::max(0., compute_k(plastic_strain));
	}

	if(max_iteration == counter) {
		suanpan_error("cannot converge in %u iterations.\n", max_iteration);
		return SUANPAN_FAIL;
	}

	trial_stress -= gamma * triple_shear / eqv_stress * dev_stress;

	trial_stiffness += triple_shear * triple_shear * (gamma / eqv_stress - 1. / (triple_shear + dk / pow_term + factor_a / denom * (eqv_stress - triple_shear * gamma))) / eqv_stress / eqv_stress * dev_stress * dev_stress.t() - triple_shear * double_shear * gamma / eqv_stress * unit_dev_tensor;

	return SUANPAN_SUCCESS;
}

int NonlinearPerzyna::clear_status() {
	current_strain.zeros();
	current_stress.zeros();
	current_history.zeros();
	current_stiffness = initial_stiffness;
	return reset_status();
}

int NonlinearPerzyna::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_history = trial_history;
	current_stiffness = trial_stiffness;
	return SUANPAN_SUCCESS;
}

int NonlinearPerzyna::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_history = current_history;
	trial_stiffness = current_stiffness;
	return SUANPAN_SUCCESS;
}

void NonlinearPerzyna::print() { suanpan_info("A 3D bilinear hardening viscoplasticity model using Perzyna rule.\n"); }
