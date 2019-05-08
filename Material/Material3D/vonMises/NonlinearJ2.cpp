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

#include "NonlinearJ2.h"
#include <Toolbox/tensorToolbox.h>

const double NonlinearJ2::two_third = 2. / 3.;
const double NonlinearJ2::root_two_third = sqrt(two_third);
const mat NonlinearJ2::unit_dev_tensor = tensor::unit_deviatoric_tensor4();
const unsigned NonlinearJ2::max_iteration = 20;

NonlinearJ2::NonlinearJ2(const unsigned T, const double E, const double V, const double R)
	: DataNonlinearJ2{E, V, E / (2. + 2. * V), E / (1. + V), pow(E / (1. + V), 2.)}
	, Material3D(T, R) {}

void NonlinearJ2::initialize(const shared_ptr<DomainBase>&) {
	trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

	trial_history = current_history.zeros(7);
}

double NonlinearJ2::get_parameter(const ParameterType P) const {
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

int NonlinearJ2::update_trial_status(const vec& t_strain) {
	trial_stress = current_stress + (trial_stiffness = initial_stiffness) * (incre_strain = (trial_strain = t_strain) - current_strain);

	trial_history = current_history;
	auto& plastic_strain = trial_history(0);
	vec back_stress(&trial_history(1), 6, false, true);

	const vec rel_stress = tensor::dev(trial_stress) - back_stress;
	const auto norm_rel_stress = tensor::stress::norm(rel_stress);

	auto yield_func = norm_rel_stress - root_two_third * std::max(0., compute_k(plastic_strain));

	if(yield_func < 0.) return SUANPAN_SUCCESS;

	const auto current_h = compute_h(plastic_strain);
	auto gamma = 0., incre_h = 0., denom = 0.;
	unsigned counter = 0;
	while(++counter < max_iteration) {
		denom = double_shear + two_third * (compute_dk(plastic_strain) + compute_dh(plastic_strain));
		const auto incre_gamma = yield_func / denom;
		const auto abs_error = fabs(incre_gamma);
		suanpan_extra_debug("NonlinearJ2 local iterative loop error: %.5E.\n", abs_error);
		if(abs_error <= tolerance) break;
		incre_h = compute_h(plastic_strain = current_history(0) + root_two_third * (gamma += incre_gamma)) - current_h;
		yield_func = norm_rel_stress - double_shear * gamma - root_two_third * (std::max(0., compute_k(plastic_strain)) + incre_h);
	}

	if(max_iteration == counter) {
		suanpan_error("cannot converge in %u iterations.\n", max_iteration);
		return SUANPAN_FAIL;
	}

	back_stress += root_two_third * incre_h / norm_rel_stress * rel_stress;

	auto t_factor = double_shear * gamma / norm_rel_stress;
	trial_stress -= t_factor * rel_stress;

	t_factor *= double_shear;
	trial_stiffness += (t_factor - square_double_shear / denom) / norm_rel_stress / norm_rel_stress * rel_stress * rel_stress.t() - t_factor * unit_dev_tensor;

	return SUANPAN_SUCCESS;
}

int NonlinearJ2::clear_status() {
	current_strain.zeros();
	current_stress.zeros();
	current_history.zeros();
	current_stiffness = initial_stiffness;
	return reset_status();
}

int NonlinearJ2::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_history = trial_history;
	current_stiffness = trial_stiffness;
	return SUANPAN_SUCCESS;
}

int NonlinearJ2::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_history = current_history;
	trial_stiffness = current_stiffness;
	return SUANPAN_SUCCESS;
}

void NonlinearJ2::print() { suanpan_info("A 3D nonlinear hardening model using von-Mises yielding criterion.\n"); }
