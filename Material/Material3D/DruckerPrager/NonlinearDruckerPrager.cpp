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

#include "NonlinearDruckerPrager.h"
#include <Toolbox/tensorToolbox.h>

const unsigned NonlinearDruckerPrager::max_iteration = 20;
const mat NonlinearDruckerPrager::unit_dev_tensor = tensor::unit_deviatoric_tensor4();
const mat NonlinearDruckerPrager::unit_x_unit = tensor::unit_tensor2 * tensor::unit_tensor2.t();

NonlinearDruckerPrager::NonlinearDruckerPrager(const unsigned T, const double E, const double V, const double ETAY, const double ETAF, const double XI, const double R)
	: DataNonlinearDruckerPrager{E, V, E / (2. + 2. * V), E / (3. - 6. * V), ETAY, ETAF, XI, E / (1. + V), 0., 0., 0., 0.}
	, Material3D(T, R) {
	access::rw(factor_a) = shear + bulk * eta_flow * eta_yield;
	access::rw(factor_b) = xi * xi / eta_flow / eta_yield;
	access::rw(factor_c) = sqrt(2.) * shear * bulk;
	access::rw(factor_d) = bulk * bulk * eta_flow * eta_yield;
}

void NonlinearDruckerPrager::initialize(const shared_ptr<DomainBase>&) {
	trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

	trial_history = current_history.zeros(1);
}

double NonlinearDruckerPrager::get_parameter(const ParameterType P) const {
	switch(P) {
	case ParameterType::E:
	case ParameterType::ELASTICMODULUS:
	case ParameterType::YOUNGSMODULUS:
		return elastic_modulus;
	case ParameterType::POISSONSRATIO:
		return poissons_ratio;
	case ParameterType::SHEARMODULUS:
	case ParameterType::G:
		return shear;
	case ParameterType::DENSITY:
		return density;
	default:
		return 0.;
	}
}

int NonlinearDruckerPrager::update_trial_status(const vec& t_strain) {
	trial_stress = current_stress + (trial_stiffness = initial_stiffness) * (incre_strain = (trial_strain = t_strain) - current_strain);

	trial_history = current_history;
	auto& plastic_strain = trial_history(0);

	const auto dev_stress = tensor::dev(trial_stress);
	const auto hydro_stress = tensor::mean(trial_stress);
	const auto sqrt_j2 = sqrt(tensor::stress::invariant2(dev_stress));

	const auto yield_const = sqrt_j2 + eta_yield * hydro_stress;

	if(yield_const - xi * compute_c(plastic_strain) <= 0.) return SUANPAN_SUCCESS;

	auto gamma = 0., denominator = 0.;

	unsigned counter = 0;
	while(++counter < max_iteration) {
		const auto incre_gamma = (yield_const - factor_a * gamma - xi * compute_c(plastic_strain)) / (denominator = factor_a + xi * xi * compute_dc(plastic_strain));
		suanpan_extra_debug("NonlinearDruckerPrager local iterative loop error: %.5E.\n", incre_gamma);
		if(fabs(incre_gamma) <= tolerance) break;
		plastic_strain = current_history(0) + xi * (gamma += incre_gamma);
	}

	if(max_iteration == counter) {
		suanpan_error("cannot converge in %u iterations.\n", max_iteration);
		return SUANPAN_FAIL;
	}

	if(sqrt_j2 - shear * gamma >= 0.) {
		const auto norm_s = tensor::stress::norm(dev_stress);

		const auto t_factor = shear / sqrt_j2 * gamma;

		trial_stress -= t_factor * dev_stress + bulk * eta_flow * gamma * tensor::unit_tensor2;

		trial_stiffness += double_shear * (t_factor - shear / denominator) / norm_s / norm_s * dev_stress * dev_stress.t() - double_shear * t_factor * unit_dev_tensor - factor_d / denominator * unit_x_unit;

		const mat t_mat = eta_yield * factor_c / denominator / norm_s * dev_stress * tensor::unit_tensor2.t();

		associated ? trial_stiffness -= t_mat + t_mat.t() : trial_stiffness -= t_mat + eta_flow / eta_yield * t_mat.t();
	} else {
		// apex return
		gamma = 0.; // volumetric strain reuse variable
		plastic_strain = current_history(0);

		counter = 0;
		while(++counter < max_iteration) {
			const auto residual = compute_c(plastic_strain) * xi / eta_flow - hydro_stress + bulk * gamma;
			const auto incre_gamma = residual / (denominator = factor_b * compute_dc(plastic_strain) + bulk);
			suanpan_extra_debug("NonlinearDruckerPrager local iterative loop error: %.5E.\n", incre_gamma);
			if(fabs(incre_gamma) <= tolerance) break;
			plastic_strain = current_history(0) + xi / eta_yield * (gamma -= incre_gamma);
		}

		if(max_iteration == counter) {
			suanpan_error("cannot converge in %u iterations.\n", max_iteration);
			return SUANPAN_FAIL;
		}

		trial_stress = (hydro_stress - bulk * gamma) * tensor::unit_tensor2;

		trial_stiffness = (bulk - bulk * bulk / denominator) * unit_x_unit;
	}

	return SUANPAN_SUCCESS;
}

int NonlinearDruckerPrager::clear_status() {
	current_strain.zeros();
	current_stress.zeros();
	current_history.zeros();
	current_stiffness = initial_stiffness;
	return reset_status();
}

int NonlinearDruckerPrager::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_history = trial_history;
	current_stiffness = trial_stiffness;
	return SUANPAN_SUCCESS;
}

int NonlinearDruckerPrager::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_history = current_history;
	trial_stiffness = current_stiffness;
	return SUANPAN_SUCCESS;
}

void NonlinearDruckerPrager::print() { suanpan_info("A 3D nonlinear model using Drucker--Prager yielding criterion.\n"); }
