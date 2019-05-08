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

#include "AFC.h"

AFC::AFC(const unsigned T, const double E, const double TYS, const double THK, const double TUK, const double CYS, const double CHK, const double CUK, const double R)
	: DataAFC{fabs(E), fabs(TYS), fabs(THK), fabs(TUK), fabs(CYS), fabs(CHK), fabs(CUK)}
	, Material1D(T, R) {}

void AFC::initialize(const shared_ptr<DomainBase>&) {
	trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;
	trial_history = current_history.zeros(1);
}

unique_ptr<Material> AFC::get_copy() { return make_unique<AFC>(*this); }

int AFC::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	if(fabs(incre_strain(0)) <= tolerance) return SUANPAN_SUCCESS;

	const auto& i_strain = incre_strain(0);
	const auto& c_strain = current_strain(0);

	trial_history = current_history;
	auto& residual = trial_history(0);

	trial_stress = current_stress;
	auto& c_stress = trial_stress(0);

	if(c_stress >= 0.) {
		if(i_strain >= 0.) {
			c_stress += elastic_modulus * i_strain;
			const auto limit = t_yield_stress + (t_strain(0) - t_yield_strain) * t_hardening;
			if(c_stress > limit) {
				c_stress = limit;
				trial_stiffness = t_hardening;
			} else trial_stiffness = elastic_modulus;
		} else if(t_strain(0) < (residual = c_strain - c_stress / t_unloading)) {
			c_stress = (t_strain(0) - residual) * elastic_modulus;
			const auto limit = (t_strain(0) + c_yield_strain) * c_hardening - c_yield_stress;
			if(c_stress < limit) {
				c_stress = limit;
				trial_stiffness = c_hardening;
			} else trial_stiffness = elastic_modulus;
		} else {
			c_stress += t_unloading * i_strain;
			trial_stiffness = t_unloading;
		}
	} else if(c_stress < 0.) {
		if(i_strain <= 0.) {
			c_stress += elastic_modulus * i_strain;
			const auto limit = (t_strain(0) + c_yield_strain) * c_hardening - c_yield_stress;
			if(c_stress < limit) {
				c_stress = limit;
				trial_stiffness = c_hardening;
			} else trial_stiffness = elastic_modulus;
		} else if(t_strain(0) > (residual = c_strain - c_stress / c_unloading)) {
			c_stress = (t_strain(0) - residual) * elastic_modulus;
			const auto limit = t_yield_stress + (t_strain(0) - t_yield_strain) * t_hardening;
			if(c_stress > limit) {
				c_stress = limit;
				trial_stiffness = t_hardening;
			} else trial_stiffness = elastic_modulus;
		} else {
			c_stress += c_unloading * i_strain;
			trial_stiffness = c_unloading;
		}
	}

	return SUANPAN_SUCCESS;
}

int AFC::clear_status() {
	current_strain.zeros();
	current_stress.zeros();
	current_stiffness = initial_stiffness;
	current_history.zeros();
	return reset_status();
}

int AFC::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_stiffness = trial_stiffness;
	current_history = trial_history;
	return SUANPAN_SUCCESS;
}

int AFC::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_stiffness = current_stiffness;
	trial_history = current_history;
	return SUANPAN_SUCCESS;
}

void AFC::print() { suanpan_info("A AFC material model.\n"); }
