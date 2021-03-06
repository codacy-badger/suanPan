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

#include "Flag.h"
#include <Toolbox/utility.h>

Flag::Flag(const unsigned T, const double E, const double YT, const double RT, const double HT, const double YC, const double RC, const double HC, const double D)
	: DataFlag{fabs(E), HT, fabs(YT), fabs(RT), HC, fabs(YC), fabs(RC)}
	, Material1D(T, D) {}

Flag::Flag(const unsigned T, const double E, const double YT, const double RT, const double HT, const double D)
	: DataFlag{fabs(E), HT, fabs(YT), fabs(RT), HT, fabs(YT), fabs(RT)}
	, Material1D(T, D) {}

void Flag::initialize(const shared_ptr<DomainBase>&) {
	trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

	trial_history = current_history.zeros(4);
}

unique_ptr<Material> Flag::get_copy() { return make_unique<Flag>(*this); }

int Flag::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

	trial_status = current_status;
	trial_history = current_history;
	auto& tr_strain = trial_history(0);
	auto& tr_low_strain = trial_history(1);
	auto& cr_strain = trial_history(2);
	auto& cr_low_strain = trial_history(3);

	const auto load_direction = suanpan::sign(incre_strain(0));

	switch(trial_status) {
	case Status::NONE:
		trial_status = load_direction > 0. ? Status::TLOAD : Status::CLOAD;
		break;
	case Status::TLOAD:
		if(load_direction < 0.) {
			tr_low_strain = ((current_stress(0) - t_residual_stress) / elastic_modulus + t_residual_strain * t_hardening_ratio - (tr_strain = current_strain(0))) / (t_hardening_ratio - 1.);
			trial_status = trial_strain(0) > tr_low_strain ? Status::TUNLOAD : trial_strain(0) > t_residual_strain ? Status::TLOW : trial_strain(0) > 0. ? Status::TLOAD : Status::CLOAD;
		}
		break;
	case Status::CLOAD:
		if(load_direction > 0.) {
			cr_low_strain = ((current_stress(0) + c_residual_stress) / elastic_modulus - c_residual_strain * c_hardening_ratio - (cr_strain = current_strain(0))) / (c_hardening_ratio - 1.);
			trial_status = trial_strain(0) < cr_low_strain ? Status::CUNLOAD : trial_strain(0) < -c_residual_strain ? Status::CLOW : trial_strain(0) < 0. ? Status::CLOAD : Status::TLOAD;
		}
		break;
	case Status::TLOW:
		if(load_direction > 0.) {
			tr_strain = (tr_low_strain = current_strain(0)) + t_yield_strain - t_residual_strain;
			trial_status = trial_strain(0) < tr_strain ? Status::TUNLOAD : Status::TLOAD;
		} else trial_status = trial_strain(0) > t_residual_strain ? Status::TLOW : trial_strain(0) > 0. ? Status::TLOAD : Status::CLOAD;
		break;
	case Status::CLOW:
		if(load_direction < 0.) {
			cr_strain = (cr_low_strain = current_strain(0)) - c_yield_strain + c_residual_strain;
			trial_status = trial_strain(0) > cr_strain ? Status::CUNLOAD : Status::CLOAD;
		} else trial_status = trial_strain(0) < -c_residual_strain ? Status::CLOW : trial_strain(0) < 0. ? Status::CLOAD : Status::TLOAD;
		break;
	case Status::TUNLOAD:
		trial_status = trial_strain(0) > tr_strain ? Status::TLOAD : trial_strain(0) > tr_low_strain ? Status::TUNLOAD : trial_strain(0) > t_residual_strain ? Status::TLOW : trial_strain(0) > 0. ? Status::TLOAD : Status::CLOAD;
		break;
	case Status::CUNLOAD:
		trial_status = trial_strain(0) < cr_strain ? Status::CLOAD : trial_strain(0) < cr_low_strain ? Status::CUNLOAD : trial_strain(0) < -c_residual_strain ? Status::CLOW : trial_strain(0) < 0. ? Status::CLOAD : Status::TLOAD;
		break;
	}

	trial_stiffness = elastic_modulus;

	if(trial_status == Status::TLOAD) trial_strain(0) > t_yield_strain ? trial_stress = t_yield_stress + (trial_stiffness *= t_hardening_ratio) * (trial_strain(0) - t_yield_strain) : trial_stress = trial_stiffness * trial_strain(0);
	else if(trial_status == Status::CLOAD) trial_strain(0) < -c_yield_strain ? trial_stress = -c_yield_stress + (trial_stiffness *= c_hardening_ratio) * (trial_strain(0) + c_yield_strain) : trial_stress = trial_stiffness * trial_strain(0);
	else if(trial_status == Status::TUNLOAD) trial_stress = current_stress + trial_stiffness * incre_strain;
	else if(trial_status == Status::CUNLOAD) trial_stress = current_stress + trial_stiffness * incre_strain;
	else if(trial_status == Status::TLOW) trial_stress = t_residual_stress + (trial_stiffness *= t_hardening_ratio) * (trial_strain - t_residual_strain);
	else if(trial_status == Status::CLOW) trial_stress = -c_residual_stress + (trial_stiffness *= c_hardening_ratio) * (trial_strain + c_residual_strain);

	return SUANPAN_SUCCESS;
}

int Flag::clear_status() {
	current_strain.zeros();
	current_stress.zeros();
	current_history.zeros();
	current_status = Status::NONE;
	current_stiffness = initial_stiffness;
	return reset_status();
}

int Flag::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_history = trial_history;
	current_status = trial_status;
	current_stiffness = trial_stiffness;
	return 0;
}

int Flag::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_history = current_history;
	trial_status = current_status;
	trial_stiffness = current_stiffness;
	return 0;
}

void Flag::print() {
	suanpan_info("A bilinear flag material model with an elastic modulus of %.3E and a hardening ratio of %.2f.\n", elastic_modulus, t_hardening_ratio);
	suanpan_info("current strain: %.3E\tcurrent stress: %.3E\n", current_strain(0), current_stress(0));
}
