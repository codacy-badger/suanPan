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

#include "AFCN.h"
#include <Toolbox/utility.h>

podarray<double> AFCN::compute_transition(const double TX, const double XS, const double YS, const double ES, const double XF, const double YF, const double EF) {
	podarray<double> response(2);

	if(fabs(TX - XS) <= datum::eps) {
		response(0) = YS;
		response(1) = ES;
		return response;
	}

	if(fabs(TX - XF) <= datum::eps) {
		response(0) = YF;
		response(1) = EF;
		return response;
	}

	const auto TA = XF - XS;
	const auto TC = TX - XS;
	const auto ESEC = (YF - YS) / TA;
	const auto TB = ESEC - ES;
	const auto R = (EF - ESEC) / TB;
	const auto TD = TB * pow(std::max(fabs(TC / TA), datum::eps), R);

	response(0) = YS + TC * (ES + TD);
	response(1) = ES + (R + 1.) * TD;

	suanpan_debug([&]() { if(!std::isfinite(response(0)) || !std::isfinite(response(1))) suanpan_fatal("infinite numbers detected.\n"); });

	return response;
}

void AFCN::compute_degradation(const double yield_strain, const double stiffness) {
	auto& max_strain = trial_history(1);
	auto& s_strain = trial_history(3);
	auto& s_stress = trial_history(4);
	auto& e_strain = trial_history(5); // limit strain
	auto& e_stress = trial_history(6); // limit stress

	max_strain = std::max(std::max(fabs(trial_strain(0)), fabs(e_strain)), max_strain);

	const auto factor = .9 * exp(degrade * (1. - max_strain / yield_strain));
	const auto s_stiffness = factor * elastic_modulus;
	const auto e_stiffness = stiffness - factor * (stiffness - elastic_modulus);

	const auto response = compute_transition(trial_strain(0), s_strain, s_stress, s_stiffness, e_strain, e_stress, e_stiffness);

	trial_stress = response(0);
	trial_stiffness = response(1);
}

AFCN::AFCN(const unsigned T, const double E, const double TYS, const double THK, const double TUK, const double CYS, const double CHK, const double CUK, const double R)
	: DataAFCN{fabs(E), fabs(TYS), fabs(THK), fabs(TUK), fabs(CYS), fabs(CHK), fabs(CUK)}
	, Material1D(T, R) {}

void AFCN::initialize(const shared_ptr<DomainBase>&) {
	trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;
	trial_history = current_history.zeros(7);
}

unique_ptr<Material> AFCN::get_copy() { return make_unique<AFCN>(*this); }

int AFCN::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	const auto& i_strain = incre_strain(0);

	if(fabs(i_strain) <= tolerance) return SUANPAN_SUCCESS;

	trial_history = current_history;
	auto& load_sign = trial_history(0);
	auto& r_strain = trial_history(2);
	auto& s_strain = trial_history(3);
	auto& s_stress = trial_history(4);
	auto& e_strain = trial_history(5); // limit strain
	auto& e_stress = trial_history(6); // limit stress

	if(current_stress(0) >= 0.)
		if(i_strain >= 0.) {
			if(load_sign <= 0.) {
				s_strain = current_strain(0), s_stress = current_stress(0);
				e_strain = s_strain + ((s_strain - t_yield_strain) * t_hardening + t_yield_stress - s_stress) / (elastic_modulus - t_hardening);
				e_stress = (e_strain - t_yield_strain) * t_hardening + t_yield_stress;
			}
			if(trial_strain(0) < e_strain) compute_degradation(t_yield_strain, t_unloading);
			else trial_stress = (trial_strain(0) - e_strain) * (trial_stiffness = t_hardening) + e_stress;
		} else {
			if(load_sign > 0.) r_strain = current_strain(0) - current_stress(0) / t_unloading;
			if(trial_strain(0) >= r_strain) trial_stress = current_stress + i_strain * (trial_stiffness = t_unloading);
			else {
				s_strain = r_strain, s_stress = 0.;
				e_strain = s_strain + ((s_strain + c_yield_strain) * c_hardening - c_yield_stress) / (elastic_modulus - c_hardening);
				e_stress = (e_strain + c_yield_strain) * c_hardening - c_yield_stress;
				if(trial_strain(0) > e_strain) compute_degradation(c_yield_strain, c_unloading);
				else trial_stress = (trial_strain(0) - e_strain) * (trial_stiffness = c_hardening) + e_stress;
			}
		}
	else if(i_strain <= 0.) {
		if(load_sign >= 0.) {
			s_strain = current_strain(0), s_stress = current_stress(0);
			e_strain = s_strain + ((s_strain + c_yield_strain) * c_hardening - c_yield_stress - s_stress) / (elastic_modulus - c_hardening);
			e_stress = (e_strain + c_yield_strain) * c_hardening - c_yield_stress;
		}
		if(trial_strain(0) > e_strain) compute_degradation(c_yield_strain, c_unloading);
		else trial_stress = (trial_strain - e_strain) * (trial_stiffness = c_hardening) + e_stress;
	} else {
		if(load_sign < 0.) r_strain = current_strain(0) - current_stress(0) / c_unloading;
		if(trial_strain(0) <= r_strain) trial_stress = current_stress + i_strain * (trial_stiffness = c_unloading);
		else {
			s_strain = r_strain, s_stress = 0.;
			e_strain = s_strain + ((s_strain - t_yield_strain) * t_hardening + t_yield_stress) / (elastic_modulus - t_hardening);
			e_stress = (e_strain - t_yield_strain) * t_hardening + t_yield_stress;
			if(trial_strain(0) < e_strain) compute_degradation(t_yield_strain, t_unloading);
			else trial_stress = (trial_strain - e_strain) * (trial_stiffness = t_hardening) + e_stress;
		}
	}

	load_sign = suanpan::sign(i_strain);

	return SUANPAN_SUCCESS;
}

int AFCN::clear_status() {
	current_strain.zeros();
	current_stress.zeros();
	current_stiffness = initial_stiffness;
	current_history.zeros();
	return reset_status();
}

int AFCN::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_stiffness = trial_stiffness;
	current_history = trial_history;
	return SUANPAN_SUCCESS;
}

int AFCN::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_stiffness = current_stiffness;
	trial_history = current_history;
	return SUANPAN_SUCCESS;
}

void AFCN::print() { suanpan_info("A AFC material model using nonlinear transition.\n"); }