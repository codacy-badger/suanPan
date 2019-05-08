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

#include "Viscosity01.h"
#include <Toolbox/utility.h>
#include <Recorder/OutputType.h>

Viscosity01::Viscosity01(const unsigned T, const double A, const double C, const double L)
	: DataViscosity01{fabs(A), fabs(C), fabs(L)}
	, Material1D(T, 0.) {}

void Viscosity01::initialize(const shared_ptr<DomainBase>&) {
	trial_damping = current_damping = initial_damping = 1. == alpha ? damping : damping * b;

	trial_strain_rate = current_strain_rate.zeros(1);
}

unique_ptr<Material> Viscosity01::get_copy() { return make_unique<Viscosity01>(*this); }

int Viscosity01::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
	incre_strain = (trial_strain = t_strain) - current_strain;
	incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

	if(norm(incre_strain) <= datum::eps && norm(incre_strain_rate) <= datum::eps) return SUANPAN_SUCCESS;

	const auto& v = trial_strain_rate(0);
	const auto abs_v = fabs(v);

	if(1. == alpha) trial_stress = (trial_damping = damping) * v;
	else if(1. < alpha || abs_v > limit) {
		const auto pow_term = damping * pow(abs_v, alpha - 1.);
		trial_stress = v * pow_term;
		trial_damping = alpha * pow_term;
	} else {
		trial_stress = damping * (a * abs_v + b) * v;
		trial_damping = damping * (2. * a * abs_v + b);
	}

	return SUANPAN_SUCCESS;
}

int Viscosity01::clear_status() {
	current_strain.zeros();
	current_strain_rate.zeros();
	current_stress.zeros();
	current_damping = initial_damping;
	return reset_status();
}

int Viscosity01::commit_status() {
	current_strain = trial_strain;
	current_strain_rate = trial_strain_rate;
	current_stress = trial_stress;
	current_damping = trial_damping;
	return SUANPAN_SUCCESS;
}

int Viscosity01::reset_status() {
	trial_strain = current_strain;
	trial_strain_rate = current_strain_rate;
	trial_stress = current_stress;
	trial_damping = current_damping;
	return SUANPAN_SUCCESS;
}

vector<vec> Viscosity01::record(const OutputType P) {
	vector<vec> data;

	if(OutputType::S == P) data.emplace_back(current_stress);
	else if(OutputType::E == P) data.emplace_back(current_strain);
	else if(OutputType::V == P) data.emplace_back(current_strain_rate);

	return data;
}

void Viscosity01::print() { suanpan_info("A 1D vicosity material %u.\n", get_tag()); }
