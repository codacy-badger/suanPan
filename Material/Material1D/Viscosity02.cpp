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

#include "Viscosity02.h"
#include <Toolbox/utility.h>
#include <Recorder/OutputType.h>

double Viscosity02::compute_dstrain(const double strain, const double strain_rate) const { return (damping_b + damping_d * atan(gap_b * strain_rate)) * gap_a / (1. + pow(gap_a * strain, 2.)); }

double Viscosity02::compute_dstrainrate(const double strain, const double strain_rate) const { return (damping_c + damping_d * atan(gap_a * strain)) * gap_b / (1. + pow(gap_b * strain_rate, 2.)); }

double Viscosity02::compute_damping_coefficient(const double strain, const double strain_rate) const {
	const auto factor_a = atan(gap_a * strain);
	const auto factor_b = atan(gap_b * strain_rate);
	return damping_a + damping_b * factor_a + damping_c * factor_b + damping_d * factor_a * factor_b;
}

Viscosity02::Viscosity02(const unsigned T, const double A, const double CA, const double CB, const double CC, const double CD, const double GA, const double GB, const double L)
	: Material1D(T, 0.)
	, alpha(fabs(A))
	, damping_a(.25 * (fabs(CA) + fabs(CB) + fabs(CC) + fabs(CD)))
	, damping_b(.5 / datum::pi * (fabs(CA) - fabs(CB) - fabs(CC) + fabs(CD)))
	, damping_c(.5 / datum::pi * (fabs(CA) + fabs(CB) - fabs(CC) - fabs(CD)))
	, damping_d((fabs(CA) - fabs(CB) + fabs(CC) - fabs(CD)) / datum::pi / datum::pi)
	, gap_a(fabs(GA))
	, gap_b(fabs(GB))
	, limit(fabs(L)) {}

void Viscosity02::initialize(const shared_ptr<DomainBase>&) {
	trial_damping = current_damping = initial_damping = 1. == alpha ? damping_a : damping_a * b;

	trial_stiffness = current_stiffness = initial_stiffness.zeros(1);

	trial_strain_rate = current_strain_rate.zeros(1);
}

unique_ptr<Material> Viscosity02::get_copy() { return make_unique<Viscosity02>(*this); }

int Viscosity02::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
	incre_strain = (trial_strain = t_strain) - current_strain;
	incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

	if(norm(incre_strain) <= datum::eps && norm(incre_strain_rate) <= datum::eps) return SUANPAN_SUCCESS;

	const auto &u = trial_strain(0), &v = trial_strain_rate(0);

	const auto abs_v = fabs(v);

	const auto eta = compute_damping_coefficient(u, v);

	double term_a, term_b;
	if(1. == alpha) {
		term_a = 1.;
		term_b = v;
	} else if(1. < alpha || limit < abs_v) {
		term_a = alpha * pow(abs_v, alpha - 1.);
		term_b = term_a * v / alpha;
	} else {
		term_a = 2. * a * abs_v + b;
		term_b = (a * abs_v + b) * v;
	}

	trial_stress = eta * term_b;
	trial_stiffness = compute_dstrain(u, v) * term_b;
	trial_damping = eta * term_a + compute_dstrainrate(u, v) * term_b;

	return SUANPAN_SUCCESS;
}

int Viscosity02::clear_status() {
	current_strain.zeros();
	current_strain_rate.zeros();
	current_stress.zeros();
	current_damping = initial_damping;
	current_stiffness = initial_stiffness;
	return reset_status();
}

int Viscosity02::commit_status() {
	current_strain = trial_strain;
	current_strain_rate = trial_strain_rate;
	current_stress = trial_stress;
	current_damping = trial_damping;
	current_stiffness = trial_stiffness;
	return SUANPAN_SUCCESS;
}

int Viscosity02::reset_status() {
	trial_strain = current_strain;
	trial_strain_rate = current_strain_rate;
	trial_stress = current_stress;
	trial_damping = current_damping;
	trial_stiffness = current_stiffness;
	return SUANPAN_SUCCESS;
}

vector<vec> Viscosity02::record(const OutputType P) {
	vector<vec> data;

	if(OutputType::S == P) data.emplace_back(current_stress);
	else if(OutputType::E == P) data.emplace_back(current_strain);
	else if(OutputType::V == P) data.emplace_back(current_strain_rate);

	return data;
}

void Viscosity02::print() { suanpan_info("A 1D vicosous damping material %u.\n", get_tag()); }
