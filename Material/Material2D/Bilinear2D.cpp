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

#include "Bilinear2D.h"
#include <Recorder/OutputType.h>
#include <Toolbox/tensorToolbox.h>

const uvec Bilinear2D::F1{0, 1, 3};
const uvec Bilinear2D::F2{2, 4, 5};

mat Bilinear2D::form_stiffness(const mat& full_stiffness) { return full_stiffness(F1, F1) - full_stiffness(F1, F2) * solve(full_stiffness(F2, F2), full_stiffness(F2, F1)); }

Bilinear2D::Bilinear2D(const unsigned T, const double E, const double V, const double Y, const double H, const double B, const PlaneType M, const double D)
	: Material2D(T, M, D)
	, elastic_modulus(E)
	, poissons_ratio(V)
	, base(0, E, V, Y, H, B, D) {}

void Bilinear2D::initialize(const shared_ptr<DomainBase>&) {
	base.Material::initialize();
	base.initialize();

	trial_full_strain = current_full_strain.zeros(6);

	trial_stiffness = current_stiffness = initial_stiffness = form_stiffness(base.get_initial_stiffness());
}

double Bilinear2D::get_parameter(const ParameterType T) const { return base.get_parameter(T); }

unique_ptr<Material> Bilinear2D::get_copy() { return make_unique<Bilinear2D>(*this); }

int Bilinear2D::update_trial_status(const vec& t_strain) {
	auto& t_stress = base.get_trial_stress();
	auto& t_stiffness = base.get_trial_stiffness();

	if(PlaneType::S == plane_type) {
		incre_strain = t_strain - trial_strain;
		trial_full_strain(F1) = trial_strain = t_strain;

		trial_full_strain(F2) -= solve(t_stiffness(F2, F2), t_stress(F2) + t_stiffness(F2, F1) * incre_strain);

		if(base.update_trial_status(trial_full_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		trial_stiffness = form_stiffness(t_stiffness);

		trial_stress = t_stress(F1) - t_stiffness(F1, F2) * solve(t_stiffness(F2, F2), t_stress(F2));
	} else {
		trial_full_strain(F1) = trial_strain = t_strain;

		if(base.update_trial_status(trial_full_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		trial_stress = t_stress(F1);
		trial_stiffness = t_stiffness(F1, F1);
	}

	return SUANPAN_SUCCESS;
}

int Bilinear2D::clear_status() {
	current_strain.zeros();
	current_stress.zeros();
	current_full_strain.zeros();
	trial_strain.zeros();
	trial_stress.zeros();
	trial_full_strain.zeros();
	trial_stiffness = current_stiffness = initial_stiffness;

	return base.clear_status();
}

int Bilinear2D::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_full_strain = trial_full_strain;
	current_stiffness = trial_stiffness;

	return base.commit_status();
}

int Bilinear2D::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_full_strain = current_full_strain;
	trial_stiffness = current_stiffness;

	return base.reset_status();
}

void Bilinear2D::print() {
	suanpan_info("A 2D bilinear plane %s material model.\n", plane_type == PlaneType::S ? "stress" : "strain");
	suanpan_info("current strain: ");
	current_strain.t().print();
	suanpan_info("current stress: ");
	current_stress.t().print();
}

vector<vec> Bilinear2D::record(const OutputType P) {
	vector<vec> output;
	output.reserve(1);

	if(P == OutputType::PE) output.emplace_back(current_strain - solve(initial_stiffness, current_stress));
	else if(P == OutputType::PEP) output.emplace_back(transform::strain::principal(current_strain - solve(initial_stiffness, current_stress)));
	else if(P == OutputType::MISES) {
		vec trial_mises(1);
		if(plane_type == PlaneType::S) trial_mises(0) = sqrt(current_stress(0) * current_stress(0) - current_stress(0) * current_stress(1) + current_stress(1) * current_stress(1) + 3. * current_stress(2) * current_stress(2));
		else if(plane_type == PlaneType::E) {
			const auto sigma_33 = elastic_modulus * poissons_ratio / (1. + poissons_ratio) / (1. - 2. * poissons_ratio) * (current_strain(0) + current_strain(1));
			const auto sigma_mean = (current_stress(0) + current_stress(1) + sigma_33) / 3.;
			const auto tmp_a = current_stress(0) - sigma_mean;
			const auto tmp_b = current_stress(1) - sigma_mean;
			const auto tmp_c = sigma_33 - sigma_mean;
			trial_mises(0) = sqrt(1.5 * (tmp_a * tmp_a + tmp_b * tmp_b + tmp_c * tmp_c + 2. * current_stress(2) * current_stress(2)));
		}
		output.emplace_back(trial_mises);
	} else if(P == OutputType::EEEQ) output.emplace_back(vec{sqrt(2. / 3.) * tensor::strain::norm(current_full_strain)});
	else return Material2D::record(P);

	return output;
}
