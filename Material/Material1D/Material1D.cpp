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

#include "Material1D.h"
#include <Recorder/OutputType.h>
#include <Toolbox/utility.h>

Material1D::Material1D(const unsigned T, const double D)
	: Material(T, MaterialType::D1, D) {}

const mat& Material1D::get_trial_stiffness() {
	// for exact zero case do not modify as the material may have no stiffness
	// for other cases using secant stiffness first
	// if still very small keep stiffness as a nonzero small number
	if(trial_stiffness(0) != 0. && fabs(trial_stiffness(0)) < 1E-8) {
		const auto denominator = trial_strain(0) - current_strain(0);
		if(fabs(denominator) > 1E-10) {
			const auto secant_stiffness = (trial_stress(0) - current_stress(0)) / denominator;
			if(fabs(secant_stiffness) < fabs(trial_stiffness(0))) trial_stiffness = secant_stiffness;
		}
		if(fabs(trial_stiffness(0)) < 1E-8) {
			const auto stiffness_sign = suanpan::sign(trial_stiffness(0));
			trial_stiffness(0) = (stiffness_sign == 0. ? 1. : stiffness_sign) * 1E-8;
		}
	}

	return trial_stiffness;
}

vector<vec> Material1D::record(const OutputType P) {
	vector<vec> data;

	if(P == OutputType::S) data.emplace_back(current_stress);
	else if(P == OutputType::E) data.emplace_back(current_strain);
	else if(P == OutputType::S11) data.emplace_back(current_stress);
	else if(P == OutputType::E11) data.emplace_back(current_strain);
	else if(P == OutputType::EE) data.emplace_back(current_stress / initial_stiffness);
	else if(P == OutputType::PE) data.emplace_back(current_strain - current_stress / initial_stiffness);

	return data;
}
