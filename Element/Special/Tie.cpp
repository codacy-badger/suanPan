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

#include "Tie.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

const mat Tie::penalty{{1., -1.}, {-1., 1.}};
const double Tie::big_num = 1E12;

Tie::Tie(const unsigned T, const unsigned NA, const unsigned DA, const unsigned NB, const unsigned DB)
	: Element(T, 2, std::max(DA, DB), {NA, NB}, {}, {}, false)
	, t_dof(2 * std::max(DA, DB))
	, t_span{DA - 1, DB - 1 + std::max(DA, DB)} {}

void Tie::initialize(const shared_ptr<DomainBase>&) {
	initial_stiffness.zeros(t_dof, t_dof);
	initial_stiffness(t_span, t_span) = big_num * penalty;
	trial_stiffness = current_stiffness = initial_stiffness;
}

int Tie::update_status() {
	const auto t_disp = get_trial_displacement();

	trial_stiffness(t_span, t_span) = big_num / std::max(fabs(t_disp(t_span(0)) - t_disp(t_span(1))), 1.) * penalty;

	trial_resistance = trial_stiffness * t_disp;

	return SUANPAN_SUCCESS;
}

int Tie::clear_status() { return SUANPAN_SUCCESS; }

int Tie::commit_status() { return SUANPAN_SUCCESS; }

int Tie::reset_status() { return SUANPAN_SUCCESS; }
