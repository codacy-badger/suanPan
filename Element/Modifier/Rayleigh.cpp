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

#include "Rayleigh.h"

Rayleigh::Rayleigh(const unsigned T, const double A, const double B, const double C, const double D, uvec&& ET)
	: Modifier(T, std::forward<uvec>(ET))
	, a(A)
	, b(B)
	, c(C)
	, d(D) {}

int Rayleigh::update_status() {
	suanpan_for_each(element_pool.cbegin(), element_pool.cend(), [&](const weak_ptr<Element>& ele_ptr) {
			if(const auto t_ptr = ele_ptr.lock()) {
				mat t_damping(t_ptr->get_total_number(), t_ptr->get_total_number(), fill::zeros);
				if(a != 0. && !t_ptr->get_trial_mass().empty()) t_damping += a * t_ptr->get_trial_mass();
				if(b != 0. && !t_ptr->get_trial_stiffness().empty()) t_damping += b * t_ptr->get_trial_stiffness();
				if(c != 0. && !t_ptr->get_current_stiffness().empty()) t_damping += c * t_ptr->get_current_stiffness();
				if(d != 0. && !t_ptr->get_initial_stiffness().empty()) t_damping += d * t_ptr->get_initial_stiffness();
				access::rw(t_ptr->get_trial_damping()) = t_damping;
			}
		});

	return SUANPAN_SUCCESS;
}
