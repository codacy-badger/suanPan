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

#include "ABM4.h"
#include "ODE.h"

ABM4::ABM4(const unsigned T, const unsigned N, const bool C)
	: ODE_Implicit(T, N, 1. / 8., C) {}

int ABM4::analyze(ODE* ode_system) {
	history_step.clear();

	const auto step_size = ode_system->get_incre_time() / double(num_step);

	auto c_time = ode_system->get_current_time();
	auto c_disp = ode_system->get_current_variable();

	history_step.emplace_back(ode_system->eval(c_time, c_disp));

	c_time += step_size;

	history_step.emplace_back(ode_system->eval(c_time, c_disp += .5 * step_size * (history_step[0] + ode_system->eval(c_time, c_disp + step_size * history_step[0]))));

	history_step.emplace_back(ode_system->eval(c_time += step_size, c_disp += step_size * (1.5 * history_step[1] - .5 * history_step[0])));

	history_step.emplace_back(ode_system->eval(c_time += step_size, c_disp += step_size / 12. * (23. * history_step[2] - 16. * history_step[1] + 5. * history_step[0])));

	vec i_disp;

	unsigned counter = 3;

	while(++counter <= num_step) {
		i_disp = step_size / 24. * (55. * history_step[3] - 59. * history_step[2] + 37. * history_step[1] - 9. * history_step[0]);
		history_step.emplace_back(ode_system->eval(c_time += step_size, c_disp += use_corrector ? step_size / 24. * (9. * ode_system->eval(c_time + step_size, c_disp + i_disp) + 19. * history_step[3] - 5. * history_step[2] + history_step[1]) : i_disp));
		history_step.pop_front();
	}

	ode_system->set_trial_variable(c_disp);

	ode_system->set_error(norm(i_disp));

	return 0;
}

void ABM4::print() { suanpan_info("A Bogacki--Shampine 2/3 pair solver object.\n"); }
