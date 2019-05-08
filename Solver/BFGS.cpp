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

#include "BFGS.h"
#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>

BFGS::BFGS(const unsigned T, const unsigned MH)
	: Solver(T)
	, max_hist(MH) {}

int BFGS::analyze() {
	auto& C = get_converger();
	auto& G = get_integrator();
	const auto& D = C->get_domain().lock();
	const auto& W = D->get_factory();

	const auto& max_iteration = C->get_max_iteration();

	suanpan_info("current analysis time: %.5f.\n", W->get_trial_time());

	// iteration counter
	unsigned counter = 0;

	// ninja alias
	auto& ninja = get_ninja(W);

	// clear container
	hist_ninja.clear();
	hist_residual.clear();
	hist_factor.clear();

	while(true) {
		// assemble resistance
		G->assemble_resistance();

		// displacement load only applies once at first iteration
		// erase for following iterations
		if(counter == 1) {
			auto& t_load = get_trial_load(W);
			for(const auto& I : D->get_constrained_dof()) t_load(I) = 0.;
		}

		if(counter == 0) {
			// asemble stiffness for the first iteration
			G->assemble_matrix();
			// process loads and constraints
			if(G->process_load() != SUANPAN_SUCCESS) return SUANPAN_FAIL;
			if(G->process_constraint() != SUANPAN_SUCCESS) return SUANPAN_FAIL;
			// commit current residual
			hist_residual.emplace_back(W->get_trial_load() - W->get_sushi());
			// solve the system and commit current displacement increment
			hist_ninja.emplace_back(W->get_stiffness()->solve(*hist_residual.crbegin()));
			// copy current displacement increment to ninja
			ninja = *hist_ninja.crbegin(); // only for updating status
		} else {
			// clear temporary factor container
			alpha.clear();
			// commit current residual
			hist_residual.emplace_back(W->get_trial_load() - W->get_sushi());
			// copy current residual to ninja
			ninja = *hist_residual.crbegin();
			// perform two-step recursive loop
			const auto S = hist_factor.size() - 1; // intermediate factor
			// right side loop
			for(unsigned I = 0; I <= S; ++I) {
				const auto SI = S - I; // intermediate factor
				// compute and commit alpha
				alpha.emplace_back(dot(hist_ninja[SI], ninja) / hist_factor[SI]);
				// update ninja
				ninja -= *alpha.crbegin() * hist_residual[SI];
			}
			// apply the Hessian from the factorazation in the first iteration
			ninja = W->get_stiffness()->solve_trs(ninja);
			// left side loop
			for(unsigned I = 0; I <= S; ++I) ninja += (alpha[S - I] - dot(hist_residual[I], ninja) / hist_factor[I]) * hist_ninja[I];
			// ninja now stores current displacement increment
			hist_ninja.emplace_back(ninja);
		}
		// commit current factor after obtaining ninja and residual
		hist_factor.emplace_back(dot(*hist_ninja.crbegin(), *hist_residual.crbegin()));

		// avoid machine error accumulation
		G->erase_machine_error();
		// update trial status for factory
		W->update_trial_displacement(W->get_trial_displacement() + W->get_ninja());
		// update for nodes and elements
		if(G->update_trial_status() != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		// exit if converged
		if(C->is_converged()) return SUANPAN_SUCCESS;
		// exit if maximum iteration is hit
		if(++counter > max_iteration) return SUANPAN_FAIL;

		// check if the maximum record number is hit (L-BFGS)
		if(counter > max_hist) {
			hist_ninja.pop_front();
			hist_residual.pop_front();
			hist_factor.pop_front();
		}
	}
}

void BFGS::print() { suanpan_info("(L-)BFGS.\n"); }
