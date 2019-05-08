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

#include "Dynamic.h"
#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Load/NodalDisplacement.h>
#include <Solver/Integrator/Newmark.h>
#include <Solver/MPDC.h>
#include <Solver/Newton.h>

Dynamic::Dynamic(const unsigned T, const double P)
	: Step(T, P) {}

int Dynamic::initialize() {
	const auto& t_domain = database.lock();

	if(modifier == nullptr) modifier = make_shared<Integrator>();

	// automatically enable displacement control solver
	if(solver == nullptr) {
		auto flag = false;
		const auto step_tag = get_tag();
		for(const auto& I : t_domain->get_load_pool())
			if(typeid(*I) == typeid(NodalDisplacement) && I->get_start_step() == step_tag) {
				flag = true;
				break;
			}
		flag ? solver = make_shared<MPDC>() : solver = make_shared<Newton>();
	}

	if(sparse_mat && symm_mat) factory->set_storage_scheme(StorageScheme::SPARSESYMM);
	else if(sparse_mat && !symm_mat) factory->set_storage_scheme(StorageScheme::SPARSE);
	else if(symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BANDSYMM);
	else if(!symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BAND);
	else if(symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::SYMMPACK);
	else if(!symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::FULL);

	factory->set_analysis_type(AnalysisType::DYNAMICS);
	factory->initialize();

	tester->set_domain(t_domain);
	modifier->set_domain(t_domain);
	solver->set_converger(tester);
	solver->set_integrator(modifier);

	return t_domain->update_current_status();
}

int Dynamic::analyze() {
	get_domain().lock()->initialize_load();

	auto& W = get_factory();

	// used in displacement controlled algorithm
	// need to initialize it for updating of load factor
	auto& ref_dof = W->get_reference_dof();
	if(!ref_dof.is_empty()) {
		W->initialize_settlement();
		W->initialize_load_factor();
		get_reference_load(W).rows(ref_dof) = eye(ref_dof.n_elem, ref_dof.n_elem);
	}

	auto& S = get_solver();
	auto& G = get_integrator();

	auto time_left = get_time_period();
	auto step = get_ini_step_size();

	unsigned num_increment = 0, num_converged_step = 0;

	while(true) {
		// check if the target time point is hit
		if(time_left <= 1E-7) return SUANPAN_SUCCESS;
		// check if the maximum substep number is hit
		if(++num_increment > get_max_substep()) {
			suanpan_warning("analyze() reaches maximum substep number %u.\n", get_max_substep());
			return SUANPAN_FAIL;
		}
		// update incremental and trial time
		G->update_incre_time(step);
		// call solver
		const auto code = S->analyze();
		if(code == SUANPAN_SUCCESS) {
			// success step
			// commit converged iteration
			G->commit_status();
			// record response
			G->record();
			// eat current increment
			time_left -= step;
			if(!is_fixed_step_size() && ++num_converged_step > 5) {
				step *= 1.2;
				num_converged_step = 0;
			}
			// check if time overflows
			if(step > time_left) step = time_left;
		} else if(code == SUANPAN_FAIL) {
			// failed step
			// reset to the start of current substep
			G->reset_status();
			// check if minimum step size is hit
			if(step <= get_min_step_size()) {
				suanpan_error("analyze() reaches minimum step size %.3E.\n", get_min_step_size());
				return SUANPAN_FAIL;
			}
			// check if fixed step size
			if(is_fixed_step_size()) {
				suanpan_error("analyze() does not converge for given fixed step size %.3E.\n", step);
				return SUANPAN_FAIL;
			}
			// step size is allowed to decrease
			step *= .5;
		} else return SUANPAN_FAIL; // positive codes are from lapack subroutines
	}
}
