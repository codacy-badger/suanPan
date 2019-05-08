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

#include "ArcLength.h"
#include <Domain/Domain.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Solver/Integrator/Integrator.h>
#include <Solver/Ramm.h>
#include <Converger/Converger.h>

ArcLength::ArcLength(const unsigned T, const unsigned NT, const unsigned DT, const double MA)
	: Step(T, 0.)
	, node(NT)
	, dof(DT)
	, maginitude(MA) {}

int ArcLength::initialize() {
	const auto& t_domain = database.lock();

	if(solver == nullptr) solver = make_shared<Ramm>();
	if(modifier == nullptr) modifier = make_shared<Integrator>();

	if(sparse_mat && symm_mat) factory->set_storage_scheme(StorageScheme::SPARSESYMM);
	else if(sparse_mat && !symm_mat) factory->set_storage_scheme(StorageScheme::SPARSE);
	else if(symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BANDSYMM);
	else if(!symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BAND);
	else if(symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::SYMMPACK);
	else if(!symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::FULL);

	factory->set_analysis_type(AnalysisType::STATICS);
	factory->initialize();

	tester->set_domain(t_domain);
	modifier->set_domain(t_domain);
	solver->set_converger(tester);
	solver->set_integrator(modifier);

	factory->set_reference_size(1);
	factory->initialize_load_factor();

	get_reference_load(factory)(t_domain->get_node(node)->get_reordered_dof().at(dof - 1)) = maginitude;

	return t_domain->update_current_status();
}

int ArcLength::analyze() {
	get_domain().lock()->initialize_load();

	auto& S = get_solver();
	auto& G = get_integrator();

	unsigned num_iteration = 0;

	while(true) {
		if(num_iteration++ > get_max_substep()) {
			suanpan_warning("analyze() reaches maximum substep number %u.\n", get_max_substep());
			return SUANPAN_FAIL;
		}
		auto code = G->process_criterion();
		if(code != SUANPAN_SUCCESS) return code;
		code = S->analyze();
		if(code == SUANPAN_SUCCESS) {
			G->commit_status();
			G->record();
		} else if(code == SUANPAN_FAIL) G->reset_status();
		else return SUANPAN_FAIL;
	}
}
