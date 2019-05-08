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

#include "Frequency.h"
#include <Domain/Factory.hpp>
#include <Domain/DomainBase.h>
#include <Toolbox/arpack_wrapper.h>

Frequency::Frequency(const unsigned T, const unsigned N)
	: Step(T, 0.)
	, eigen_number(N) {}

int Frequency::initialize() {
	const auto& t_domain = database.lock();

	if(sparse_mat && symm_mat) factory->set_storage_scheme(StorageScheme::SPARSESYMM);
	else if(sparse_mat && !symm_mat) factory->set_storage_scheme(StorageScheme::SPARSE);
	else if(symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BANDSYMM);
	else if(!symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BAND);
	else if(symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::SYMMPACK);
	else if(!symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::FULL);

	factory->set_analysis_type(AnalysisType::EIGEN);
	factory->initialize();

	return t_domain->update_current_status();
}

int Frequency::analyze() {
	const auto& D = database.lock();
	auto& W = D->get_factory();

	D->initialize_load();
	D->assemble_mass();
	D->assemble_stiffness();

	if(D->process_load() != SUANPAN_SUCCESS || D->process_constraint() != SUANPAN_SUCCESS) return SUANPAN_FAIL;

	if(eig_solve(get_eigenvalue(W), get_eigenvector(W), W->get_stiffness(), W->get_mass(), eigen_number, "LM") != SUANPAN_SUCCESS) {
		suanpan_warning("fail to decompose the system, try to increase the number of eigen values.\n");
		return SUANPAN_SUCCESS;
	}

	D->record();

	return SUANPAN_SUCCESS;
}

void Frequency::set_eigen_number(const unsigned N) const { access::rw(eigen_number) = N; }

unsigned Frequency::get_eigen_number() const { return eigen_number; }
