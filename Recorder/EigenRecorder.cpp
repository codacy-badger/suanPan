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

#include "EigenRecorder.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

#ifdef SUANPAN_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

EigenRecorder::EigenRecorder(const unsigned T, const bool H)
	: Recorder(T, {}, OutputType::NL, 1, false, H) {}

void EigenRecorder::record(const shared_ptr<DomainBase>& D) {
	auto& W = D->get_factory();

	if(W->get_analysis_type() != AnalysisType::EIGEN) return;

	eigen_value = W->get_eigenvalue();
	const auto& eig_vector = W->get_eigenvector();

	eigen_pool = vector<std::map<unsigned, vec>>(eigen_value.n_elem);

	const auto& node_pool = D->get_node_pool();

	for(uword I = 0; I < eigen_value.n_elem; ++I) {
		const vec t_eigen_vector = eig_vector.col(I);
		for(const auto& J : node_pool) eigen_pool[I].emplace(J->get_tag(), t_eigen_vector(J->get_reordered_dof()));
	}
}

void EigenRecorder::save() {
	if(eigen_value.is_empty()) return;

#ifdef SUANPAN_HDF5
	if(if_hdf5()) {
		string file_name = "Eigenvalue.h5";

		const auto file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		hsize_t dimention[2] = {eigen_value.n_elem, 1};
		H5LTmake_dataset(file_id, "Eigenvalue", 2, dimention, H5T_NATIVE_DOUBLE, eigen_value.mem);

		for(uword I = 0; I < eigen_value.n_elem; ++I) {
			const auto group_name = "Eigenvalue " + std::to_string(I + 1);
			const auto group_id = H5Gcreate(file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

			for(const auto& J : eigen_pool[I]) {
				const auto dataset_name = "N" + std::to_string(J.first);
				dimention[0] = J.second.n_cols;
				dimention[1] = J.second.n_rows;
				H5LTmake_dataset(group_id, dataset_name.c_str(), 2, dimention, H5T_NATIVE_DOUBLE, J.second.mem);
			}

			H5Gclose(group_id);
		}

		H5Fclose(file_id);
	}
#endif
}

void EigenRecorder::print() { suanpan_info("A recorder to record eigen values and eigen vectors.\n"); }
