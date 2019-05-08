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

#include "FrameRecorder.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Element/Element.h>

#ifdef SUANPAN_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

FrameRecorder::FrameRecorder(const unsigned T, const OutputType L, const unsigned I)
	: Recorder(T, {}, L, I, false, true) {
#ifdef SUANPAN_HDF5
	ostringstream file_name;
	file_name << to_char(get_variable_type());
	file_name << ".h5";
	file_id = H5Fcreate(file_name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#endif
}

FrameRecorder::~FrameRecorder() {
#ifdef SUANPAN_HDF5
	H5Fclose(file_id);
#endif
}

void FrameRecorder::record(const shared_ptr<DomainBase>& D) {
#ifdef SUANPAN_HDF5
	if(1 == interval || counter++ == interval) {
		counter = 1;

		ostringstream group_name;
		group_name << "/";
		group_name << D->get_factory()->get_current_time();

		const auto group_id = H5Gcreate(file_id, group_name.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		for(const auto& I : D->get_element_pool()) {
			const auto data = I->record(get_variable_type());

			if(!data.empty()) {
				mat data_to_write(data[0].n_elem, data.size());

				uword idx = 0;
				for(const auto& J : data) data_to_write.col(idx++) = J;

				hsize_t dimention[2] = {data_to_write.n_cols, data_to_write.n_rows};

				H5LTmake_dataset(group_id, std::to_string(I->get_tag()).c_str(), 2, dimention, H5T_NATIVE_DOUBLE, data_to_write.mem);
			}
		}

		H5Gclose(group_id);
	}
#endif
}

void FrameRecorder::save() {}

void FrameRecorder::print() { suanpan_info("An Element Recorder.\n"); }