/*******************************************************************************
 * Copyright (C) 2017-2019 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
/**
 * @class SparseSymmMat
 * @brief A SparseSymmMat class that holds matrices.
 *
 * @author tlc
 * @date 14/01/2019
 * @version 0.1.0
 * @file SparseSymmMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSSYMMEMAT_HPP
#define SPARSSYMMEMAT_HPP

template<typename T> class SparseSymmMat final : public SparseMat<T> {
	superlu_opts symm_opts;
public:
	using SparseMat<T>::SparseMat;
	using SparseMat<T>::solve;

	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> int SparseSymmMat<T>::solve(Mat<T>& out_mat, const Mat<T>& in_mat) {
#ifdef SUANPAN_MAGMA
	SparseMat<T>::solve(out_mat, in_mat);
	return 0;
#else
	this->triplet_to_csc();

	const uvec row_idx(this->csc_mat.row_idx, this->csc_mat.c_size, false, false);
	const uvec col_ptr(this->csc_mat.col_ptr, this->csc_mat.n_cols + 1, false, false);
	const Col<T> val_idx(this->csc_mat.val_idx, this->csc_mat.c_size, false, false);

	this->arma_mat = SpMat<T>(row_idx, col_ptr, val_idx, this->csc_mat.n_rows, this->csc_mat.n_cols);

	symm_opts.symmetric = true;
	return spsolve(out_mat, this->arma_mat, in_mat, "s", symm_opts) ? 0 : -1;
#endif
}

#endif

//! @}
