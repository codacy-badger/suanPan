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
 * @class SparseMat
 * @brief A SparseMat class that holds matrices.
 *
 * @author tlc
 * @date 06/05/2018
 * @version 0.1.0
 * @file SparseMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSEMAT_HPP
#define SPARSEMAT_HPP

/*
#ifndef SUANPAN_MAGMA
#include <dmumps_c.h>
#endif
*/

template<typename T> class SparseMat final : public MetaMat<T> {
	superlu_opts symm_opts;
#ifdef SUANPAN_MAGMA
	magma_dopts option;
	magma_queue_t queue;
#endif
public:
	using MetaMat<T>::triplet_mat;
	using MetaMat<T>::csc_mat;
	using MetaMat<T>::arma_mat;

	SparseMat();
	SparseMat(uword, uword);
	~SparseMat();

	void unify(uword) override;

	void init_magma();

	bool is_empty() const override;

	void zeros() override;
	void reset() override;
	T max() const override;

	const T& operator()(uword, uword) const override;
	T& at(uword, uword) override;

	void triplet_to_csc();

	const T* memptr() const override { return nullptr; }

	T* memptr() override { throw invalid_argument("not supported"); }

	MetaMat<T> operator+(const MetaMat<T>&) override;
	MetaMat<T> operator-(const MetaMat<T>&) override;
	MetaMat<T>& operator+=(const MetaMat<T>&) override;
	MetaMat<T>& operator-=(const MetaMat<T>&) override;

	MetaMat<T> operator*(T) override;
	Mat<T> operator*(const Mat<T>&) override;

	MetaMat<T>& operator*=(T) override;

	Mat<T> solve(const Mat<T>&) override;
	int solve(Mat<T>&, const Mat<T>&) override;

	Mat<T> solve_trs(const Mat<T>&) override;
	int solve_trs(Mat<T>&, const Mat<T>&) override;

	MetaMat<T> factorize() override { throw invalid_argument("not supported"); }

	MetaMat<T> i() override { throw invalid_argument("not supported"); }

	MetaMat<T> inv() override { throw invalid_argument("not supported"); }

	void print() override;
	void save(const char*) override;
};

template<typename T> SparseMat<T>::SparseMat() { init_magma(); }

template<typename T> SparseMat<T>::SparseMat(const uword in_row, const uword in_col)
	: MetaMat<T>(in_row, in_col, 0) { init_magma(); }

template<typename T> SparseMat<T>::~SparseMat() {
#ifdef SUANPAN_MAGMA
	magma_queue_destroy(queue);
#endif
}

template<typename T> void SparseMat<T>::unify(const uword idx) {
	const auto t_idx = index_t(idx);
	for(index_t I = 0; I < triplet_mat.c_size; ++I) if(triplet_mat.row_idx[I] == t_idx || triplet_mat.col_idx[I] == t_idx) triplet_mat.val_idx[I] = 0.;
	triplet_mat.at(t_idx, t_idx) = 1.;
}

template<typename T> void SparseMat<T>::init_magma() {
#ifdef SUANPAN_MAGMA
	magma_queue_create(0, &queue);

	option.alignment = 1;
	option.blocksize = 32;
	option.input_format = Magma_CSR;
	option.input_location = Magma_CPU;
	option.output_format = Magma_CSR;
	option.output_location = Magma_CPU;
	option.precond_par.atol = 1E-10;
	option.precond_par.levels = 0;
	option.precond_par.maxiter = 500;
	option.precond_par.pattern = 1;
	option.precond_par.restart = 100;
	option.precond_par.rtol = 1E-14;
	option.precond_par.solver = Magma_ILU;
	option.precond_par.sweeps = 5;
	option.precond_par.trisolver = Magma_CUSOLVE;
	option.scaling = Magma_NOSCALE;
	option.solver_par.atol = 1E-10;
	option.solver_par.maxiter = 1000;
	option.solver_par.num_eigenvalues = 0;
	option.solver_par.restart = 100;
	option.solver_par.rtol = 1E-14;
	option.solver_par.solver = Magma_CGMERGE;
	option.solver_par.verbose = 0;
	option.solver_par.version = 0;

	magma_dsolverinfo_init(&option.solver_par, &option.precond_par, queue);
#endif
}

template<typename T> bool SparseMat<T>::is_empty() const { return triplet_mat.is_empty(); }

template<typename T> void SparseMat<T>::zeros() {
	triplet_mat.zeros();
	csc_mat.zeros();
	arma_mat.zeros();
}

template<typename T> void SparseMat<T>::reset() {
	triplet_mat.reset();
	csc_mat.reset();
	arma_mat.reset();
}

template<typename T> T SparseMat<T>::max() const { return triplet_mat.max(); }

template<typename T> const T& SparseMat<T>::operator()(const uword in_row, const uword in_col) const { return csc_mat(index_t(in_row), index_t(in_col)); }

template<typename T> T& SparseMat<T>::at(const uword in_row, const uword in_col) { return triplet_mat.at(index_t(in_row), index_t(in_col)); }

template<typename T> void SparseMat<T>::triplet_to_csc() { csc_mat = triplet_mat; }

template<typename T> MetaMat<T> SparseMat<T>::operator+(const MetaMat<T>& in_mat) {
	auto N = *this;
	N.triplet_mat += in_mat.triplet_mat;
	return N;
}

template<typename T> MetaMat<T> SparseMat<T>::operator-(const MetaMat<T>& in_mat) {
	auto N = *this;
	N.triplet_mat -= in_mat.triplet_mat;
	return N;
}

template<typename T> MetaMat<T>& SparseMat<T>::operator+=(const MetaMat<T>& in_mat) {
	triplet_mat += in_mat.triplet_mat;
	return dynamic_cast<MetaMat<T>&>(*this);
}

template<typename T> MetaMat<T>& SparseMat<T>::operator-=(const MetaMat<T>& in_mat) {
	triplet_mat -= in_mat.triplet_mat;
	return dynamic_cast<MetaMat<T>&>(*this);
}

template<typename T> MetaMat<T> SparseMat<T>::operator*(const T scalar) {
	auto N = *this;
	N.triplet_mat *= scalar;
	return N;
}

template<typename T> Mat<T> SparseMat<T>::operator*(const Mat<T>& in_mat) {
	Mat<T> out_mat(in_mat.n_rows, in_mat.n_cols, fill::zeros);

	if(in_mat.n_cols == 1) for(auto I = 0; I < triplet_mat.c_size; ++I) out_mat(triplet_mat.row_idx[I]) += triplet_mat.val_idx[I] * in_mat(triplet_mat.col_idx[I]);
	else for(auto I = 0; I < triplet_mat.c_size; ++I) for(auto J = 0; J < in_mat.n_cols; ++J) out_mat(triplet_mat.row_idx[I], J) += triplet_mat.val_idx[I] * in_mat(triplet_mat.col_idx[I], J);

	return out_mat;
}

template<typename T> MetaMat<T>& SparseMat<T>::operator*=(const T scalar) {
	triplet_mat *= scalar;
	return dynamic_cast<MetaMat<T>&>(*this);
}

template<typename T> Mat<T> SparseMat<T>::solve(const Mat<T>& in_mat) {
	Mat<T> out_mat;
	if(solve(out_mat, in_mat) != 0) out_mat.reset();
	return out_mat;
}

template<typename T> int SparseMat<T>::solve(Mat<T>& out_mat, const Mat<T>& in_mat) {
#ifdef SUANPAN_MAGMA
	out_mat = in_mat;

	csr_form<T> csr_mat(triplet_mat);

	magma_d_matrix A{Magma_CSR}, dA{Magma_CSR}, b{Magma_DENSE}, db{Magma_DENSE}, dx{Magma_DENSE};

	A.storage_type = Magma_CSR;
	A.memory_location = Magma_CPU;
	A.num_rows = magma_int_t(csr_mat.n_rows);
	A.num_cols = magma_int_t(csr_mat.n_cols);
	A.nnz = magma_int_t(csr_mat.c_size);
	A.true_nnz = magma_int_t(csr_mat.c_size);
	A.fill_mode = MagmaFull;
	A.sym = Magma_GENERAL;

	// save time
	A.val = csr_mat.val_idx;
	A.col = csr_mat.col_idx;
	A.row = csr_mat.row_ptr;

	b.ownership = MagmaFalse;
	b.memory_location = Magma_CPU;
	b.sym = Magma_GENERAL;
	b.diagorder_type = Magma_VALUE;
	b.fill_mode = MagmaFull;
	b.num_rows = magma_int_t(out_mat.n_rows);
	b.num_cols = 1;
	b.nnz = magma_int_t(out_mat.n_rows);
	b.max_nnz_row = 1;
	b.diameter = 0;
	b.blocksize = 1;
	b.numblocks = 1;
	b.alignment = 1;
	b.major = MagmaColMajor;
	b.ld = magma_int_t(out_mat.n_rows);
	b.val = out_mat.memptr();

	magma_d_precondsetup(A, b, &option.solver_par, &option.precond_par, queue);

	magma_dmtransfer(A, &dA, Magma_CPU, Magma_DEV, queue);
	magma_dmtransfer(b, &db, Magma_CPU, Magma_DEV, queue);
	magma_dmtransfer(b, &dx, Magma_CPU, Magma_DEV, queue);

	const auto info = magma_d_solver(dA, db, &dx, &option, queue);

	magma_dmtransfer(dx, &b, Magma_DEV, Magma_CPU, queue);

	memory::release(access::rw(out_mat.mem));
	access::rw(out_mat.mem) = b.val;

	printf("convergence = [\n");
	magma_dsolverinfo(&option.solver_par, &option.precond_par, queue);
	printf("];\n\n");

	magma_dmfree(&dA, queue);
	magma_dmfree(&dx, queue);
	magma_dmfree(&db, queue);

	return int(info);
#else
	triplet_to_csc();

	const uvec row_idx(csc_mat.row_idx, csc_mat.c_size, false, false);
	const uvec col_ptr(csc_mat.col_ptr, csc_mat.n_cols + 1, false, false);
	const Col<T> val_idx(csc_mat.val_idx, csc_mat.c_size, false, false);

	arma_mat = SpMat<T>(row_idx, col_ptr, val_idx, csc_mat.n_rows, csc_mat.n_cols);

	return spsolve(out_mat, arma_mat, in_mat, "s", symm_opts) ? 0 : -1;

	/*
	triplet_mat.csc_condense();

	out_mat = in_mat;

	DMUMPS_STRUC_C mumps_job;

	mumps_job.comm_fortran = -987654;
	mumps_job.par = 1;
	mumps_job.sym = 0;
	mumps_job.job = -1;
	dmumps_c(&mumps_job);

	mumps_job.n = int(triplet_mat.n_rows);
	mumps_job.nnz = int(triplet_mat.c_size);

	const auto l_irn = new int[mumps_job.nnz], l_jrn = new int[mumps_job.nnz];
	for(auto I = 0; I < mumps_job.nnz; ++I) {
		l_irn[I] = static_cast<int>(triplet_mat.row_idx[I] + 1);
		l_jrn[I] = static_cast<int>(triplet_mat.col_idx[I] + 1);
	}

	mumps_job.irn = l_irn;
	mumps_job.jcn = l_jrn;
	mumps_job.a = triplet_mat.val_idx;
	mumps_job.rhs = out_mat.memptr();

	mumps_job.icntl[0] = -1;
	mumps_job.icntl[1] = -1;
	mumps_job.icntl[2] = -1;
	mumps_job.icntl[3] = -1;

	mumps_job.job = 6;
	dmumps_c(&mumps_job);

	mumps_job.job = -2;
	dmumps_c(&mumps_job);

	delete[]l_irn;
	delete[]l_jrn;

	return mumps_job.info[0];
	*/
#endif
}

template<typename T> Mat<T> SparseMat<T>::solve_trs(const Mat<T>& in_mat) { return solve(in_mat); }

template<typename T> int SparseMat<T>::solve_trs(Mat<T>& out_mat, const Mat<T>& in_mat) { return solve(out_mat, in_mat); }

template<typename T> void SparseMat<T>::print() { arma_mat.print(); }

template<typename T> void SparseMat<T>::save(const char* name) { arma_mat.save(name, coord_ascii); }

#endif

//! @}
