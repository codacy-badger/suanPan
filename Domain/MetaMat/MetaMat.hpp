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
 * @class MetaMat
 * @brief A MetaMat class that holds matrices.
 *
 * @author tlc
 * @date 08/09/2017
 * @version 0.2.0
 * @file MetaMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef METAMAT_HPP
#define METAMAT_HPP

#include <Toolbox/debug.h>
#include <Domain/MetaMat/csc_form.hpp>
#include <Domain/MetaMat/csr_form.hpp>
#include <Domain/MetaMat/triplet_form.hpp>

template<typename T> class MetaMat {
protected:
	static const char TRAN;
public:
	triplet_form<T> triplet_mat;
	csc_form<T> csc_mat;
	SpMat<T> arma_mat;

	Col<int> IPIV;
	bool factored = false;
	const uword n_rows;
	const uword n_cols;
	const uword n_elem;

	const T* const memory = nullptr;

	MetaMat();
	MetaMat(uword, uword, uword);
	MetaMat(const MetaMat&);
	MetaMat(MetaMat&&) noexcept;
	MetaMat& operator=(const MetaMat&);
	MetaMat& operator=(MetaMat&&) noexcept;
	virtual ~MetaMat();

	virtual bool is_empty() const;
	virtual void init();
	virtual void zeros();
	virtual void reset();

	virtual void unify(uword);

	virtual T max() const;

	virtual const T& operator()(uword, uword) const;
	virtual T& at(uword, uword);

	virtual const T* memptr() const;
	virtual T* memptr();

	virtual MetaMat operator+(const MetaMat&);
	virtual MetaMat operator-(const MetaMat&);
	virtual MetaMat& operator+=(const MetaMat&);
	virtual MetaMat& operator-=(const MetaMat&);

	virtual MetaMat operator*(T);
	virtual Mat<T> operator*(const Mat<T>&);

	virtual MetaMat& operator*=(T);

	virtual Mat<T> solve(const Mat<T>&);
	virtual int solve(Mat<T>&, const Mat<T>&);

	virtual Mat<T> solve_trs(const Mat<T>&);
	virtual int solve_trs(Mat<T>&, const Mat<T>&);

	virtual MetaMat factorize();

	virtual MetaMat i();
	virtual MetaMat inv();

	virtual void print();
	virtual void save(const char*);
};

template<typename T> const char MetaMat<T>::TRAN = 'N';

template<typename T> Mat<T> to_mat(const MetaMat<T>& in_mat) {
	Mat<T> out_mat(in_mat.n_rows, in_mat.n_cols, fill::zeros);

	for(uword I = 0; I < in_mat.n_rows; ++I) for(uword J = 0; J < in_mat.n_cols; ++J) out_mat(I, J) = in_mat(I, J);

	return out_mat;
}

template<typename T> MetaMat<T>::MetaMat()
	: n_rows(0)
	, n_cols(0)
	, n_elem(0) {}

template<typename T> MetaMat<T>::MetaMat(const uword in_rows, const uword in_cols, const uword in_elem)
	: triplet_mat(index_t(in_rows), index_t(in_cols))
	, n_rows(in_rows)
	, n_cols(in_cols)
	, n_elem(in_elem) {
	MetaMat<T>::init();
	MetaMat<T>::zeros();
}

template<typename T> MetaMat<T>::MetaMat(const MetaMat& old_mat)
	: triplet_mat(old_mat.triplet_mat)
	, factored(old_mat.factored)
	, n_rows(old_mat.n_rows)
	, n_cols(old_mat.n_cols)
	, n_elem(old_mat.n_elem) {
	MetaMat<T>::init();
	if(old_mat.memptr() != nullptr) memcpy(MetaMat<T>::memptr(), old_mat.memptr(), old_mat.n_elem * sizeof(T));
}

template<typename T> MetaMat<T>::MetaMat(MetaMat&& old_mat) noexcept
	: triplet_mat(std::forward<triplet_form<T>>(old_mat.triplet_mat))
	, factored(old_mat.factored)
	, n_rows(old_mat.n_rows)
	, n_cols(old_mat.n_cols)
	, n_elem(old_mat.n_elem) {
	access::rw(memory) = old_mat.memory;
	access::rw(old_mat.memory) = nullptr;
}

template<typename T> MetaMat<T>& MetaMat<T>::operator=(const MetaMat& old_mat) {
	if(this != &old_mat) {
		factored = old_mat.factored;
		triplet_mat = old_mat.triplet_mat;
		access::rw(n_rows) = old_mat.n_rows;
		access::rw(n_cols) = old_mat.n_cols;
		access::rw(n_elem) = old_mat.n_elem;
		init();
		if(old_mat.memptr() != nullptr) memcpy(memptr(), old_mat.memptr(), old_mat.n_elem * sizeof(T));
	}
	return *this;
}

template<typename T> MetaMat<T>& MetaMat<T>::operator=(MetaMat&& old_mat)
noexcept {
	if(this != &old_mat) {
		factored = old_mat.factored;
		triplet_mat = old_mat.triplet_mat;
		access::rw(n_rows) = old_mat.n_rows;
		access::rw(n_cols) = old_mat.n_cols;
		access::rw(n_elem) = old_mat.n_elem;
		access::rw(memory) = old_mat.memory;
		access::rw(old_mat.memory) = nullptr;
	}
	return *this;
}

template<typename T> MetaMat<T>::~MetaMat() { if(memory != nullptr) memory::release(access::rw(memory)); }

template<typename T> bool MetaMat<T>::is_empty() const { return 0 == n_elem; }

template<typename T> void MetaMat<T>::init() {
	if(memory != nullptr) memory::release(access::rw(memory));
	access::rw(memory) = is_empty() ? nullptr : memory::acquire<T>(n_elem);
}

template<typename T> void MetaMat<T>::zeros() {
	arrayops::fill_zeros(memptr(), n_elem);
	factored = false;
}

template<typename T> void MetaMat<T>::reset() {
	access::rw(n_rows) = 0;
	access::rw(n_cols) = 0;
	access::rw(n_elem) = 0;
	if(memory != nullptr) memory::release(access::rw(memory));
	access::rw(memory) = nullptr;
	factored = false;
}

template<typename T> void MetaMat<T>::unify(const uword idx) {
	for(uword I = 0; I < n_rows; ++I) at(I, idx) = 0.;
	for(uword I = 0; I < n_cols; ++I) at(idx, I) = 0.;
	at(idx, idx) = 1.;
}

template<typename T> T MetaMat<T>::max() const { return op_max::direct_max(memptr(), n_elem); }

template<typename T> const T& MetaMat<T>::operator()(const uword in_row, const uword in_col) const { return memory[in_row + in_col * n_rows]; }

template<typename T> T& MetaMat<T>::at(const uword in_row, const uword in_col) { return access::rw(memory[in_row + in_col * n_rows]); }

template<typename T> const T* MetaMat<T>::memptr() const { return memory; }

template<typename T> T* MetaMat<T>::memptr() { return const_cast<T*>(memory); }

template<typename T> MetaMat<T> MetaMat<T>::operator+(const MetaMat& M) {
	auto N = *this;
	N += M;
	return N;
}

template<typename T> MetaMat<T> MetaMat<T>::operator-(const MetaMat& M) {
	auto N = *this;
	N -= M;
	return N;
}

template<typename T> MetaMat<T>& MetaMat<T>::operator+=(const MetaMat& M) {
	if(n_rows == M.n_rows && n_cols == M.n_cols && n_elem == M.n_elem) {
		arrayops::inplace_plus(memptr(), M.memptr(), n_elem);
		factored = false;
	}
	return *this;
}

template<typename T> MetaMat<T>& MetaMat<T>::operator-=(const MetaMat& M) {
	if(n_rows == M.n_rows && n_cols == M.n_cols && n_elem == M.n_elem) {
		arrayops::inplace_minus(memptr(), M.memptr(), n_elem);
		factored = false;
	}
	return *this;
}

template<typename T> MetaMat<T> MetaMat<T>::operator*(const T value) {
	auto new_mat = *this;
	arrayops::inplace_mul(new_mat.memptr(), value, new_mat.n_elem);
	return new_mat;
}

template<typename T> Mat<T> MetaMat<T>::operator*(const Mat<T>& B) {
	auto C = B;

	if(1 == B.n_cols) {
		auto M = static_cast<int>(n_rows);
		auto N = static_cast<int>(n_cols);
		T ALPHA = 1.;
		auto LDA = M;
		auto INCX = 1;
		T BETA = 0.;
		auto INCY = 1;

		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_sgemv)(&TRAN, &M, &N, (E*)(&ALPHA), (E*)(memptr()), &LDA, (E*)(B.memptr()), &INCX, (E*)(&BETA), (E*)(C.memptr()), &INCY);
		} else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgemv)(&TRAN, &M, &N, (E*)(&ALPHA), (E*)(memptr()), &LDA, (E*)(B.memptr()), &INCX, (E*)(&BETA), (E*)(C.memptr()), &INCY);
		}
	} else {
		auto M = static_cast<int>(n_rows);
		auto N = static_cast<int>(B.n_cols);
		auto K = static_cast<int>(n_cols);
		T ALPHA = 1.;
		auto LDA = M;
		auto LDB = K;
		T BETA = 0.;
		auto LDC = M;

		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_sgemm)(&TRAN, &TRAN, &M, &N, &K, (E*)(&ALPHA), (E*)(memptr()), &LDA, (E*)(B.memptr()), &LDB, (E*)(&BETA), (E*)(C.memptr()), &LDC);
		} else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgemm)(&TRAN, &TRAN, &M, &N, &K, (E*)(&ALPHA), (E*)(memptr()), &LDA, (E*)(B.memptr()), &LDB, (E*)(&BETA), (E*)(C.memptr()), &LDC);
		}
	}

	return C;
}

template<typename T> MetaMat<T>& MetaMat<T>::operator*=(const T value) {
	arrayops::inplace_mul(memptr(), value, n_elem);
	return *this;
}

template<typename T> Mat<T> MetaMat<T>::solve(const Mat<T>& B) {
	Mat<T> X;
	if(solve(X, B) != 0) X.reset();
	return X;
}

template<typename T> int MetaMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
	if(factored) {
		suanpan_debug("the matrix is factored.\n");
		return solve_trs(X, B);
	}

	auto INFO = 0;

#ifdef SUANPAN_MAGMA
	X.set_size(B.n_rows, B.n_cols);

	const auto N = magma_int_t(n_rows);
	const auto NRHS = magma_int_t(B.n_cols);
	const auto LDA = N;
	const auto LDB = magma_int_t(B.n_rows);
	IPIV.zeros(N);

	magma_queue_t QUEUE;
	magma_queue_create(0, &QUEUE);
	const auto LDDA = magma_roundup(N, 64);

	if(std::is_same<T, float>::value) {
		using E = float;
		magmaFloat_ptr GPUA, GPUB;
		magma_smalloc(&GPUA, LDDA * N);
		magma_smalloc(&GPUB, LDDA * NRHS);
		magma_ssetmatrix(N, N, (E*)(this->memptr()), LDA, GPUA, LDDA, QUEUE);
		magma_ssetmatrix(N, NRHS, (E*)(B.memptr()), LDB, GPUB, LDDA, QUEUE);
		magma_sgesv_gpu(N, NRHS, GPUA, LDDA, IPIV.memptr(), GPUB, LDDA, &INFO);
		magma_sgetmatrix(N, NRHS, GPUB, LDDA, (E*)(X.memptr()), LDB, QUEUE);
		magma_free(GPUA);
		magma_free(GPUB);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		magmaDouble_ptr GPUA, GPUB;
		magma_dmalloc(&GPUA, LDDA * N);
		magma_dmalloc(&GPUB, LDDA * NRHS);
		magma_dsetmatrix(N, N, (E*)(this->memptr()), LDA, GPUA, LDDA, QUEUE);
		magma_dsetmatrix(N, NRHS, (E*)(B.memptr()), LDB, GPUB, LDDA, QUEUE);
		magma_dgesv_gpu(N, NRHS, GPUA, LDDA, IPIV.memptr(), GPUB, LDDA, &INFO);
		magma_dgetmatrix(N, NRHS, GPUB, LDDA, (E*)(X.memptr()), LDB, QUEUE);
		magma_free(GPUA);
		magma_free(GPUB);
	}
#else
	X = B;

	auto N = int(n_rows);
	auto NRHS = int(B.n_cols);
	auto LDA = N;
	auto LDB = int(B.n_rows);
	IPIV.zeros(N);

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgesv)(&N, &NRHS, (E*)(this->memptr()), &LDA, IPIV.memptr(), (E*)(X.memptr()), &LDB, &INFO);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgesv)(&N, &NRHS, (E*)(this->memptr()), &LDA, IPIV.memptr(), (E*)(X.memptr()), &LDB, &INFO);
	}

	if(INFO == 0) factored = true;
#endif

	if(INFO != 0) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;

}

template<typename T> Mat<T> MetaMat<T>::solve_trs(const Mat<T>& B) {
	Mat<T> X;
	if(solve_trs(X, B) != 0) X.reset();
	return X;
}

template<typename T> int MetaMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	if(!factored) {
		suanpan_debug("the matrix is not factored.\n");
		return solve(X, B);
	}

	if(IPIV.is_empty()) return -1;

	X = B;

	auto N = static_cast<int>(n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDA = N;
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgetrs)(const_cast<char*>(&TRAN), &N, &NRHS, (E*)(this->memptr()), &LDA, IPIV.memptr(), (E*)(X.memptr()), &LDB, &INFO);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgetrs)(const_cast<char*>(&TRAN), &N, &NRHS, (E*)(this->memptr()), &LDA, IPIV.memptr(), (E*)(X.memptr()), &LDB, &INFO);
	}

	return INFO;
}

template<typename T> MetaMat<T> MetaMat<T>::factorize() {
	auto X = *this;

	if(factored) {
		suanpan_warning("the matrix is factored.\n");
		return X;
	}

	arma_debug_check(X.n_rows != X.n_cols, "i() only accepts sqaure matrix.");

	auto M = static_cast<int>(X.n_rows);
	auto N = M;
	auto LDA = M;
	X.IPIV.zeros(N);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgetrf)(&M, &N, (E*)(X.memptr()), &LDA, X.IPIV.memptr(), &INFO);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgetrf)(&M, &N, (E*)(X.memptr()), &LDA, X.IPIV.memptr(), &INFO);
	}

	if(INFO != 0) {
		suanpan_error("factorize() fails.\n");
		X.reset();
	} else X.factored = true;

	return X;
}

template<typename T> MetaMat<T> MetaMat<T>::i() {
	auto X = *this;

	arma_debug_check(X.n_rows != X.n_cols, "i() only accepts sqaure matrix.");

	auto M = static_cast<int>(X.n_rows);
	auto N = M;
	auto LDA = M;
	X.IPIV.zeros(N);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgetrf)(&M, &N, (E*)(X.memptr()), &LDA, X.IPIV.memptr(), &INFO);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgetrf)(&M, &N, (E*)(X.memptr()), &LDA, X.IPIV.memptr(), &INFO);
	}

	if(INFO != 0) {
		X.reset();
		return X;
	}

	auto LWORK = 8 * M;
	const auto WORK = new T[LWORK];

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgetri)(&N, (E*)(X.memptr()), &LDA, X.IPIV.memptr(), (E*)(WORK), &LWORK, &INFO);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgetri)(&N, (E*)(X.memptr()), &LDA, X.IPIV.memptr(), (E*)(WORK), &LWORK, &INFO);
	}

	delete[] WORK;

	if(INFO != 0) X.reset();

	return X;
}

template<typename T> MetaMat<T> MetaMat<T>::inv() { return i(); }

template<typename T> void MetaMat<T>::print() {
	for(uword I = 0; I < n_rows; ++I) {
		for(uword J = 0; J < n_cols; ++J) suanpan_info("%+.3E\t", operator()(I, J));
		suanpan_info("\n");
	}
	suanpan_info("\n");
}

template<typename T> void MetaMat<T>::save(const char* name) { if(!to_mat(*this).save(name, raw_ascii)) suanpan_error("cannot save matrix to file.\n"); }

#endif

//! @}