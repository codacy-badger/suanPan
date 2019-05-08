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
 * @class BandMat
 * @brief A BandMat class that holds matrices.
 *
 * @author tlc
 * @date 06/09/2017
 * @version 0.1.0
 * @file BandMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef BANDMAT_HPP
#define BANDMAT_HPP

#include <Toolbox/debug.h>

template<typename T> class BandMat final : public MetaMat<T> {
	using MetaMat<T>::TRAN;

	static T bin;
	using MetaMat<T>::i;
	using MetaMat<T>::inv;

	const uword l_band;
	const uword u_band;
	const uword s_band;
	const uword m_rows; // memory block layout
public:
	using MetaMat<T>::IPIV;
	using MetaMat<T>::factored;
	using MetaMat<T>::n_cols;
	using MetaMat<T>::n_rows;
	using MetaMat<T>::n_elem;
	using MetaMat<T>::memory;
	using MetaMat<T>::solve;
	using MetaMat<T>::solve_trs;
	using MetaMat<T>::factorize;

	BandMat();
	BandMat(uword, uword, uword);
	BandMat(const BandMat&);
	~BandMat() = default;

	const T& operator()(uword, uword) const override;
	T& at(uword, uword) override;

	Mat<T> operator*(const Mat<T>&) override;

	int solve(Mat<T>&, const Mat<T>&) override;
	int solve_trs(Mat<T>&, const Mat<T>&) override;

	MetaMat<T> factorize() override;
};

template<typename T> T BandMat<T>::bin = 0.;

template<typename T> BandMat<T>::BandMat()
	: MetaMat<T>()
	, l_band(0)
	, u_band(0)
	, s_band(0)
	, m_rows(0) {}

template<typename T> BandMat<T>::BandMat(const uword in_size, const uword in_l, const uword in_u)
	: MetaMat<T>(in_size, in_size, (2 * in_l + in_u + 1) * in_size)
	, l_band(in_l)
	, u_band(in_u)
	, s_band(l_band + u_band)
	, m_rows(2 * in_l + in_u + 1) {}

template<typename T> BandMat<T>::BandMat(const BandMat& old_mat)
	: MetaMat<T>(old_mat)
	, l_band(old_mat.l_band)
	, u_band(old_mat.u_band)
	, s_band(old_mat.s_band)
	, m_rows(old_mat.m_rows) {}

template<typename T> const T& BandMat<T>::operator()(const uword in_row, const uword in_col) const {
	const auto n_bw = int(in_row - in_col);
	if(n_bw > int(l_band) || n_bw < -int(u_band)) return bin = 0.;

	return memory[in_row - in_col + s_band + in_col * m_rows];
}

template<typename T> T& BandMat<T>::at(const uword in_row, const uword in_col) {
	const auto n_bw = int(in_row) - int(in_col);
	if(n_bw > int(l_band) || n_bw < -int(u_band)) return bin = 0.;

	return access::rw(memory[in_row - in_col + s_band + in_col * m_rows]);
}

template<typename T> Mat<T> BandMat<T>::operator*(const Mat<T>& X) {
	if(!X.is_colvec()) throw invalid_argument("requires a coolumn vector");

	auto Y = X;

	auto M = static_cast<int>(n_rows);
	auto N = static_cast<int>(n_cols);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	T ALPHA = 1.;
	auto LDA = static_cast<int>(m_rows);
	auto INC = 1;
	T BETA = 0.;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)(&ALPHA), (E*)(this->memptr() + l_band), &LDA, (E*)(X.memptr()), &INC, (E*)(&BETA), (E*)(Y.memptr()), &INC);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)(&ALPHA), (E*)(this->memptr() + l_band), &LDA, (E*)(X.memptr()), &INC, (E*)(&BETA), (E*)(Y.memptr()), &INC);
	}

	return Y;
}

template<typename T> int BandMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
	if(factored) {
		suanpan_warning("the matrix is factored.\n");
		return this->solve_trs(X, B);
	}

	suanpan_debug([&]() { if(n_rows != n_cols) throw invalid_argument("requires a square matrix"); });

	X = B;

	auto N = static_cast<int>(n_rows);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);
	IPIV.zeros(N);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgbsv)(&N, &KL, &KU, &NRHS, (E*)(this->memptr()), &LDAB, IPIV.memptr(), (E*)(X.memptr()), &LDB, &INFO);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgbsv)(&N, &KL, &KU, &NRHS, (E*)(this->memptr()), &LDAB, IPIV.memptr(), (E*)(X.memptr()), &LDB, &INFO);
	}

	if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);
	else factored = true;

	return INFO;
}

template<typename T> int BandMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	if(!factored) {
		suanpan_warning("the matrix is not factored.\n");
		return this->solve(X, B);
	}

	if(IPIV.is_empty()) return -1;

	suanpan_debug([&]() { if(n_rows != n_cols) throw invalid_argument("requires a square matrix");; });

	X = B;

	auto N = static_cast<int>(n_rows);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)(this->memptr()), &LDAB, IPIV.memptr(), (E*)(X.memptr()), &LDB, &INFO);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)(this->memptr()), &LDAB, IPIV.memptr(), (E*)(X.memptr()), &LDB, &INFO);
	}

	if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> MetaMat<T> BandMat<T>::factorize() {
	auto X = *this;

	if(factored) {
		suanpan_warning("the matrix is factored.\n");
		return X;
	}

	suanpan_debug([&]() { if(n_rows != n_cols) throw invalid_argument("requires a square matrix"); });

	auto M = static_cast<int>(n_rows);
	auto N = static_cast<int>(n_cols);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	auto LDAB = static_cast<int>(m_rows);
	X.IPIV.zeros(N);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgbtrf)(&M, &N, &KL, &KU, (E*)(X.memptr()), &LDAB, X.IPIV.memptr(), &INFO);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgbtrf)(&M, &N, &KL, &KU, (E*)(X.memptr()), &LDAB, X.IPIV.memptr(), &INFO);
	}

	if(INFO != 0) {
		suanpan_error("factorize() fails.\n");
		X.reset();
	} else X.factored = true;

	return X;
}

#endif

//! @}