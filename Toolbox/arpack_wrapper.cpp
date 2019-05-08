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

#include "arpack_wrapper.h"

/*
int eig_solve(cx_vec& eigval, cx_mat& eigvec, mat& K, const unsigned num, const char* form) {
	auto IDO = 0;
	auto BMAT = 'I'; // standard eigenvalue problem A*x=lambda*x
	auto N = static_cast<int>(K.n_rows);
	char WHICH[2];
	for(auto I = 0; I < 2; ++I) WHICH[I] = form[I];
	auto NEV = std::min(static_cast<int>(num), N - 2);
	auto TOL = 0.;
	auto NCV = std::min(std::min(NEV + 8, 2 * NEV + 1), N) + 2;
	auto LDV = N;
	auto LWORKL = 3 * NCV * (NCV + 2);
	auto INFO = 0;

	podarray<int> IPARAM(11), IPNTR(14);
	podarray<double> RESID(N), V(N * NCV), WORKD(3 * N), WORKL(LWORKL);

	IPARAM(0) = 1;   // exact shift
	IPARAM(2) = 500; // maximum iteration
	IPARAM(6) = 1;   // mode 1: A*x=lambda*x

	while(IDO != 99) {
		arma_fortran(arma_dnaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
		if(IDO == 1 || IDO == -1) {
			const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
			vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false);
			Y = K * X;
		}
	}

	if(INFO != 0) return INFO;

	auto RVEC = 1;
	auto HOWMNY = 'A';
	auto LDZ = N;
	auto SIGMAR = 0.;
	auto SIGMAI = 0.;

	podarray<int> SELECT(NCV);
	podarray<double> DR(NEV + 1), DI(NEV + 1), Z(N * (NEV + 1)), WORKEV(3 * NCV);

	arma_fortran(arma_dneupd)(&RVEC, &HOWMNY, SELECT.memptr(), DR.memptr(), DI.memptr(), Z.memptr(), &LDZ, &SIGMAR, &SIGMAI, WORKEV.memptr(), &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

	eigval.set_size(NEV);
	eigvec.set_size(N, NEV);
	using T = double;
	for(auto I = 0; I < NEV; ++I) eigval(I) = std::complex<T>(DR(I), DI(I));
	for(auto I = 0; I < NEV; ++I) {
		if(I < NEV - 1 && eigval[I] == std::conj(eigval[I + 1])) {
			for(auto J = 0; J < N; ++J) {
				eigvec.at(J, I) = std::complex<T>(Z[N * I + J], Z[N * (I + 1) + J]);
				eigvec.at(J, I + 1) = std::complex<T>(Z[N * I + J], -Z[N * (I + 1) + J]);
			}
			++I;
		} else if(I == NEV - 1 && std::complex<T>(eigval[I]).imag() != 0.) for(auto J = 0; J < N; ++J) eigvec.at(J, I) = std::complex<T>(Z[N * I + J], Z[N * (I + 1) + J]);
		else for(auto J = 0; J < N; ++J) eigvec.at(J, I) = std::complex<T>(Z[N * I + J], 0.);
	}

	return INFO;
}

int eig_solve(cx_vec& eigval, cx_mat& eigvec, mat& K, mat& M, const unsigned num, const char* form) {
	auto IDO = 0;
	auto BMAT = 'G'; // generalized eigenvalue problem A*x=lambda*B*x
	auto N = static_cast<int>(K.n_rows);
	char WHICH[2];
	for(auto I = 0; I < 2; ++I) WHICH[I] = form[I];
	auto NEV = std::min(static_cast<int>(num), N - 2);
	auto TOL = 0.;
	auto NCV = std::min(std::min(NEV + 8, 2 * NEV + 1), N) + 2;
	auto LDV = N;
	auto LWORKL = 3 * NCV * (NCV + 2);
	auto INFO = 0;

	podarray<int> IPARAM(11), IPNTR(14);
	podarray<double> RESID(N), V(N * NCV), WORKD(3 * N), WORKL(LWORKL);

	IPARAM(0) = 1;   // exact shift
	IPARAM(2) = 500; // maximum iteration
	IPARAM(6) = 2;   // mode 2: A*x=lambda*M*x

	while(IDO != 99) {
		arma_fortran(arma_dnaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
		vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
		vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false);
		if(IDO == 1 || IDO == -1) {
			X = K * X;
			Y = solve(M, X);
		} else if(IDO == 2) Y = M * X;
	}

	if(INFO != 0) return INFO;

	auto RVEC = 1;
	auto HOWMNY = 'A';
	auto LDZ = N;
	auto SIGMAR = 0.;
	auto SIGMAI = 0.;

	podarray<int> SELECT(NCV);
	podarray<double> DR(NEV + 1), DI(NEV + 1), Z(N * (NEV + 1)), WORKEV(3 * NCV);

	arma_fortran(arma_dneupd)(&RVEC, &HOWMNY, SELECT.memptr(), DR.memptr(), DI.memptr(), Z.memptr(), &LDZ, &SIGMAR, &SIGMAI, WORKEV.memptr(), &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

	eigval.set_size(NEV);
	eigvec.set_size(N, NEV);
	using T = double;
	for(auto I = 0; I < NEV; ++I) eigval(I) = std::complex<T>(DR(I), DI(I));
	for(auto I = 0; I < NEV; ++I) {
		if(I < NEV - 1 && eigval[I] == std::conj(eigval[I + 1])) {
			for(auto J = 0; J < N; ++J) {
				eigvec.at(J, I) = std::complex<T>(Z[N * I + J], Z[N * (I + 1) + J]);
				eigvec.at(J, I + 1) = std::complex<T>(Z[N * I + J], -Z[N * (I + 1) + J]);
			}
			++I;
		} else if(I == NEV - 1 && std::complex<T>(eigval[I]).imag() != 0.) for(auto J = 0; J < N; ++J) eigvec.at(J, I) = std::complex<T>(Z[N * I + J], Z[N * (I + 1) + J]);
		else for(auto J = 0; J < N; ++J) eigvec.at(J, I) = std::complex<T>(Z[N * I + J], 0.);
	}

	return INFO;
}

int eig_solve(vec& eigval, mat& eigvec, mat& K, const unsigned num, const char* form) {
	auto IDO = 0;
	auto BMAT = 'I'; // standard eigenvalue problem A*x=lambda*x
	auto N = static_cast<int>(K.n_rows);
	char WHICH[2];
	for(auto I = 0; I < 2; ++I) WHICH[I] = form[I];
	auto NEV = std::min(static_cast<int>(num), N - 1);
	auto TOL = 1E-10;
	auto NCV = std::min(std::max(NEV + 3, 2 * NEV + 1), N) + 2;
	auto LDV = N;
	auto LWORKL = 3 * NCV * (NCV + 2);
	auto INFO = 0;

	podarray<int> IPARAM(11), IPNTR(14);
	podarray<double> RESID(N), V(N * NCV), WORKD(3 * N), WORKL(LWORKL);

	IPARAM.zeros();
	IPARAM(0) = 1;    // exact shift
	IPARAM(2) = 1000; // maximum iteration
	IPARAM(6) = 1;    // mode 1: A*x=lambda*x

	while(IDO != 99) {
		arma_fortran(arma_dsaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
		if(IDO == 1 || IDO == -1) {
			const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
			vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false);
			Y = K * X;
			Y.print();
			(K * X).print();
		}
	}

	if(INFO != 0) return INFO;

	auto RVEC = 1;
	auto HOWMNY = 'A';
	auto LDZ = N;
	auto SIGMA = 0.;

	podarray<int> SELECT(NCV);

	eigval.set_size(NEV);
	eigvec.set_size(N, NEV);

	arma_fortran(arma_dseupd)(&RVEC, &HOWMNY, SELECT.memptr(), eigval.memptr(), eigvec.memptr(), &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

	return INFO;
}

int eig_solve(vec& eigval, mat& eigvec, mat& K, mat& M, const unsigned num, const char* form) {
	auto IDO = 0;
	auto BMAT = 'G'; // generalized eigenvalue problem A*x=lambda*M*x
	auto N = static_cast<int>(K.n_rows);
	char WHICH[2];
	for(auto I = 0; I < 2; ++I) WHICH[I] = form[I];
	auto NEV = std::min(static_cast<int>(num), N - 1);
	auto TOL = 0.;
	auto NCV = std::min(std::min(NEV + 8, 2 * NEV + 1), N) + 2;
	auto LDV = N;
	auto LWORKL = NCV * (NCV + 8);
	auto INFO = 0;

	podarray<int> IPARAM(11), IPNTR(14);
	podarray<double> RESID(N), V(N * NCV), WORKD(3 * N), WORKL(LWORKL);

	IPARAM(0) = 1;   // exact shift
	IPARAM(2) = 500; // maximum iteration
	IPARAM(6) = 2;   // mode 2: A*x=lambda*M*x

	while(IDO != 99) {
		arma_fortran(arma_dsaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
		vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
		vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false);
		if(IDO == 1 || IDO == -1) {
			X = K * X;
			Y = solve(M, X);
		} else if(IDO == 2) Y = M * X;
	}

	if(INFO != 0) return INFO;

	auto RVEC = 1;
	auto HOWMNY = 'A';
	auto LDZ = N;
	auto SIGMA = 0.;

	podarray<int> SELECT(NCV);

	eigval.set_size(NEV);
	eigvec.set_size(N, NEV);

	arma_fortran(arma_dseupd)(&RVEC, &HOWMNY, SELECT.memptr(), eigval.memptr(), eigvec.memptr(), &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

	return INFO;
}

int eig_solve(vec& eigval, mat& eigvec, const std::shared_ptr<MetaMat<double>>& K, const unsigned num, const char* form) {
	auto IDO = 0;
	auto BMAT = 'I'; // standard eigenvalue problem A*x=lambda*x
	auto N = static_cast<int>(K->n_cols);
	char WHICH[2];
	for(auto I = 0; I < 2; ++I) WHICH[I] = form[I];
	auto NEV = std::min(static_cast<int>(num), N - 1);
	auto TOL = 0.;
	auto NCV = std::min(std::min(NEV + 8, 2 * NEV + 1), N) + 2;
	auto LDV = N;
	auto LWORKL = NCV * (NCV + 8);
	auto INFO = 0;

	podarray<int> IPARAM(11), IPNTR(14);
	podarray<double> RESID(N), V(N * NCV), WORKD(3 * N), WORKL(LWORKL);

	IPARAM(0) = 1;   // exact shift
	IPARAM(2) = 500; // maximum iteration
	IPARAM(6) = 1;   // mode 1: A*x=lambda*x

	while(IDO != 99) {
		arma_fortran(arma_dsaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
		if(IDO == 1 || IDO == -1) {
			const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
			vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false);
			Y = *K * X;
		}
	}

	if(INFO != 0) return INFO;

	auto RVEC = 1;
	auto HOWMNY = 'A';
	auto LDZ = N;
	auto SIGMA = 0.;

	podarray<int> SELECT(NCV);

	eigval.set_size(NEV);
	eigvec.set_size(N, NEV);

	arma_fortran(arma_dseupd)(&RVEC, &HOWMNY, SELECT.memptr(), eigval.memptr(), eigvec.memptr(), &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

	return INFO;
}

int eig_solve(vec& eigval, mat& eigvec, const std::shared_ptr<MetaMat<double>>& K, const std::shared_ptr<MetaMat<double>>& M, const unsigned num, const char* form) {
	auto IDO = 0;
	auto BMAT = 'G'; // generalized eigenvalue problem A*x=lambda*M*x
	auto N = static_cast<int>(K->n_cols);
	char WHICH[2];
	for(auto I = 0; I < 2; ++I) WHICH[I] = form[I];
	auto NEV = std::min(static_cast<int>(num), N - 1);
	auto TOL = 0.;
	auto NCV = std::min(std::min(NEV + 8, 2 * NEV + 1), N) + 2;
	auto LDV = N;
	auto LWORKL = NCV * (NCV + 8);
	auto INFO = 0;

	podarray<int> IPARAM(11), IPNTR(14);
	podarray<double> RESID(N), V(N * NCV), WORKD(3 * N), WORKL(LWORKL);

	IPARAM(0) = 1;   // exact shift
	IPARAM(2) = 500; // maximum iteration
	IPARAM(6) = 2;   // mode 2: A*x=lambda*M*x

	auto MCOPY = *M;
	auto counter = 0;
	while(IDO != 99) {
		arma_fortran(arma_dsaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
		vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
		vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false);
		if(IDO == 1 || IDO == -1) {
			X = *K * X;
			Y = counter++ == 0 ? MCOPY.solve(X) : MCOPY.solve_trs(X);
		} else if(IDO == 2) Y = *M * X;
	}

	if(INFO != 0) return INFO;

	auto RVEC = 1;
	auto HOWMNY = 'A';
	auto LDZ = N;
	auto SIGMA = 0.;

	podarray<int> SELECT(NCV);

	eigval.set_size(NEV);
	eigvec.set_size(N, NEV);

	arma_fortran(arma_dseupd)(&RVEC, &HOWMNY, SELECT.memptr(), eigval.memptr(), eigvec.memptr(), &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

	return INFO;
}

int eig_solve2(vec& eigval, mat& eigvec, const std::shared_ptr<MetaMat<double>>& K, const std::shared_ptr<MetaMat<double>>& M, const unsigned num, const char* form) {
	auto IDO = 0;
	auto BMAT = 'G'; // generalized eigenvalue problem A*x=lambda*M*x
	auto N = static_cast<int>(K->n_cols);
	char WHICH[2];
	for(auto I = 0; I < 2; ++I) WHICH[I] = form[I];
	auto NEV = std::min(static_cast<int>(num), N - 1);
	auto TOL = 0.;
	auto NCV = std::min(std::min(NEV + 8, 2 * NEV + 1), N) + 2;
	auto LDV = N;
	auto LWORKL = NCV * (NCV + 8);
	auto INFO = 0;

	podarray<int> IPARAM(11), IPNTR(14);
	podarray<double> RESID(N), V(N * NCV), WORKD(3 * N), WORKL(LWORKL);

	IPARAM(0) = 1;   // exact shift
	IPARAM(2) = 500; // maximum iteration
	IPARAM(6) = 3;   // mode 3: K*x=lambda*M*x

	auto KCOPY = *K;
	while(IDO != 99) {
		arma_fortran(arma_dsaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
		vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false);
		if(IDO == -1) {
			vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
			X = *M * X;
			Y = KCOPY.solve(X);
		} else if(IDO == 1) {
			const vec X(WORKD.memptr() + IPNTR[2] - 1, N, false);
			Y = KCOPY.solve_trs(X);
		} else if(IDO == 2) {
			const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
			Y = *M * X;
		}
	}

	if(INFO != 0) return INFO;

	auto RVEC = 1;
	auto HOWMNY = 'A';
	auto LDZ = N;
	auto SIGMA = 0.;

	podarray<int> SELECT(NCV);

	eigval.set_size(NEV);
	eigvec.set_size(N, NEV);

	arma_fortran(arma_dseupd)(&RVEC, &HOWMNY, SELECT.memptr(), eigval.memptr(), eigvec.memptr(), &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

	return INFO;
}
*/

int eig_solve(vec& eigval, mat& eigvec, const std::shared_ptr<MetaMat<double>>& K, const std::shared_ptr<MetaMat<double>>& KG, const unsigned num, const char* form) {
	auto IDO = 0;
	auto BMAT = 'G'; // generalized eigenvalue problem A*x=lambda*M*x
	auto N = static_cast<int>(K->n_cols);
	char WHICH[2];
	for(auto I = 0; I < 2; ++I) WHICH[I] = form[I];
	auto NEV = std::min(static_cast<int>(num), N - 1);
	auto TOL = 1E-10;
	auto NCV = std::min(std::max(NEV + 2, 2 * NEV) + 5, N);
	auto LDV = N;
	auto LWORKL = 2 * NCV * (NCV + 8);
	auto INFO = 0;

	podarray<int> IPARAM(11), IPNTR(14);
	podarray<double> RESID(N), V(uword(N) * uword(NCV)), WORKD(5 * uword(N)), WORKL(LWORKL);

	IPARAM(0) = 1;    // exact shift
	IPARAM(2) = 2000; // maximum iteration
	IPARAM(6) = 4;    // mode 4: K*x=lambda*KG*x

	auto SIGMA = 1.;

	*KG *= -SIGMA;
	*KG += *K;

	while(99 != IDO) {
		arma_fortran(arma_dsaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
		// ReSharper disable once CppEntityAssignedButNoRead
		vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false);
		if(-1 == IDO) {
			vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
			X = *K * X;
			Y = KG->solve(X);
		} else if(1 == IDO) {
			const vec X(WORKD.memptr() + IPNTR[2] - 1, N, false);
			Y = KG->solve_trs(X);
		} else if(2 == IDO) {
			const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
			Y = *K * X;
		}
	}

	if(0 != INFO) return INFO;

	suanpan_debug("Arnoldi iteraiton counter: %d.\n", IPARAM(2));

	auto RVEC = 1;
	auto HOWMNY = 'A';
	auto LDZ = N;

	podarray<int> SELECT(NCV);

	eigval.set_size(NEV);
	eigvec.set_size(N, NEV);

	arma_fortran(arma_dseupd)(&RVEC, &HOWMNY, SELECT.memptr(), eigval.memptr(), eigvec.memptr(), &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

	return INFO;
}

int eig_solve(vec& eigval, mat& eigvec, const std::shared_ptr<MetaMat<double>>& K, const std::shared_ptr<MetaMat<double>>& KG) {
	auto IDO = 0;
	auto BMAT = 'G'; // generalized eigenvalue problem A*x=lambda*M*x
	auto N = static_cast<int>(K->n_cols);
	char WHICH[2] = {'L', 'M'};
	auto NEV = 1;
	auto TOL = 1E-10;
	auto NCV = std::min(std::max(NEV + 2, 2 * NEV) + 5, N);
	auto LDV = N;
	auto LWORKL = 2 * NCV * (NCV + 8);
	auto INFO = 0;

	podarray<int> IPARAM(11), IPNTR(14);
	podarray<double> RESID(N), V(uword(N) * uword(NCV)), WORKD(5 * uword(N)), WORKL(LWORKL);

	IPARAM(0) = 1;    // exact shift
	IPARAM(2) = 2000; // maximum iteration
	IPARAM(6) = 4;    // mode 4: K*x=lambda*KG*x

	auto SIGMA = 1.;

	// for buckling analysis do not multiply one.
	*KG += *K;

	while(99 != IDO) {
		arma_fortran(arma_dsaupd)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);
		// ReSharper disable once CppInitializedValueIsAlwaysRewritten
		// ReSharper disable once CppEntityAssignedButNoRead
		vec Y(WORKD.memptr() + IPNTR[1] - 1, N, false);
		if(-1 == IDO) {
			vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
			X = *K * X;
			Y = KG->solve(X);
		} else if(1 == IDO) {
			const vec X(WORKD.memptr() + IPNTR[2] - 1, N, false);
			Y = KG->solve_trs(X);
		} else if(2 == IDO) {
			const vec X(WORKD.memptr() + IPNTR[0] - 1, N, false);
			Y = *K * X;
		}
	}

	if(0 != INFO) return INFO;

	suanpan_debug("Arnoldi iteraiton counter: %d.\n", IPARAM(2));

	auto RVEC = 1;
	auto HOWMNY = 'A';
	auto LDZ = N;

	podarray<int> SELECT(NCV);

	eigval.set_size(NEV);
	eigvec.set_size(N, NEV);

	arma_fortran(arma_dseupd)(&RVEC, &HOWMNY, SELECT.memptr(), eigval.memptr(), eigvec.memptr(), &LDZ, &SIGMA, &BMAT, &N, WHICH, &NEV, &TOL, RESID.memptr(), &NCV, V.memptr(), &LDV, IPARAM.memptr(), IPNTR.memptr(), WORKD.memptr(), WORKL.memptr(), &LWORKL, &INFO);

	return INFO;
}
