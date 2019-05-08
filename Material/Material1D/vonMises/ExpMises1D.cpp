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

#include "ExpMises1D.h"

double ExpMises1D::compute_k(const double p_strain) const {
	const auto b_strain = -b * p_strain;
	return yield_stress * ((1. + a) * exp(b_strain) - a * exp(2. * b_strain));
}

double ExpMises1D::compute_dk(const double p_strain) const {
	const auto b_strain = -b * p_strain;
	return yield_stress * b * (2. * a * exp(2. * b_strain) - (1. + a) * exp(b_strain));
}

double ExpMises1D::compute_h(const double) const { return 0.; }

double ExpMises1D::compute_dh(const double) const { return 0.; }

ExpMises1D::ExpMises1D(const unsigned T, const double E, const double Y, const double A, const double B, const double R)
	: DataExpMises1D{fabs(Y), fabs(A), fabs(B)}
	, NonlinearMises1D(T, E, R) {}

unique_ptr<Material> ExpMises1D::get_copy() { return make_unique<ExpMises1D>(*this); }
