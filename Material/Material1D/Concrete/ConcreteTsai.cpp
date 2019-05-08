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

#include "ConcreteTsai.h"

podarray<double> ConcreteTsai::compute_compression_backbone(const double t_strain) const {
	podarray<double> response(2);

	const auto normal_strain = std::max(datum::eps, t_strain / peak_strain);
	const auto tmp_a = pow(normal_strain, NC);
	const auto tmp_b = NC == 1. ? 1. + (MC - 1. + log(normal_strain)) * normal_strain : 1. + (MC - NC / (NC - 1.)) * normal_strain + tmp_a / (NC - 1.);
	response(0) = peak_stress * MC * normal_strain / tmp_b;
	response(1) = initial_stiffness(0) * (1. - tmp_a) / tmp_b / tmp_b;

	return response;
}

podarray<double> ConcreteTsai::compute_tension_backbone(const double t_strain) const {
	podarray<double> response(2);

	const auto normal_strain = std::max(datum::eps, t_strain / crack_strain);
	const auto tmp_a = pow(normal_strain, NT);
	const auto tmp_b = NT == 1. ? 1. + (MT - 1. + log(normal_strain)) * normal_strain : 1. + (MT - NT / (NT - 1.)) * normal_strain + tmp_a / (NT - 1.);
	response(0) = crack_stress * MT * normal_strain / tmp_b;
	response(1) = initial_stiffness(0) * (1. - tmp_a) / tmp_b / tmp_b;

	return response;
}

double ConcreteTsai::compute_compression_residual(const double reverse_c_strain, const double reverse_c_stress) const { return std::min(0., reverse_c_strain - reverse_c_stress * (reverse_c_strain / peak_strain + .57) / (reverse_c_stress / peak_strain + .57 * initial_stiffness(0))); }

double ConcreteTsai::compute_tension_residual(const double reverse_t_strain, const double reverse_t_stress) const { return std::max(0., reverse_t_strain - reverse_t_stress * (reverse_t_strain / crack_strain + .67) / (reverse_t_stress / crack_strain + .67 * initial_stiffness(0))); }

ConcreteTsai::ConcreteTsai(const unsigned T, const double CS, const double TS, const double MCC, const double NCC, const double MTT, const double NTT, const double MP, const double CE, const double TE, const double R)
	: ConcreteSimple(T, CS, TS, CE, TE, MP, R)
	, MC(std::max(1.1, MCC))
	, MT(std::max(1.1, MTT))
	, NC(std::max(1.1 * MC / (MC - 1.), NCC))
	, NT(std::max(1.1 * MT / (MT - 1.), NTT)) {}

ConcreteTsai::ConcreteTsai(const unsigned T, const double CS, const double TS, const double MCT, const double NCT, const double MP, const double CE, const double TE, const double R)
	: ConcreteTsai(T, CS, TS, MCT, NCT, MCT, NCT, MP, CE, TE, R) {}

void ConcreteTsai::initialize(const shared_ptr<DomainBase>&) {
	trial_stiffness = current_stiffness = initial_stiffness = std::max(crack_stress / crack_strain * MT, peak_stress / peak_strain * MC);
	access::rw(MC) = initial_stiffness(0) * peak_strain / peak_stress;
	access::rw(MT) = initial_stiffness(0) * crack_strain / crack_stress;
	access::rw(NC) = std::max(1.1 * MC / (MC - 1.), NC);
	access::rw(NT) = std::max(1.1 * MT / (MT - 1.), NT);
}

unique_ptr<Material> ConcreteTsai::get_copy() { return make_unique<ConcreteTsai>(*this); }
