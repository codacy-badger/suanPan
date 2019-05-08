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

#include "ConcreteCM.h"
#include <Toolbox/utility.h>
#include <Recorder/OutputType.h>

/*
 *  auto& max_c_strain = trial_history(0);               // maximum compression strain recorded
 *  auto& max_t_strain = trial_history(1);               // maximum tension strain recorded
 *  auto& reverse_c_strain = trial_history(2);           // unloading point strain compression
 *  auto& reverse_c_stress = trial_history(3);           // unloading point stress compression
 *  auto& reverse_t_strain = trial_history(4);           // unloading point strain tension
 *  auto& reverse_t_stress = trial_history(5);           // unloading point stress tension
 *  auto& residual_c_strain = trial_history(6);          // residual strain in compression unloading path
 *  auto& residual_c_stiffness = trial_history(7);       // tangent stiffness at residual
 *  auto& residual_t_strain = trial_history(8);          // residual strain in compression unloading path
 *  auto& residual_t_stiffness = trial_history(9);       // tangent stiffness at residual
 *  auto& reload_c_strain = trial_history(10);           // reload point strain compression
 *  auto& reload_c_stress = trial_history(11);           // reload point stress compression
 *  auto& reload_c_stiffness = trial_history(12);        // reload point stiffness compression
 *  auto& reload_t_strain = trial_history(13);           // reload point strain tension
 *  auto& reload_t_stress = trial_history(14);           // reload point stress tension
 *  auto& reload_t_stiffness = trial_history(15);        // reload point stiffness tension
 *  auto& reverse_below_c_stress = trial_history(16);    // unloading point stress compression
 *  auto& reverse_below_c_stiffness = trial_history(17); // unloading point stiffness compression
 *  auto& reverse_below_t_stress = trial_history(18);    // unloading point stress tension
 *  auto& reverse_below_t_stiffness = trial_history(19); // unloading point stiffness tension
 *  auto& trans_c_strain = trial_history(20);            // transition strain compression
 *  auto& trans_c_stress = trial_history(21);            // transition stress compression
 *  auto& trans_c_stiffness = trial_history(22);         // transition stiffness compression
 *  auto& trans_t_strain = trial_history(23);            // transition strain tension
 *  auto& trans_t_stress = trial_history(24);            // transition stress tension
 *  auto& trans_t_stiffness = trial_history(25);         // transition stiffness tension
 *  auto& new_reverse_c_strain = trial_history(26);      // new reverse strain induced by partial reloading compression
 *  auto& new_reverse_c_stress = trial_history(27);      // new reverse stress induced by partial reloading compression
 *  auto& new_reverse_t_strain = trial_history(28);      // new reverse strain induced by partial reloading tension
 *  auto& new_reverse_t_stress = trial_history(29);      // new reverse stress induced by partial reloading tension
 *  auto& new_reverse_c_stiffness = trial_history(30);   // new reverse stiffness induced by partial reloading compression
 *  auto& new_reverse_t_stiffness = trial_history(31);   // new reverse stiffness induced by partial reloading tension
 */

ConcreteCM::ConcreteCM(const unsigned T, const double SC, const double ST, const double MCC, const double NCC, const double MTT, const double NTT, const double EC, const double ET, const bool LT, const double R)
	: Material1D(T, R)
	, peak_stress(-fabs(SC))
	, peak_strain(-fabs(EC))
	, crack_stress(fabs(ST))
	, crack_strain(fabs(ET))
	, MC(MCC)
	, NC(NCC)
	, MT(MTT)
	, NT(NTT)
	, linear_trans(LT) {}

void ConcreteCM::initialize(const shared_ptr<DomainBase>&) {
	trial_stiffness = current_stiffness = initial_stiffness = std::max(crack_stress / crack_strain * MT, peak_stress / peak_strain * MC);

	access::rw(MT) = crack_strain / crack_stress * initial_stiffness(0);
	access::rw(MC) = peak_strain / peak_stress * initial_stiffness(0);

	trial_history = current_history.zeros(32);
}

unique_ptr<Material> ConcreteCM::get_copy() { return make_unique<ConcreteCM>(*this); }

double ConcreteCM::get_parameter(const ParameterType P) const {
	switch(P) {
	case ParameterType::DENSITY:
		return density;
	case ParameterType::ELASTICMODULUS:
	case ParameterType::YOUNGSMODULUS:
	case ParameterType::E:
		return initial_stiffness(0);
	case ParameterType::PEAKSTRAIN:
		return peak_strain;
	case ParameterType::CRACKSTRAIN:
		return crack_strain;
	default:
		return 0.;
	}
}

int ConcreteCM::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

	trial_load_status = current_load_status;

	trial_history = current_history;
	const auto& reverse_c_strain = trial_history(2);           // unloading point strain compression
	const auto& reverse_t_strain = trial_history(4);           // unloading point strain tension
	const auto& residual_c_strain = trial_history(6);          // residual strain in compression unloading path
	const auto& residual_t_strain = trial_history(8);          // residual strain in compression unloading path
	const auto& reverse_below_c_stress = trial_history(16);    // unloading point stress compression
	const auto& reverse_below_c_stiffness = trial_history(17); // unloading point stiffness compression
	const auto& reverse_below_t_stress = trial_history(18);    // unloading point stress tension
	const auto& reverse_below_t_stiffness = trial_history(19); // unloading point stiffness tension
	auto& max_c_strain = trial_history(0);                     // maximum compression strain recorded
	auto& max_t_strain = trial_history(1);                     // maximum tension strain recorded
	auto& residual_c_stiffness = trial_history(7);             // tangent stiffness at residual
	auto& residual_t_stiffness = trial_history(9);             // tangent stiffness at residual
	auto& trans_c_strain = trial_history(20);                  // transition strain compression
	auto& trans_c_stress = trial_history(21);                  // transition stress compression
	auto& trans_c_stiffness = trial_history(22);               // transition stiffness compression
	auto& trans_t_strain = trial_history(23);                  // transition strain tension
	auto& trans_t_stress = trial_history(24);                  // transition stress tension
	auto& trans_t_stiffness = trial_history(25);               // transition stiffness tension
	auto& new_reverse_c_strain = trial_history(26);            // new reverse strain induced by partial reloading compression
	auto& new_reverse_c_stress = trial_history(27);            // new reverse stress induced by partial reloading compression
	auto& new_reverse_t_strain = trial_history(28);            // new reverse strain induced by partial reloading tension
	auto& new_reverse_t_stress = trial_history(29);            // new reverse stress induced by partial reloading tension
	auto& new_reverse_c_stiffness = trial_history(30);         // new reverse stiffness induced by partial reloading compression
	auto& new_reverse_t_stiffness = trial_history(31);         // new reverse stiffness induced by partial reloading tension

	const auto load_direction = int(suanpan::sign(incre_strain(0)));

	podarray<double> response(2);

	switch(trial_load_status) {
	case Status::NONE:
		if(-1 == load_direction) {
			// compression backbone
			trial_load_status = Status::CBACKBONE;
			response = compute_compression_backbone(max_c_strain = trial_strain(0));
		} else {
			// tension backbone
			trial_load_status = Status::TBACKBONE;
			response = compute_tension_backbone(max_t_strain = trial_strain(0));
		}
		break;
	case Status::CBACKBONE:
		if(1 == load_direction) {
			// unload from compression backbone
			update_compression_reverse(current_strain(0));
			response = compute_compression_unload(trial_strain(0));
		} else response = compute_compression_backbone(max_c_strain = trial_strain(0));
		break;
	case Status::TBACKBONE:
		if(-1 == load_direction) {
			// unload from tension backbone
			update_tension_reverse(current_strain(0));
			response = compute_tension_unload(trial_strain(0));
		} else response = compute_tension_backbone(max_t_strain = trial_strain(0));
		break;
	case Status::CUNLOAD:
		response = -1 == load_direction ? compute_compression_reload(trial_strain(0)) : compute_compression_unload(trial_strain(0));
		break;
	case Status::TUNLOAD:
		response = 1 == load_direction ? compute_tension_reload(trial_strain(0)) : compute_tension_unload(trial_strain(0));
		break;
	case Status::CSUBUNLOAD:
		response = -1 == load_direction ? compute_compression_reload(trial_strain(0)) : compute_compression_sub_unload(trial_strain(0));
		break;
	case Status::TSUBUNLOAD:
		response = 1 == load_direction ? compute_tension_reload(trial_strain(0)) : compute_tension_sub_unload(trial_strain(0));
		break;
	case Status::CRELOAD:
		if(1 == load_direction) {
			// new unload from compression
			if(current_strain(0) < reverse_c_strain) update_compression_reverse(current_strain(0));
			new_reverse_c_strain = current_strain(0);
			new_reverse_c_stress = current_stress(0);
			const auto secant_stiffness = new_reverse_c_stress / (new_reverse_c_strain - residual_c_strain);
			residual_c_stiffness = std::min(residual_c_stiffness, .9 * secant_stiffness);
			new_reverse_c_stiffness = std::max(initial_stiffness(0), 1.1 * secant_stiffness);
			response = compute_compression_sub_unload(trial_strain(0));
		} else response = compute_compression_reload(trial_strain(0));
		break;
	case Status::TRELOAD:
		if(-1 == load_direction) {
			// new unload from tension
			if(current_strain(0) > reverse_t_strain) update_tension_reverse(current_strain(0));
			new_reverse_t_strain = current_strain(0);
			new_reverse_t_stress = current_stress(0);
			const auto secant_stiffness = new_reverse_t_stress / (new_reverse_t_strain - residual_t_strain);
			residual_t_stiffness = std::min(residual_t_stiffness, .9 * secant_stiffness);
			new_reverse_t_stiffness = std::max(initial_stiffness(0), 1.1 * secant_stiffness);
			response = compute_tension_sub_unload(trial_strain(0));
		} else response = compute_tension_reload(trial_strain(0));
		break;
	case Status::CTRANS:
		if(trial_strain(0) < trans_c_strain) response = compute_tension_unload(trial_strain(0));
		else if(trial_strain(0) > trans_t_strain) response = compute_compression_unload(trial_strain(0));
		else if(-1 == load_direction) response = compute_transition(trial_strain(0), trans_c_strain, trans_c_stress, trans_c_stiffness, trans_t_strain, trans_t_stress, trans_t_stiffness, linear_trans);
		else {
			trial_load_status = Status::TTRANS;
			trans_c_strain = current_strain(0);
			trans_c_stress = current_stress(0);
			response = compute_transition(trans_t_strain = residual_t_strain, residual_c_strain, 0., residual_c_stiffness, reverse_t_strain, reverse_below_t_stress, reverse_below_t_stiffness, linear_trans);
			trans_t_stress = response(0);
			const auto secant_stiffness = (trans_t_stress - trans_c_stress) / (trans_t_strain - trans_c_strain);
			trans_t_stiffness = std::min(response(1), .9 * secant_stiffness);
			trans_c_stiffness = std::max(initial_stiffness(0), 1.1 * secant_stiffness);
			response = compute_transition(trial_strain(0), trans_t_strain, trans_t_stress, trans_t_stiffness, trans_c_strain, trans_c_stress, trans_c_stiffness, linear_trans);
		}
		break;
	case Status::TTRANS:
		if(trial_strain(0) < trans_c_strain) response = compute_tension_unload(trial_strain(0));
		else if(trial_strain(0) > trans_t_strain) response = compute_compression_unload(trial_strain(0));
		else if(1 == load_direction) response = compute_transition(trial_strain(0), trans_t_strain, trans_t_stress, trans_t_stiffness, trans_c_strain, trans_c_stress, trans_c_stiffness, linear_trans);
		else {
			trial_load_status = Status::CTRANS;
			trans_t_strain = current_strain(0);
			trans_t_stress = current_stress(0);
			response = compute_transition(trans_c_strain = residual_c_strain, residual_t_strain, 0., residual_t_stiffness, reverse_c_strain, reverse_below_c_stress, reverse_below_c_stiffness, linear_trans);
			trans_c_stress = response(0);
			const auto secant_stiffness = (trans_t_stress - trans_c_stress) / (trans_t_strain - trans_c_strain);
			trans_c_stiffness = std::min(response(1), .9 * secant_stiffness);
			trans_t_stiffness = std::max(initial_stiffness(0), 1.1 * secant_stiffness);
			response = compute_transition(trial_strain(0), trans_c_strain, trans_c_stress, trans_c_stiffness, trans_t_strain, trans_t_stress, trans_t_stiffness, linear_trans);
		}
		break;
	}

	trial_stress = response(0);
	trial_stiffness = response(1);

	suanpan_debug([&]() { if(!trial_stress.is_finite() || !trial_stiffness.is_finite()) throw invalid_argument("infinite number detected.\n"); });

	return SUANPAN_SUCCESS;
}

int ConcreteCM::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
	if(update_trial_status(t_strain) == SUANPAN_SUCCESS) {
		const auto factor = 35E-3 * peak_stress * peak_stress;
		const auto magnification_factor = (1. + pow(fabs(t_strain_rate(0) / factor), 1. / 6.)) / (1. + pow(fabs(1E-5 / factor), 1. / 6.));
		if(std::isfinite(magnification_factor)) {
			trial_stress *= magnification_factor;
			trial_stiffness *= magnification_factor;
		}
		return SUANPAN_SUCCESS;
	}
	return SUANPAN_FAIL;
}

int ConcreteCM::clear_status() {
	current_strain.zeros();
	current_stress.zeros();
	current_history.zeros();
	current_load_status = Status::NONE;
	current_stiffness = initial_stiffness;
	return reset_status();
}

int ConcreteCM::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_history = trial_history;
	current_load_status = trial_load_status;
	current_stiffness = trial_stiffness;
	return SUANPAN_SUCCESS;
}

int ConcreteCM::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_history = current_history;
	trial_load_status = current_load_status;
	trial_stiffness = current_stiffness;
	return SUANPAN_SUCCESS;
}

vector<vec> ConcreteCM::record(const OutputType P) { return Material1D::record(P); }

void ConcreteCM::print() { suanpan_info("A concrete model based on Chang & Mander's concrete model.\n"); }

podarray<double> ConcreteCM::compute_compression_backbone(const double t_strain) const {
	podarray<double> status(2);

	const auto normal_strain = std::max(datum::eps, t_strain / peak_strain);

	suanpan_debug([&]() { if(normal_strain < 0.) throw invalid_argument("need positive normalised strain"); });

	const auto tmp_a = pow(normal_strain, NC);
	const auto tmp_b = NC == 1. ? 1. + (MC - 1. + log(normal_strain)) * normal_strain : 1. + (MC - NC / (NC - 1.)) * normal_strain + tmp_a / (NC - 1.);
	status(0) = peak_stress * MC * normal_strain / tmp_b;
	status(1) = initial_stiffness(0) * (1. - tmp_a) / tmp_b / tmp_b;

	return status;
}

podarray<double> ConcreteCM::compute_tension_backbone(const double t_strain) const {
	podarray<double> status(2);

	const auto normal_strain = std::max(datum::eps, t_strain / crack_strain);

	suanpan_debug([&]() { if(normal_strain < 0.) throw invalid_argument("need positive normalised strain"); });

	const auto tmp_a = pow(normal_strain, NT);
	const auto tmp_b = NT == 1. ? 1. + (MT - 1. + log(normal_strain)) * normal_strain : 1. + (MT - NT / (NT - 1.)) * normal_strain + tmp_a / (NT - 1.);
	status(0) = crack_stress * MT * normal_strain / tmp_b;
	status(1) = initial_stiffness(0) * (1. - tmp_a) / tmp_b / tmp_b;

	return status;
}

podarray<double> ConcreteCM::compute_compression_unload(const double t_strain) {
	const auto& reverse_c_strain = trial_history(2);           // unloading point strain compression
	const auto& reverse_c_stress = trial_history(3);           // unloading point stress compression
	const auto& reverse_t_strain = trial_history(4);           // unloading point strain tension
	const auto& residual_c_strain = trial_history(6);          // residual strain in compression unloading path
	const auto& residual_c_stiffness = trial_history(7);       // tangent stiffness at residual
	const auto& reload_t_strain = trial_history(13);           // reload point strain tension
	const auto& reload_t_stress = trial_history(14);           // reload point stress tension
	const auto& reload_t_stiffness = trial_history(15);        // reload point stiffness tension
	const auto& reverse_below_t_stress = trial_history(18);    // unloading point stress tension
	const auto& reverse_below_t_stiffness = trial_history(19); // unloading point stiffness tension
	auto& max_t_strain = trial_history(1);                     // maximum tension strain recorded

	podarray<double> response;

	trial_load_status = Status::CUNLOAD;

	if(t_strain <= residual_c_strain)
		// rule 3
		response = compute_transition(t_strain, residual_c_strain, 0., residual_c_stiffness, reverse_c_strain, reverse_c_stress, initial_stiffness(0), linear_trans);
	else if(t_strain <= reverse_t_strain)
		// rule 9
		response = compute_transition(t_strain, residual_c_strain, 0., residual_c_stiffness, reverse_t_strain, reverse_below_t_stress, reverse_below_t_stiffness, linear_trans);
	else if(t_strain <= reload_t_strain)
		// rule 8
		response = compute_transition(t_strain, reload_t_strain, reload_t_stress, reload_t_stiffness, reverse_t_strain, reverse_below_t_stress, reverse_below_t_stiffness, linear_trans);
	else {
		// tension backbone
		trial_load_status = Status::TBACKBONE;
		response = compute_tension_backbone(max_t_strain = t_strain);
	}

	return response;
}

podarray<double> ConcreteCM::compute_tension_unload(const double t_strain) {
	const auto& reverse_c_strain = trial_history(2);           // unloading point strain compression
	const auto& reverse_t_strain = trial_history(4);           // unloading point strain tension
	const auto& reverse_t_stress = trial_history(5);           // unloading point stress tension
	const auto& residual_t_strain = trial_history(8);          // residual strain in compression unloading path
	const auto& residual_t_stiffness = trial_history(9);       // tangent stiffness at residual
	const auto& reload_c_strain = trial_history(10);           // reload point strain compression
	const auto& reload_c_stress = trial_history(11);           // reload point stress compression
	const auto& reload_c_stiffness = trial_history(12);        // reload point stiffness compression
	const auto& reverse_below_c_stress = trial_history(16);    // unloading point stress compression
	const auto& reverse_below_c_stiffness = trial_history(17); // unloading point stiffness compression
	auto& max_c_strain = trial_history(0);                     // maximum compression strain recorded

	podarray<double> response;

	trial_load_status = Status::TUNLOAD;

	if(t_strain >= residual_t_strain)
		// rule 4
		response = compute_transition(t_strain, residual_t_strain, 0., residual_t_stiffness, reverse_t_strain, reverse_t_stress, initial_stiffness(0), linear_trans);
	else if(t_strain >= reverse_c_strain)
		// rule 10
		response = compute_transition(t_strain, residual_t_strain, 0., residual_t_stiffness, reverse_c_strain, reverse_below_c_stress, reverse_below_c_stiffness, linear_trans);
	else if(t_strain >= reload_c_strain) {
		// rule 7
		response = compute_transition(t_strain, reload_c_strain, reload_c_stress, reload_c_stiffness, reverse_c_strain, reverse_below_c_stress, reverse_below_c_stiffness, linear_trans);
	} else {
		// compression backbone
		trial_load_status = Status::CBACKBONE;
		response = compute_compression_backbone(max_c_strain = t_strain);
	}

	return response;
}

podarray<double> ConcreteCM::compute_compression_sub_unload(const double t_strain) {
	const auto& residual_c_strain = trial_history(6);        // residual strain in compression unloading path
	const auto& residual_c_stiffness = trial_history(7);     // tangent stiffness at residual
	const auto& new_reverse_c_strain = trial_history(26);    // new reverse strain induced by partial reloading compression
	const auto& new_reverse_c_stress = trial_history(27);    // new reverse stress induced by partial reloading compression
	const auto& new_reverse_c_stiffness = trial_history(30); // new reverse stiffness induced by partial reloading compression

	trial_load_status = Status::CSUBUNLOAD;

	return t_strain <= residual_c_strain ? compute_transition(t_strain, residual_c_strain, 0., residual_c_stiffness, new_reverse_c_strain, new_reverse_c_stress, new_reverse_c_stiffness, linear_trans) : compute_compression_unload(t_strain);
}

podarray<double> ConcreteCM::compute_tension_sub_unload(const double t_strain) {
	const auto& residual_t_strain = trial_history(8);        // residual strain in compression unloading path
	const auto& residual_t_stiffness = trial_history(9);     // tangent stiffness at residual
	const auto& new_reverse_t_strain = trial_history(28);    // new reverse strain induced by partial reloading tension
	const auto& new_reverse_t_stress = trial_history(29);    // new reverse stress induced by partial reloading tension
	const auto& new_reverse_t_stiffness = trial_history(31); // new reverse stiffness induced by partial reloading tension

	trial_load_status = Status::TSUBUNLOAD;

	return t_strain >= residual_t_strain ? compute_transition(t_strain, residual_t_strain, 0., residual_t_stiffness, new_reverse_t_strain, new_reverse_t_stress, new_reverse_t_stiffness, linear_trans) : compute_tension_unload(t_strain);
}

podarray<double> ConcreteCM::compute_compression_reload(const double t_strain) {
	const auto& reverse_c_strain = trial_history(2);        // unloading point strain compression
	const auto& reverse_t_strain = trial_history(4);        // unloading point strain tension
	const auto& residual_c_strain = trial_history(6);       // residual strain in compression unloading path
	const auto& residual_t_strain = trial_history(8);       // residual strain in compression unloading path
	const auto& reload_c_strain = trial_history(10);        // reload point strain compression
	const auto& reload_c_stress = trial_history(11);        // reload point stress compression
	const auto& reload_c_stiffness = trial_history(12);     // reload point stiffness compression
	const auto& reverse_below_c_stress = trial_history(16); // unloading point stress compression
	auto& max_c_strain = trial_history(0);                  // maximum compression strain recorded
	auto& residual_t_stiffness = trial_history(9);          // tangent stiffness at residual
	auto& reverse_below_c_stiffness = trial_history(17);    // unloading point stiffness compression
	auto& trans_c_strain = trial_history(20);               // transition strain compression
	auto& trans_c_stress = trial_history(21);               // transition stress compression
	auto& trans_c_stiffness = trial_history(22);            // transition stiffness compression
	auto& trans_t_strain = trial_history(23);               // transition strain tension
	auto& trans_t_stress = trial_history(24);               // transition stress tension
	auto& trans_t_stiffness = trial_history(25);            // transition stiffness tension
	auto& new_reverse_t_strain = trial_history(28);         // new reverse strain induced by partial reloading tension
	auto& new_reverse_t_stress = trial_history(29);         // new reverse stress induced by partial reloading tension
	auto& new_reverse_t_stiffness = trial_history(31);      // new reverse stiffness induced by partial reloading tension

	podarray<double> status(2);

	if(current_strain(0) > .5 * (residual_t_strain + reverse_t_strain)) {
		if(current_strain(0) > reverse_t_strain) update_tension_reverse(current_strain(0));
		new_reverse_t_strain = current_strain(0);
		new_reverse_t_stress = current_stress(0);
		const auto secant_stiffness = new_reverse_t_stress / (new_reverse_t_strain - residual_t_strain);
		residual_t_stiffness = std::min(residual_t_stiffness, .9 * secant_stiffness);
		new_reverse_t_stiffness = std::max(initial_stiffness(0), 1.1 * secant_stiffness);
		status = compute_tension_sub_unload(t_strain);
	} else if(current_strain(0) > residual_c_strain) {
		trans_t_strain = current_strain(0);
		trans_t_stress = current_stress(0);
		auto secant_stiffness = reverse_below_c_stress / (reverse_c_strain - residual_t_strain);
		residual_t_stiffness = std::min(residual_t_stiffness, .9 * secant_stiffness);
		reverse_below_c_stiffness = std::max(reverse_below_c_stiffness, 1.1 * secant_stiffness);
		status = compute_transition(trans_c_strain = residual_c_strain, residual_t_strain, 0., residual_t_stiffness, reverse_c_strain, reverse_below_c_stress, reverse_below_c_stiffness, linear_trans);
		secant_stiffness = (trans_t_stress - (trans_c_stress = status(0))) / (trans_t_strain - trans_c_strain);
		trans_c_stiffness = std::min(status(1), .9 * secant_stiffness);
		trans_t_stiffness = std::max(initial_stiffness(0), 1.1 * secant_stiffness);
		if(trial_strain(0) > residual_c_strain) {
			trial_load_status = Status::CTRANS;
			status = compute_transition(t_strain, trans_c_strain, trans_c_stress, trans_c_stiffness, trans_t_strain, trans_t_stress, trans_t_stiffness, linear_trans);
		} else status = compute_tension_unload(t_strain);
	} else if(current_strain(0) > reverse_c_strain) {
		if(trial_load_status != Status::CRELOAD) update_compression_new_reverse(current_strain(0), current_stress(0));
		if(trial_strain(0) > reverse_c_strain) {
			// linear reloading
			trial_load_status = Status::CRELOAD;
			status(0) = reverse_below_c_stress + (status(1) = reverse_below_c_stiffness) * (t_strain - reverse_c_strain);
		} else if(trial_strain(0) > reload_c_strain) {
			// rule 7 (new)
			trial_load_status = Status::TUNLOAD;
			status = compute_transition(t_strain, reverse_c_strain, reverse_below_c_stress, reverse_below_c_stiffness, reload_c_strain, reload_c_stress, reload_c_stiffness, linear_trans);
		} else {
			// compression backbone
			trial_load_status = Status::CBACKBONE;
			status = compute_compression_backbone(max_c_strain = t_strain);
		}
	}

	return status;
}

podarray<double> ConcreteCM::compute_tension_reload(const double t_strain) {
	const auto& reverse_c_strain = trial_history(2);        // unloading point strain compression
	const auto& reverse_t_strain = trial_history(4);        // unloading point strain tension
	const auto& residual_c_strain = trial_history(6);       // residual strain in compression unloading path
	const auto& residual_t_strain = trial_history(8);       // residual strain in compression unloading path
	const auto& reload_t_strain = trial_history(13);        // reload point strain tension
	const auto& reload_t_stress = trial_history(14);        // reload point stress tension
	const auto& reload_t_stiffness = trial_history(15);     // reload point stiffness tension
	const auto& reverse_below_t_stress = trial_history(18); // unloading point stress tension
	auto& max_t_strain = trial_history(1);                  // maximum tension strain recorded
	auto& residual_c_stiffness = trial_history(7);          // tangent stiffness at residual
	auto& reverse_below_t_stiffness = trial_history(19);    // unloading point stiffness tension
	auto& trans_c_strain = trial_history(20);               // transition strain compression
	auto& trans_c_stress = trial_history(21);               // transition stress compression
	auto& trans_c_stiffness = trial_history(22);            // transition stiffness compression
	auto& trans_t_strain = trial_history(23);               // transition strain tension
	auto& trans_t_stress = trial_history(24);               // transition stress tension
	auto& trans_t_stiffness = trial_history(25);            // transition stiffness tension
	auto& new_reverse_c_strain = trial_history(26);         // new reverse strain induced by partial reloading compression
	auto& new_reverse_c_stress = trial_history(27);         // new reverse stress induced by partial reloading compression
	auto& new_reverse_c_stiffness = trial_history(30);      // new reverse stiffness induced by partial reloading compression

	podarray<double> status(2);

	if(current_strain(0) < .5 * (residual_c_strain + reverse_c_strain)) {
		if(current_strain(0) < reverse_c_strain) update_compression_reverse(current_strain(0));
		new_reverse_c_strain = current_strain(0);
		new_reverse_c_stress = current_stress(0);
		const auto secant_stiffness = new_reverse_c_stress / (new_reverse_c_strain - residual_c_strain);
		residual_c_stiffness = std::min(residual_c_stiffness, .9 * secant_stiffness);
		new_reverse_c_stiffness = std::max(initial_stiffness(0), 1.1 * secant_stiffness);
		status = compute_compression_sub_unload(t_strain);
	} else if(current_strain(0) < residual_t_strain) {
		trans_c_strain = current_strain(0);
		trans_c_stress = current_stress(0);
		auto secant_stiffness = reverse_below_t_stress / (reverse_t_strain - residual_c_strain);
		residual_c_stiffness = std::max(residual_c_stiffness, 1.1 * secant_stiffness);
		reverse_below_t_stiffness = std::min(reverse_below_t_stiffness, .9 * secant_stiffness);
		status = compute_transition(trans_t_strain = residual_t_strain, residual_c_strain, 0., residual_c_stiffness, reverse_t_strain, reverse_below_t_stress, reverse_below_t_stiffness, linear_trans);
		secant_stiffness = ((trans_t_stress = status(0)) - trans_c_stress) / (trans_t_strain - trans_c_strain);
		trans_t_stiffness = std::min(status(1), .9 * secant_stiffness);
		trans_c_stiffness = std::max(initial_stiffness(0), 1.1 * secant_stiffness);
		if(trial_strain(0) < residual_t_strain) {
			trial_load_status = Status::TTRANS;
			status = compute_transition(t_strain, trans_t_strain, trans_t_stress, trans_t_stiffness, trans_c_strain, trans_c_stress, trans_c_stiffness, linear_trans);
		} else status = compute_compression_unload(t_strain);
	} else if(current_strain(0) < reverse_t_strain) {
		if(trial_load_status != Status::TRELOAD) update_tension_new_reverse(current_strain(0), current_stress(0));
		if(trial_strain(0) < reverse_t_strain) {
			trial_load_status = Status::TRELOAD;
			status(0) = reverse_below_t_stress + (status(1) = reverse_below_t_stiffness) * (t_strain - reverse_t_strain);
		} else if(trial_strain(0) < reload_t_strain) {
			trial_load_status = Status::CUNLOAD;
			status = compute_transition(t_strain, reverse_t_strain, reverse_below_t_stress, reverse_below_t_stiffness, reload_t_strain, reload_t_stress, reload_t_stiffness, linear_trans);
		} else {
			trial_load_status = Status::TBACKBONE;
			status = compute_tension_backbone(max_t_strain = t_strain);
		}
	}

	return status;
}

void ConcreteCM::update_compression_reverse(const double c_strain) {
	const auto& reverse_t_strain = trial_history(4);           // unloading point strain tension
	const auto& reverse_below_t_stress = trial_history(18);    // unloading point stress tension
	const auto& reverse_below_t_stiffness = trial_history(19); // unloading point stiffness tension
	auto& reverse_c_strain = trial_history(2);                 // unloading point strain compression
	auto& reverse_c_stress = trial_history(3);                 // unloading point stress compression
	auto& residual_c_strain = trial_history(6);                // residual strain in compression unloading path
	auto& residual_c_stiffness = trial_history(7);             // tangent stiffness at residual
	auto& residual_t_strain = trial_history(8);                // residual strain in compression unloading path
	auto& residual_t_stiffness = trial_history(9);             // tangent stiffness at residual
	auto& reload_c_strain = trial_history(10);                 // reload point strain compression
	auto& reload_c_stress = trial_history(11);                 // reload point stress compression
	auto& reload_c_stiffness = trial_history(12);              // reload point stiffness compression
	auto& reverse_below_c_stress = trial_history(16);          // unloading point stress compression
	auto& reverse_below_c_stiffness = trial_history(17);       // unloading point stiffness compression

	auto response = compute_compression_backbone(reverse_c_strain = c_strain);
	reverse_c_stress = response(0);
	const auto normalized_strain = fabs(reverse_c_strain / peak_strain);
	const auto secant_stiffness = (fabs(reverse_c_stress / peak_strain) + .57 * initial_stiffness(0)) / (normalized_strain + .57);
	residual_c_strain = std::min(0., reverse_c_strain - reverse_c_stress / secant_stiffness);
	residual_c_stiffness = std::min(.9 * secant_stiffness, .5 * initial_stiffness(0) * exp(-2. * normalized_strain));

	response = compute_compression_backbone(reload_c_strain = std::min(reverse_c_strain + reverse_c_strain / (1.15 + 2.75 * normalized_strain), 1.05 * peak_strain));
	reload_c_stress = response(0);
	reload_c_stiffness = response(1);
	reverse_below_c_stress = std::max(reverse_c_stress * (1. - .09 * sqrt(normalized_strain)), reload_c_stress + .9 * reload_c_stiffness * (reverse_c_strain - reload_c_strain));
	reverse_below_c_stiffness = reverse_below_c_stress / (reverse_c_strain - residual_c_strain);

	if(reverse_below_t_stiffness == 0.) update_tension_reverse(1.05 * crack_strain);

	residual_t_stiffness = std::min(residual_t_stiffness, .9 * reverse_below_c_stress / (reverse_c_strain - residual_t_strain));
	residual_c_stiffness = std::min(residual_c_stiffness, .9 * reverse_below_t_stress / (reverse_t_strain - residual_c_strain));
}

void ConcreteCM::update_tension_reverse(const double c_strain) {
	const auto& reverse_c_strain = trial_history(2);           // unloading point strain compression
	const auto& reverse_below_c_stress = trial_history(16);    // unloading point stress compression
	const auto& reverse_below_c_stiffness = trial_history(17); // unloading point stiffness compression
	auto& reverse_t_strain = trial_history(4);                 // unloading point strain tension
	auto& reverse_t_stress = trial_history(5);                 // unloading point stress tension
	auto& residual_c_strain = trial_history(6);                // residual strain in compression unloading path
	auto& residual_c_stiffness = trial_history(7);             // tangent stiffness at residual
	auto& residual_t_strain = trial_history(8);                // residual strain in compression unloading path
	auto& residual_t_stiffness = trial_history(9);             // tangent stiffness at residual
	auto& reload_t_strain = trial_history(13);                 // reload point strain tension
	auto& reload_t_stress = trial_history(14);                 // reload point stress tension
	auto& reload_t_stiffness = trial_history(15);              // reload point stiffness tension
	auto& reverse_below_t_stress = trial_history(18);          // unloading point stress tension
	auto& reverse_below_t_stiffness = trial_history(19);       // unloading point stiffness tension

	auto response = compute_tension_backbone(reverse_t_strain = c_strain);
	reverse_t_stress = response(0);
	const auto normalized_strain = fabs(reverse_t_strain / crack_strain);
	const auto secant_stiffness = (fabs(reverse_t_stress / crack_strain) + .67 * initial_stiffness(0)) / (normalized_strain + .67);
	residual_t_strain = std::max(0., reverse_t_strain - reverse_t_stress / secant_stiffness);
	residual_t_stiffness = std::min(.9 * secant_stiffness, initial_stiffness(0) / (pow(normalized_strain, 1.1) + 1.));

	response = compute_tension_backbone(reload_t_strain = std::max(1.05 * crack_strain, 1.22 * reverse_t_strain));
	reload_t_stress = response(0);
	reload_t_stiffness = response(1);
	reverse_below_t_stress = std::min(.85 * reverse_t_stress, reload_t_stress + .9 * reload_t_stiffness * (reverse_t_strain - reload_t_strain));
	reverse_below_t_stiffness = reverse_below_t_stress / (reverse_t_strain - residual_t_strain);

	if(reverse_below_c_stiffness == 0.) update_compression_reverse(1.05 * peak_strain);

	residual_c_stiffness = std::min(residual_c_stiffness, .9 * reverse_below_t_stress / (reverse_t_strain - residual_c_strain));
	residual_t_stiffness = std::min(residual_t_stiffness, .9 * reverse_below_c_stress / (reverse_c_strain - residual_t_strain));
}

void ConcreteCM::update_compression_new_reverse(const double c_strain, const double c_stress) {
	const auto& reverse_c_strain = trial_history(2);     // unloading point strain compression
	const auto& reverse_c_stress = trial_history(3);     // unloading point stress compression
	const auto& residual_c_strain = trial_history(6);    // residual strain in compression unloading path
	const auto& reload_c_strain = trial_history(10);     // reload point strain compression
	const auto& reload_c_stress = trial_history(11);     // reload point stress compression
	const auto& reload_c_stiffness = trial_history(12);  // reload point stiffness compression
	auto& reverse_below_c_stress = trial_history(16);    // unloading point stress compression
	auto& reverse_below_c_stiffness = trial_history(17); // unloading point stiffness compression

	const auto diff_strain = reverse_c_strain - c_strain;
	const auto ratio = diff_strain / (reverse_c_strain - residual_c_strain);
	reverse_below_c_stress = std::max(reverse_c_stress - (reverse_c_stress - reverse_below_c_stress) * ratio, reload_c_stress + .9 * reload_c_stiffness * (reverse_c_strain - reload_c_strain));

	reverse_below_c_stiffness = (reverse_below_c_stress - c_stress) / diff_strain;
}

void ConcreteCM::update_tension_new_reverse(const double c_strain, const double c_stress) {
	const auto& reverse_t_strain = trial_history(4);     // unloading point strain tension
	const auto& reverse_t_stress = trial_history(5);     // unloading point stress tension
	const auto& residual_t_strain = trial_history(8);    // residual strain in compression unloading path
	const auto& reload_t_strain = trial_history(13);     // reload point strain tension
	const auto& reload_t_stress = trial_history(14);     // reload point stress tension
	const auto& reload_t_stiffness = trial_history(15);  // reload point stiffness tension
	auto& reverse_below_t_stress = trial_history(18);    // unloading point stress tension
	auto& reverse_below_t_stiffness = trial_history(19); // unloading point stiffness tension

	const auto diff_strain = reverse_t_strain - c_strain;
	const auto ratio = diff_strain / (reverse_t_strain - residual_t_strain);
	reverse_below_t_stress = std::min(reverse_t_stress - (reverse_t_stress - reverse_below_t_stress) * ratio, reload_t_stress + .9 * reload_t_stiffness * (reverse_t_strain - reload_t_strain));
	reverse_below_t_stiffness = (reverse_below_t_stress - c_stress) / diff_strain;
}

podarray<double> ConcreteCM::compute_transition(const double TX, const double XS, const double YS, const double ES, const double XF, const double YF, const double EF, const bool LT) {
	podarray<double> response(2);

	if(TX == XS) {
		response(0) = YS;
		response(1) = ES;
		return response;
	}

	if(TX == XF) {
		response(0) = YF;
		response(1) = EF;
		return response;
	}

	if(LT) return compute_linear_transition(TX, XS, YS, XF, YF);

	const auto TA = XF - XS;
	const auto TC = TX - XS;
	const auto ESEC = (YF - YS) / TA; // eq. 3-119
	const auto TB = ESEC - ES;
	const auto R = (EF - ESEC) / TB; // eq. 3-122
	const auto TD = TB * pow(std::max(fabs(TC / TA), datum::eps), R);

	response(0) = YS + TC * (ES + TD); // eq. 3-120
	response(1) = ES + (R + 1.) * TD;  // eq. 3-121

	suanpan_debug([&]() { if(!std::isfinite(response(0)) || !std::isfinite(response(1))) suanpan_fatal("infinite numbers detected.\n"); });

	return response;
}

podarray<double> ConcreteCM::compute_linear_transition(const double TX, const double XS, const double YS, const double XF, const double YF) {
	podarray<double> response(2);

	response(0) = (response(1) = (YF - YS) / (XF - XS)) * (TX - XS) + YS;

	return response;
}
