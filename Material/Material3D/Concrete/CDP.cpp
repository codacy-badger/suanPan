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

#include "CDP.h"
#include <Recorder/OutputType.h>
#include <Toolbox/tensorToolbox.h>

const unsigned CDP::max_iteration = 20;
const double CDP::root_three_two = sqrt(1.5);
const mat CDP::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

double CDP::compute_r(const vec& in) {
	const auto r = .5 + .5 * accu(in) / accu(abs(in));
	return !std::isfinite(r) || r < 0. ? .0 : r > 1. ? 1. : r;
}

vec CDP::compute_dr(const vec& in) {
	const auto g = accu(abs(in));

	vec out = .5 * (vec(3).fill(g) - accu(in) * sign(in)) / (g * g);

	if(!out.is_finite()) out.zeros();

	return out;
}

double CDP::compute_s(const double r) const { return s0 + r - s0 * r; }

podarray<double> CDP::compute_tension_backbone(const double kappa) {
	podarray<double> out(6);

	const auto s_phi = sqrt(1. + a_t * (a_t + 2.) * kappa);
	const auto t_phi = (1. + .5 * a_t) / s_phi;
	const auto b_phi = (1. + a_t - s_phi) / a_t;
	const auto p_phi = pow(b_phi, cb_t);

	out(0) = 1. - p_phi;                                    // d
	out(1) = f_t * s_phi * b_phi;                           // f
	out(2) = out(1) / p_phi;                                // \bar{f}
	out(3) = cb_t * t_phi * p_phi / b_phi;                  // \md{d}
	out(4) = f_t * t_phi * (1. + a_t - 2. * s_phi);         // \md{f}
	out(5) = (out(4) + f_t * t_phi * cb_t * s_phi) / p_phi; // \md{\bar{f}}

	return out;
}

podarray<double> CDP::compute_compression_backbone(const double kappa) {
	podarray<double> out(6);

	const auto s_phi = sqrt(1. + a_c * (a_c + 2.) * kappa);
	const auto t_phi = (1. + .5 * a_c) / s_phi;
	const auto b_phi = (1. + a_c - s_phi) / a_c;
	const auto p_phi = pow(b_phi, cb_c);

	out(0) = 1. - p_phi;                                    // d
	out(1) = f_c * s_phi * b_phi;                           // f
	out(2) = out(1) / p_phi;                                // \bar{f}
	out(3) = cb_c * t_phi * p_phi / b_phi;                  // \md{d}
	out(4) = f_c * t_phi * (1. + a_c - 2. * s_phi);         // \md{f}
	out(5) = (out(4) + f_c * t_phi * cb_c * s_phi) / p_phi; // \md{\bar{f}}

	return out;
}

CDP::CDP(const unsigned T, const double E, const double V, const double ST, const double SC, const double GT, const double GC, const double AT, const double AC, const double DT, const double DC, const double AP, const double BC, const double S, const double R)
	: Material3D(T, R)
	, elastic_modulus(fabs(E))
	, poissons_ratio(V < .5 ? V : .2)
	, a_t(AT < 1. ? AT : .5)
	, cb_t(log(DT < 1. ? DT : .95) / log(.5 * (1. + a_t - sqrt(1. + a_t * a_t)) / a_t))
	, g_t(fabs(GT))
	, f_t(fabs(ST))
	, a_c(AC > 1. ? AC : 4.)
	, cb_c(log(DC < 1. ? DC : .95) / log(.5 + .5 / a_c))
	, g_c(fabs(GC))
	, f_c(-fabs(SC) * 4. * a_c * pow(1. + a_c, -2.))
	, alpha((fabs(BC) - 1.) / (2. * fabs(BC) - 1.))
	, alpha_p(fabs(AP))
	, s0(fabs(S)) {
	// tension
	const auto half_stress = .5 * f_t;
	const auto half_strain = log(1. + a_t + sqrt(1. + a_t * a_t)) / f_t / (1. + .5 * a_t) * g_t + half_stress / elastic_modulus;
	const auto ratio_t = half_stress / half_strain / elastic_modulus;
	if(ratio_t >= DT) {
		suanpan_warning("CDP model requires a minimum tension degradation of %.2f, now reset it.\n", ratio_t);
		access::rw(cb_t) = log(ratio_t >= .9 ? ratio_t : ratio_t + .05) / log(.5 * (1. + a_t - sqrt(1. + a_t * a_t)) / a_t);
	}
	// compression
	const auto peak_stress = -fabs(SC);
	const auto peak_strain = log(2. * a_c / (1. + a_c)) / f_c / (1. + .5 * a_c) * g_c + peak_stress / elastic_modulus;
	const auto ratio_c = peak_stress / peak_strain / elastic_modulus;
	if(ratio_c >= DC) {
		suanpan_warning("CDP model requires a minimum compression degradation of %.2f, now reset it.\n", ratio_c);
		access::rw(cb_c) = log(ratio_c >= .9 ? ratio_c : ratio_c + .05) / log(.5 + .5 / a_c);
	}

	access::rw(tolerance) = 1E-13;
}

void CDP::initialize(const shared_ptr<DomainBase>&) {
	access::rw(inv_stiffness) = tensor::isotropic_flexibility(elastic_modulus, poissons_ratio);

	trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

	trial_history = current_history.zeros(10);
}

unique_ptr<Material> CDP::get_copy() { return make_unique<CDP>(*this); }

double CDP::get_parameter(const ParameterType P) const {
	switch(P) {
	case ParameterType::DENSITY:
		return density;
	case ParameterType::ELASTICMODULUS:
	case ParameterType::YOUNGSMODULUS:
		return elastic_modulus;
	case ParameterType::POISSONSRATIO:
		return poissons_ratio;
	default:
		return 0.;
	}
}

int CDP::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

	trial_history = current_history;
	auto& d_t = trial_history(0);
	auto& d_c = trial_history(1);
	auto& kappa_t = trial_history(2);
	auto& kappa_c = trial_history(3);
	vec plastic_strain(&trial_history(4), 6, false, true);

	const auto& current_kappa_t = current_history(2);
	const auto& current_kappa_c = current_history(3);

	trial_stress = (trial_stiffness = initial_stiffness) * (trial_strain - plastic_strain); // 6

	vec principal_stress;    // 3
	mat principal_direction; // 3x3
	if(!eig_sym(principal_stress, principal_direction, tensor::stress::to_tensor(trial_stress), "std")) return SUANPAN_FAIL;

	const auto s = tensor::dev(trial_stress);    // 6
	const auto norm_s = tensor::stress::norm(s); // 1
	vec n = s / norm_s;                          // 6
	if(!n.is_finite()) n.zeros();

	const auto ps = tensor::dev(principal_stress); // 3
	const vec pn = normalise(ps);                  // 3

	const vec dsigmadlambda = -double_shear * pn - three_alpha_p_bulk; // 6

	const auto dgdsigma_t = (pn(2) + alpha_p) / g_t;
	const auto dgdsigma_c = (pn(0) + alpha_p) / g_c;

	auto new_stress = principal_stress;     // converged principal stress
	const auto& max_stress = new_stress(2); // algebraically maximum principal stress

	const auto const_yield = alpha * accu(principal_stress) + root_three_two * norm_s;

	vec residual(3), incre;
	mat jacobian(3, 3, fill::zeros);
	mat left(3, 6);

	podarray<double> t_para, c_para;

	auto lambda = 0., ref_error = 0.;
	double r, beta;
	vec dr;

	unsigned counter = 0;
	while(true) {
		if(max_iteration == ++counter) {
			suanpan_debug("CDP cannot converge within %u iterations.\n", max_iteration);
			return SUANPAN_FAIL;
		}

		t_para = compute_tension_backbone(kappa_t);
		c_para = compute_compression_backbone(kappa_c);

		const auto tension_flag = max_stress > 0.;

		beta = -one_minus_alpha * c_para(2) / t_para(2) - alpha - 1.;

		residual(0) = const_yield + pfplambda * lambda + one_minus_alpha * c_para(2);

		if(tension_flag) residual(0) += beta * max_stress;

		r = compute_r(new_stress);

		if(1 == counter && residual(0) < 0.) {
			const auto damage = (1. - d_c) * (1. - compute_s(r) * d_t);
			trial_stress *= damage;
			trial_stiffness *= damage;
			return SUANPAN_SUCCESS;
		}

		const auto t_term = t_para(1) * dgdsigma_t;
		const auto c_term = c_para(1) * dgdsigma_c;

		residual(1) = r * t_term * lambda + current_kappa_t - kappa_t;
		residual(2) = (c_term - r * c_term) * lambda + current_kappa_c - kappa_c;

		if(tension_flag) {
			jacobian(0, 0) = pfplambda + beta * dsigmadlambda(2);
			const auto tmp_term = one_minus_alpha * max_stress / t_para(2);
			jacobian(0, 1) = tmp_term * c_para(2) / t_para(2) * t_para(5);
			jacobian(0, 2) = (one_minus_alpha - tmp_term) * c_para(5);
		} else {
			jacobian(0, 0) = pfplambda;
			jacobian(0, 1) = 0.;
			jacobian(0, 2) = one_minus_alpha * c_para(5);
		}

		const auto dlambda = r + lambda * dot(dr = compute_dr(new_stress), dsigmadlambda);
		jacobian(1, 0) = t_term * dlambda;
		jacobian(2, 0) = c_term - c_term * dlambda;
		jacobian(1, 1) = r * lambda * dgdsigma_t * t_para(4) - 1.;
		jacobian(2, 2) = (lambda - r * lambda) * dgdsigma_c * c_para(4) - 1.;

		if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

		auto error = norm(residual);
		if(1 == counter) ref_error = std::max(1., error);
		suanpan_debug("CDP local iteraton error: %.5E.\n", error /= ref_error);
		if(error <= tolerance && norm(incre) <= tolerance) break;

		lambda -= incre(0);
		kappa_t -= incre(1);
		kappa_c -= incre(2);
		new_stress -= dsigmadlambda * incre(0);

		if(kappa_t > 1.) kappa_t = .95; // avoid overshoot
		if(kappa_c > 1.) kappa_c = .95; // avoid overshoot
	}

	// update damage indices
	d_t = t_para(0), d_c = c_para(0);
	// update plastic strain
	plastic_strain += lambda * (n + unit_alpha_p);

	const auto recovery = compute_s(r);
	const auto damage_c = d_c - 1.;
	const auto damage_t = recovery * d_t - 1.;
	const auto damage = damage_c * damage_t;

	// update trial stress
	trial_stress = transform::compute_jacobian_principal_to_nominal(principal_direction) * new_stress;

	const mat dnde = double_shear / norm_s * (unit_dev_tensor - n * n.t());

	// \dfrac{\partial\bar{\sigma}}{\partial\varepsilon^{tr}}
	trial_stiffness -= double_shear * lambda * dnde;

	const auto trans = transform::compute_jacobian_nominal_to_principal(principal_direction);

	const rowvec drdsigma = dr.t() * trans;
	const rowvec prpe = drdsigma * trial_stiffness;

	// compute local derivatives
	left.row(0) = 3. * alpha * bulk * tensor::unit_tensor2.t() + root_three_two * double_shear * n.t();
	left.row(1) = t_para(1) * lambda * (r / g_t * trans.row(2) * dnde + dgdsigma_t * prpe);
	left.row(2) = c_para(1) * lambda * ((1. - r) / g_c * trans.row(0) * dnde - dgdsigma_c * prpe);

	if(max_stress > 0.) left.row(0) += beta * trans.row(2) * trial_stiffness;

	const mat right = -solve(jacobian, left);
	const auto& dlambdade = right.row(0);
	const auto& dkappade = right.rows(1, 2);

	// \dfrac{\mathrm{d}\bar{\sigma}}{\mathrm{d}\varepsilon^{tr}}
	trial_stiffness -= (double_shear * n + three_alpha_p_bulk * tensor::unit_tensor2) * dlambdade;

	trial_stiffness = (damage * eye(6, 6) + d_t * damage_c * (1. - s0) * trial_stress * drdsigma) * trial_stiffness + trial_stress * rowvec{recovery * damage_c * t_para(3), damage_t * c_para(3)} * dkappade;

	trial_stress *= damage;

	return SUANPAN_SUCCESS;
}

int CDP::clear_status() {
	current_strain.zeros();
	current_stress.zeros();
	current_history.zeros();
	current_stiffness = initial_stiffness;
	return reset_status();
}

int CDP::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_history = trial_history;
	current_stiffness = trial_stiffness;
	return SUANPAN_SUCCESS;
}

int CDP::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_history = current_history;
	trial_stiffness = current_stiffness;
	return SUANPAN_SUCCESS;
}

vector<vec> CDP::record(const OutputType T) {
	vector<vec> data;

	if(T == OutputType::DT) data.emplace_back(vec{current_history(0)});
	else if(T == OutputType::DC) data.emplace_back(vec{current_history(1)});
	else if(T == OutputType::KAPPAT) data.emplace_back(vec{current_history(2)});
	else if(T == OutputType::KAPPAC) data.emplace_back(vec{current_history(3)});
	else return Material3D::record(T);

	return data;
}

void CDP::print() { suanpan_info("A concrete damage plasticity model.\n"); }
