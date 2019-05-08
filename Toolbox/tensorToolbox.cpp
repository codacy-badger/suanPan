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

#include "tensorToolbox.h"

mat tensor::isotropic_stiffness(const double modulus, const double poissons_ratio) {
	const auto shear_modulus = modulus / (2. + 2. * poissons_ratio);
	const auto lambda = shear_modulus * poissons_ratio / (.5 - poissons_ratio);

	mat stiffness(6, 6, fill::zeros);

	for(auto I = 0; I < 3; ++I) for(auto J = 0; J < 3; ++J) stiffness(I, J) = lambda;
	for(auto I = 0; I < 3; ++I) stiffness(I, I) += 2. * shear_modulus;

	stiffness(3, 3) = stiffness(4, 4) = stiffness(5, 5) = shear_modulus;

	return stiffness;
}

mat tensor::isotropic_flexibility(const double modulus, const double poissons_ratio) {
	mat flexibility(6, 6, fill::zeros);

	flexibility(0, 1) = flexibility(1, 2) = flexibility(2, 0) = flexibility(0, 2) = flexibility(1, 0) = flexibility(2, 1) = -poissons_ratio / modulus;

	flexibility(0, 0) = flexibility(1, 1) = flexibility(2, 2) = 1. / modulus;

	flexibility(3, 3) = flexibility(4, 4) = flexibility(5, 5) = (2. * poissons_ratio + 2.) / modulus;

	return flexibility;
}

mat tensor::orthotropic_stiffness(const vec& modulus, const vec& poissons_ratio) {
	mat t_mat(3, 3);
	t_mat(0, 0) = 1. / modulus(0);
	t_mat(1, 0) = t_mat(0, 1) = -poissons_ratio(0) * (t_mat(1, 1) = 1. / modulus(1));
	t_mat(2, 1) = t_mat(1, 2) = -poissons_ratio(1) * (t_mat(2, 2) = 1. / modulus(2));
	t_mat(0, 2) = t_mat(2, 0) = -poissons_ratio(2) * t_mat(2, 2);

	mat stiffness(6, 6, fill::zeros);
	stiffness(span(0, 2), span(0, 2)) = inv(t_mat);
	stiffness(3, 3) = modulus(3);
	stiffness(4, 4) = modulus(4);
	stiffness(5, 5) = modulus(5);

	return stiffness;
}

mat tensor::orthotropic_flexibility(const vec&, const vec&) { throw invalid_argument("to be implemented"); }

mat tensor::unit_deviatoric_tensor4() {
	mat T = zeros(6, 6);

	T(3, 3) = T(4, 4) = T(5, 5) = .5;

	for(auto I = 0; I < 3; ++I) for(auto J = 0; J < 3; ++J) T(I, J) = I == J ? 2. / 3. : -1. / 3.;

	return T;
}

mat tensor::unit_deviatoric_tensor4v2() {
	mat T = eye(6, 6);

	T(span(0, 2), span(0, 2)) -= 1. / 3.;

	return T;
}

mat tensor::unit_symmetric_tensor4() {
	mat T = zeros(6, 6);

	for(auto I = 0; I < 3; ++I) T(I, I) = 1.;
	for(auto I = 3; I < 6; ++I) T(I, I) = .5;

	return T;
}

/**
 * \brief compute the first invariant of the given 3D strain tensor, could be either normal or deviatoric strain
 * \param E 3D strain tensor in Voigt notation
 * \return the first invariant trace(E)
 */
double tensor::strain::invariant1(const vec& E) {
	if(E.n_elem == 3 || E.n_elem == 6) return E(0) + E(1) + E(2);

	throw invalid_argument("need a valid strain vector");
}

/**
 * \brief compute the second invariant of the given 3D strain tensor, could be either normal or deviatoric strain
 * \param E 3D strain tensor in Voigt notation
 * \return the second invariant 0.5*(trace(E^2)-trace(E)^2)
 */
double tensor::strain::invariant2(const vec& E) {
	if(E.n_elem == 3) return -E(0) * E(1) - E(1) * E(2) - E(2) * E(0);
	if(E.n_elem == 6) return -E(0) * E(1) - E(1) * E(2) - E(2) * E(0) + .25 * (E(3) * E(3) + E(4) * E(4) + E(5) * E(5));

	throw invalid_argument("need a valid strain vector");
}

/**
 * \brief compute the third invariant of the given 3D strain tensor, could be either normal or deviatoric strain
 * \param E 3D strain tensor in Voigt notation
 * \return the third invariant det(E)
 */
double tensor::strain::invariant3(const vec& E) {
	if(E.n_elem == 3) return prod(E);
	if(E.n_elem == 6) return E(0) * E(1) * E(2) + .25 * (E(3) * E(4) * E(5) - E(0) * E(4) * E(4) - E(1) * E(5) * E(5) - E(2) * E(3) * E(3));

	throw invalid_argument("need a valid strain vector");
}

/**
 * \brief compute the first invariant of the given 3D stress tensor, could be either normal or deviatoric stress
 * \param S 3D stress tensor in Voigt notation
 * \return the first invariant trace(S)
 */
double tensor::stress::invariant1(const vec& S) {
	if(S.n_elem == 3 || S.n_elem == 6) return S(0) + S(1) + S(2);

	throw invalid_argument("need a valid stress vector");
}

/**
 * \brief compute the second invariant of the given 3D stress tensor, could be either normal or deviatoric stress
 * \param S 3D stress tensor in Voigt notation
 * \return the second invariant 0.5*(trace(S^2)-trace(S)^2)
 */
double tensor::stress::invariant2(const vec& S) {
	if(S.n_elem == 3) return -S(0) * S(1) - S(1) * S(2) - S(2) * S(0);
	if(S.n_elem == 6) return -S(0) * S(1) - S(1) * S(2) - S(2) * S(0) + S(3) * S(3) + S(4) * S(4) + S(5) * S(5);

	throw invalid_argument("need a valid stress vector");
}

/**
 * \brief compute the third invariant of the given 3D stress tensor, could be either normal or deviatoric stress
 * \param S 3D stress tensor in Voigt notation
 * \return the third invariant det(S)
 */
double tensor::stress::invariant3(const vec& S) {
	if(S.n_elem == 3) return prod(S);
	if(S.n_elem == 6) return S(0) * S(1) * S(2) + 2. * S(3) * S(4) * S(5) - S(0) * S(4) * S(4) - S(1) * S(5) * S(5) - S(2) * S(3) * S(3);

	throw invalid_argument("need a valid stress vector");
}

double tensor::strain::lode(const vec& E) {
	suanpan_debug([&]() { if(E.n_elem != 3) throw invalid_argument("need principal strain"); });
	const auto J2 = invariant2(E);
	if(J2 <= datum::eps) return 0.;
	auto value = sqrt(6.75) * invariant3(E) * pow(J2, -1.5);
	if(value < -1.) value = -1.;
	if(value > 1.) value = 1.;
	return acos(value) / 3.;
}

double tensor::stress::lode(const vec& S) {
	suanpan_debug([&]() { if(S.n_elem != 3) throw invalid_argument("need principal stress"); });
	const auto J2 = invariant2(S);
	if(J2 <= datum::eps) return 0.;
	auto value = sqrt(6.75) * invariant3(S) * pow(std::max(datum::eps, J2), -1.5);
	if(value < -1.) value = -1.;
	if(value > 1.) value = 1.;
	return acos(value) / 3.;
}

double tensor::trace(const vec& S) { return S(0) + S(1) + S(2); }

double tensor::mean(const vec& S) { return trace(S) / 3.; }

vec tensor::dev(const vec& S) {
	auto D = S;
	D(span(0, 2)) -= mean(D);
	return D;
}

mat tensor::dev(const mat& in) {
	auto out = in;
	out.diag() -= mean(out.diag());
	return out;
}

// transform deformation gradient to green strain
mat tensor::strain::to_green(mat&& gradient) {
	if(gradient.n_elem == 9) {
		const mat t_gradient(gradient.memptr(), 3, 3);
		return .5 * (t_gradient.t() * t_gradient - eye(size(t_gradient)));
	}
	if(gradient.n_elem == 4) {
		const mat t_gradient(gradient.memptr(), 2, 2);
		return .5 * (t_gradient.t() * t_gradient - eye(size(t_gradient)));
	}
	throw invalid_argument("need a valid strain vector");
}

mat tensor::strain::to_green(const mat& gradient) {
	if(gradient.n_elem == 9) {
		const mat t_gradient(gradient.memptr(), 3, 3);
		return .5 * (t_gradient.t() * t_gradient - eye(size(t_gradient)));
	}
	if(gradient.n_elem == 4) {
		const mat t_gradient(gradient.memptr(), 2, 2);
		return .5 * (t_gradient.t() * t_gradient - eye(size(t_gradient)));
	}
	throw invalid_argument("need a valid strain vector");
}

mat tensor::strain::to_tensor(const vec& in_strain) {
	mat out_strain;

	if(in_strain.n_elem == 3) {
		out_strain.set_size(2, 2);
		out_strain(0, 0) = in_strain(0);
		out_strain(1, 1) = in_strain(1);
		out_strain(0, 1) = out_strain(1, 0) = .5 * in_strain(2);
		return out_strain;
	}

	if(in_strain.n_elem == 6) {
		out_strain.set_size(3, 3);
		out_strain(0, 0) = in_strain(0);
		out_strain(1, 1) = in_strain(1);
		out_strain(2, 2) = in_strain(2);
		out_strain(0, 1) = out_strain(1, 0) = .5 * in_strain(3);
		out_strain(1, 2) = out_strain(2, 1) = .5 * in_strain(4);
		out_strain(2, 0) = out_strain(0, 2) = .5 * in_strain(5);
		return out_strain;
	}

	throw invalid_argument("need a valid stress vector");
}

vec tensor::strain::to_voigt(const mat& in_strain) {
	vec out_strain;

	if(in_strain.n_elem == 9) {
		out_strain.set_size(6);

		out_strain(0) = in_strain(0);
		out_strain(1) = in_strain(4);
		out_strain(2) = in_strain(8);
		out_strain(3) = in_strain(1) + in_strain(3);
		out_strain(4) = in_strain(5) + in_strain(7);
		out_strain(5) = in_strain(2) + in_strain(6);

		return out_strain;
	}

	if(in_strain.n_elem == 4) {
		out_strain.set_size(3);

		out_strain(0) = in_strain(0);
		out_strain(1) = in_strain(3);
		out_strain(2) = in_strain(1) + in_strain(2);

		return out_strain;
	}

	throw invalid_argument("need a valid strain tensor");
}

double tensor::strain::norm(const vec& in) {
	if(in.n_elem == 6) return sqrt(dot(norm_weight, square(in)));
	if(in.n_elem == 3) return arma::norm(in);
	throw invalid_argument("need a valid strain vector");
}

mat tensor::stress::to_tensor(const vec& in_stress) {
	mat out_stress;

	if(in_stress.n_elem == 3) {
		out_stress.set_size(2, 2);
		out_stress(0, 0) = in_stress(0);
		out_stress(1, 1) = in_stress(1);
		out_stress(0, 1) = out_stress(1, 0) = in_stress(2);
		return out_stress;
	}

	if(in_stress.n_elem == 6) {
		out_stress.set_size(3, 3);
		out_stress(0, 0) = in_stress(0);
		out_stress(1, 1) = in_stress(1);
		out_stress(2, 2) = in_stress(2);
		out_stress(0, 1) = out_stress(1, 0) = in_stress(3);
		out_stress(1, 2) = out_stress(2, 1) = in_stress(4);
		out_stress(2, 0) = out_stress(0, 2) = in_stress(5);
		return out_stress;
	}

	throw invalid_argument("need a valid stress vector");
}

vec tensor::stress::to_voigt(const mat& in_stress) {
	vec out_stress;

	if(in_stress.n_elem == 9) {
		out_stress.set_size(6);

		out_stress(0) = in_stress(0);
		out_stress(1) = in_stress(4);
		out_stress(2) = in_stress(8);
		out_stress(3) = in_stress(1);
		out_stress(4) = in_stress(2);
		out_stress(5) = in_stress(5);

		return out_stress;
	}

	if(in_stress.n_elem == 4) {
		out_stress.set_size(3);

		out_stress(0) = in_stress(0);
		out_stress(1) = in_stress(3);
		out_stress(2) = in_stress(1);

		return out_stress;
	}

	throw invalid_argument("need a valid stress tensor");
}

double tensor::stress::norm(const vec& in) {
	if(in.n_elem == 6) return sqrt(dot(norm_weight, square(in)));
	if(in.n_elem == 3) return arma::norm(in);
	throw invalid_argument("need a valid stress vector");
}

double transform::atan2(const vec& direction_cosine) { return std::atan2(direction_cosine(1), direction_cosine(0)); }

mat transform::compute_jacobian_nominal_to_principal(const mat& in) {
	mat out(3, 6);

	out(span(0, 2), span(0, 2)) = square(in).t();

	out(0, 3) = 2. * in(0, 0) * in(1, 0);
	out(0, 4) = 2. * in(1, 0) * in(2, 0);
	out(0, 5) = 2. * in(2, 0) * in(0, 0);

	out(1, 3) = 2. * in(0, 1) * in(1, 1);
	out(1, 4) = 2. * in(1, 1) * in(2, 1);
	out(1, 5) = 2. * in(2, 1) * in(0, 1);

	out(2, 3) = 2. * in(0, 2) * in(1, 2);
	out(2, 4) = 2. * in(1, 2) * in(2, 2);
	out(2, 5) = 2. * in(2, 2) * in(0, 2);

	return out;
}

mat transform::compute_jacobian_principal_to_nominal(const mat& in) {
	mat out(6, 3);

	out(span(0, 2), span(0, 2)) = square(in);

	out(3, 0) = in(0, 0) * in(1, 0);
	out(4, 0) = in(1, 0) * in(2, 0);
	out(5, 0) = in(2, 0) * in(0, 0);

	out(3, 1) = in(0, 1) * in(1, 1);
	out(4, 1) = in(1, 1) * in(2, 1);
	out(5, 1) = in(2, 1) * in(0, 1);

	out(3, 2) = in(0, 2) * in(1, 2);
	out(4, 2) = in(1, 2) * in(2, 2);
	out(5, 2) = in(2, 2) * in(0, 2);

	return out;
}

mat transform::compute_differential_nominal_to_principal(const mat& in) {
	const auto I1 = in(0) + in(1) + in(2);
	const auto I2 = in(0) * in(1) + in(1) * in(2) + in(2) * in(0) - in(3) * in(3) - in(4) * in(4) - in(5) * in(5);
	const auto I3 = in(0) * in(1) * in(2) + 2. * in(3) * in(4) * in(5) - in(0) * in(4) * in(4) - in(1) * in(5) * in(5) - in(2) * in(3) * in(3);

	const rowvec di1ds{1., 1., 1., 0., 0., 0.};
	const rowvec di2ds{in(1) + in(2), in(0) + in(2), in(1) + in(0), -2. * in(3), -2. * in(4), -2. * in(5)};
	const rowvec di3ds{in(1) * in(2) - in(4) * in(4), in(0) * in(2) - in(5) * in(5), in(0) * in(1) - in(3) * in(3), 2. * in(4) * in(5) - 2. * in(2) * in(3), 2. * in(3) * in(5) - 2. * in(0) * in(4), 2. * in(3) * in(4) - 2. * in(1) * in(5)};
	// const rowvec di2ds{in(1) + in(2), in(0) + in(2), in(1) + in(0), -in(3), -in(4), -in(5)};
	// const rowvec di3ds{in(1) * in(2) - in(4) * in(4), in(0) * in(2) - in(5) * in(5), in(0) * in(1) - in(3) * in(3), in(4) * in(5) - in(2) * in(3), in(3) * in(5) - in(0) * in(4), in(3) * in(4) - in(1) * in(5)};

	const auto Q = (I1 * I1 - 3. * I2) / 9.;
	const auto R = I1 * I2 / 6. - pow(I1 / 3., 3.) - .5 * I3;

	const rowvec dqds = 2. / 9. * I1 * di1ds - di2ds / 3.;
	const rowvec drds = (I2 / 6. - I1 * I1 / 9.) * di1ds + I1 / 6. * di2ds - .5 * di3ds;

	const auto part = pow(Q, -1.5);
	const auto value = std::max(-1., std::min(1., R * part));

	// 2./3.*Q*dthetads
	const rowvec dthetads = 2. / 3. * pow(1. - value * value + datum::eps, -.5) * part * (1.5 * R * dqds - Q * drds);

	// theta/3.
	const auto theta = std::acos(value) / 3.;
	const auto twp_third_pi = 2. / 3. * datum::pi;

	const vec cos_theta{cos(theta), cos(theta - twp_third_pi), cos(theta + twp_third_pi)};
	const vec sin_theta{sin(theta), sin(theta - twp_third_pi), sin(theta + twp_third_pi)};

	mat out = (sin_theta * dthetads - cos_theta * dqds) / sqrt(Q);

	out.cols(0, 2) += 1. / 3.;

	suanpan_debug([&]() { if(!out.is_finite()) throw logic_error("nan detected"); });

	return out;
}

mat transform::compute_differential_principal_to_nominal(const mat&) {
	mat out(6, 3);

	return out;
}

vec transform::haigh_westergaard_to_principal(const double sigma, const double rho, const double theta) {
	const auto two_over_three_pi = 2. / 3. * datum::pi;
	const auto sqrt_three_over_two = sqrt(1.5);

	vec p_stress(3);

	p_stress(2) = cos(theta);
	p_stress(1) = cos(theta - two_over_three_pi);
	p_stress(0) = cos(theta + two_over_three_pi);
	p_stress *= rho / sqrt_three_over_two;
	p_stress += sigma;

	return p_stress;
}

vec transform::triangle::to_area_coordinate(const vec& g_coord, const mat& nodes) {
	if(nodes.n_cols != 2 || nodes.n_rows != 3) throw invalid_argument("need 3 by 2 mat");

	vec a_coord(3);
	a_coord(0) = 1.;
	a_coord(1) = g_coord(0);
	a_coord(2) = g_coord(1);

	mat element_coor(3, 3);
	element_coor.row(0).fill(1.);
	element_coor.rows(1, 2) = nodes.t();

	return solve(element_coor, a_coord);
}

double transform::strain::angle(const vec& strain) { return .5 * std::atan2(strain(2), strain(0) - strain(1)); }

mat transform::strain::trans(const double angle) {
	const auto sin_angle = sin(2. * angle);
	const auto cos_angle = cos(2. * angle);

	mat trans(3, 3);
	trans(0, 0) = trans(1, 1) = .5 + .5 * cos_angle;
	trans(0, 1) = trans(1, 0) = .5 - .5 * cos_angle;
	trans(1, 2) = -(trans(0, 2) = .5 * sin_angle);
	trans(2, 0) = -(trans(2, 1) = sin_angle);
	trans(2, 2) = cos_angle;

	return trans;
}

vec transform::strain::principal(const vec& strain) {
	const auto tmp_a = .5 * (strain(0) + strain(1));
	const auto tmp_b = .5 * sqrt(pow(strain(0) - strain(1), 2.) + pow(strain(2), 2.));

	vec p_strain(3);
	p_strain(0) = tmp_a + tmp_b;
	p_strain(1) = tmp_a - tmp_b;
	p_strain(2) = 0.;

	return p_strain;
}

vec transform::strain::rotate(const vec& strain, const double theta) { return trans(theta) * strain; }

double transform::stress::angle(const vec& stress) { return .5 * std::atan2(2. * stress(2), stress(0) - stress(1)); }

mat transform::stress::trans(const double angle) {
	const auto sin_angle = sin(2. * angle);
	const auto cos_angle = cos(2. * angle);

	mat trans(3, 3);
	trans(0, 0) = trans(1, 1) = .5 + .5 * cos_angle;
	trans(0, 1) = trans(1, 0) = .5 - .5 * cos_angle;
	trans(1, 2) = -(trans(0, 2) = sin_angle);
	trans(2, 0) = -(trans(2, 1) = .5 * sin_angle);
	trans(2, 2) = cos_angle;

	return trans;
}

vec transform::stress::principal(const vec& stress) {
	const auto tmp_a = .5 * (stress(0) + stress(1));
	const auto tmp_b = .5 * sqrt(pow(stress(0) - stress(1), 2.) + pow(2. * stress(2), 2.));

	vec p_stress(3);
	p_stress(0) = tmp_a + tmp_b;
	p_stress(1) = tmp_a - tmp_b;
	p_stress(2) = 0.;

	return p_stress;
}

vec transform::stress::rotate(const vec& stress, const double theta) { return trans(theta) * stress; }

mat transform::beam::global_to_local(const double cos, const double sin, const double length) {
	mat trans_mat(3, 6, fill::zeros);

	trans_mat(0, 0) = -(trans_mat(0, 3) = cos);
	trans_mat(0, 1) = -(trans_mat(0, 4) = sin);
	trans_mat(1, 0) = trans_mat(2, 0) = -(trans_mat(1, 3) = trans_mat(2, 3) = sin / length);
	trans_mat(1, 4) = trans_mat(2, 4) = -(trans_mat(1, 1) = trans_mat(2, 1) = cos / length);
	trans_mat(1, 2) = trans_mat(2, 5) = 1.;

	return trans_mat;
}

/**
 * \brief a subroutine to get transformation matrix between global and local coordinate system
 * \param direction_cosine the direction cosine consists of either two or three elements
 * \param length the length of the member
 * \return a matrix that transforms global quantities to local quantities with rigid body motions removed
 */
mat transform::beam::global_to_local(const vec& direction_cosine, const double length) {
	if(direction_cosine.n_elem == 2) return global_to_local(direction_cosine(0), direction_cosine(1), length);
	throw logic_error("direction cosine must contains two or three elements.\n");
}
