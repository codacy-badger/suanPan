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

#ifndef TENSORTOOLBOX_H
#define TENSORTOOLBOX_H

#include <suanPan.h>

namespace tensor {
	mat isotropic_stiffness(double, double);
	mat isotropic_flexibility(double, double);
	mat orthotropic_stiffness(const vec&, const vec&);
	mat orthotropic_flexibility(const vec&, const vec&);

	mat unit_deviatoric_tensor4();
	mat unit_deviatoric_tensor4v2();
	mat unit_symmetric_tensor4();

	static const vec unit_tensor2{1., 1., 1., 0., 0., 0.};

	namespace stress {
		// applies to 3D tensor only, either principal or not
		double invariant1(const vec&);
		// applies to 3D tensor only, either principal or not
		double invariant2(const vec&);
		// applies to 3D tensor only, either principal or not
		double invariant3(const vec&);

		double lode(const vec&);

		static const vec norm_weight{1., 1., 1., 2., 2., 2.};
	} // namespace stress
	namespace strain {
		// applies to 3D tensor only, either principal or not
		double invariant1(const vec&);
		// applies to 3D tensor only, either principal or not
		double invariant2(const vec&);
		// applies to 3D tensor only, either principal or not
		double invariant3(const vec&);

		double lode(const vec&);

		static const vec norm_weight{1., 1., 1., .5, .5, .5};
	} // namespace strain
	double trace(const vec&);
	double mean(const vec&);
	vec dev(const vec&);

	mat dev(const mat&);

	namespace strain {
		mat to_green(mat&&);
		mat to_green(const mat&);
		mat to_tensor(const vec&);
		vec to_voigt(const mat&);
		double norm(const vec&);
	} // namespace strain
	namespace stress {
		mat to_tensor(const vec&);
		vec to_voigt(const mat&);
		double norm(const vec&);
	} // namespace stress
}     // namespace tensor

namespace transform {
	double atan2(const vec&);
	mat compute_jacobian_nominal_to_principal(const mat&);
	mat compute_jacobian_principal_to_nominal(const mat&);
	mat compute_differential_nominal_to_principal(const mat&);
	mat compute_differential_principal_to_nominal(const mat&);
	vec haigh_westergaard_to_principal(double, double, double);

	namespace strain {
		double angle(const vec&);
		mat trans(double);
		vec principal(const vec&);
		vec rotate(const vec&, double);
	} // namespace strain
	namespace stress {
		double angle(const vec&);
		mat trans(double);
		vec principal(const vec&);
		vec rotate(const vec&, double);
	} // namespace stress
	namespace beam {
		mat global_to_local(double, double, double);
		mat global_to_local(const vec&, double);
	} // namespace beam
	namespace triangle {
		vec to_area_coordinate(const vec&, const mat&);
	}
} // namespace transform

namespace suanpan {
	template<typename T> T ramp(const T in) { return in > T(0) ? in : T(0); }
} // namespace suanpan

#endif
