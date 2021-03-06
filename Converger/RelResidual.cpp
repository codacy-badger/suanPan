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

#include "RelResidual.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

RelResidual::RelResidual(const unsigned T, const double E, const unsigned M, const bool P)
	: Converger(T, E, M, P) {}

bool RelResidual::is_converged() {
	const auto& D = get_domain().lock();
	auto& W = D->get_factory();

	vec residual = W->get_trial_load() - W->get_trial_resistance();

	if(!W->get_reference_load().is_empty() && !W->get_trial_load_factor().is_empty()) residual += W->get_reference_load() * W->get_trial_load_factor();

	for(const auto& t_dof : D->get_restrained_dof()) residual(t_dof) = 0.;

	if(-1. == ref_residual) {
		ref_residual = norm(residual);
		set_error(1.);
	} else set_error(norm(residual) / ref_residual);

	if(get_tolerance() > get_error()) {
		set_conv_flag(true);
		ref_residual = -1.;
	} else set_conv_flag(false);

	if(is_print()) suanpan_info("relative residual: %.5E.\n", get_error());

	return get_conv_flag();
}
