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

#include "PlaneStress.h"
#include <Domain/DomainBase.h>

const uvec PlaneStress::F1{0, 1, 3};
const uvec PlaneStress::F2{2, 4, 5};

mat PlaneStress::form_stiffness(const mat& full_stiffness) {
	const auto& A = full_stiffness(F1, F1);
	const auto& U = full_stiffness(F1, F2);
	const auto& C = full_stiffness(F2, F2);
	const auto& V = full_stiffness(F2, F1);

	return A - U * solve(C, V);

	/*
	const auto S = inv(A);
	const auto VS = V * S;

	return inv(S * (eye(3, 3) - U * solve(VS * U - C, VS)));
	*/
}

PlaneStress::PlaneStress(const unsigned T, const unsigned BT, const unsigned MI, const bool FM)
	: Material2D(T, PlaneType::S, 0.)
	, base_tag(BT)
	, max_iteration(MI)
	, use_full_matrix(FM) { access::rw(tolerance) = 1E2 * tolerance; }

PlaneStress::PlaneStress(const PlaneStress& old_obj)
	: Material2D(old_obj)
	, base_tag(old_obj.base_tag)
	, max_iteration(old_obj.max_iteration)
	, use_full_matrix(old_obj.use_full_matrix)
	, base(old_obj.base == nullptr ? nullptr : old_obj.base->get_copy())
	, trial_full_strain(old_obj.trial_full_strain)
	, current_full_strain(old_obj.current_full_strain) {}

void PlaneStress::initialize(const shared_ptr<DomainBase>& D) {
	if(D == nullptr) return;

	if(!D->find_material(base_tag)) {
		D->disable_material(get_tag());
		return;
	}

	base = get_material(D, base_tag)->get_copy();

	trial_full_strain = current_full_strain.zeros(6);

	current_stiffness = trial_stiffness = initial_stiffness = form_stiffness(base->get_initial_stiffness());
}

double PlaneStress::get_parameter(const ParameterType P) const { return base->get_parameter(P); }

unique_ptr<Material> PlaneStress::get_copy() { return make_unique<PlaneStress>(*this); }

int PlaneStress::update_trial_status(const vec& t_strain) {
	auto& t_stress = base->get_trial_stress();
	auto& t_stiffness = base->get_trial_stiffness();

	if(1 == max_iteration) {
		incre_strain = t_strain - trial_strain;

		trial_full_strain(F1) = trial_strain = t_strain;

		if(use_full_matrix) {
			trial_full_strain(F2) -= solve(t_stiffness(F2, F2), t_stress(F2) + t_stiffness(F2, F1) * incre_strain);

			if(base->update_trial_status(trial_full_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

			trial_stress = t_stress(F1) - t_stiffness(F1, F2) * solve(t_stiffness(F2, F2), t_stress(F2));
		} else {
			trial_full_strain(F2) -= (t_stress(F2) + t_stiffness(F2, F1) * incre_strain) / vec(t_stiffness.diag())(F2);

			if(base->update_trial_status(trial_full_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

			trial_stress = t_stress(F1) - t_stiffness(F1, F2) * t_stress(F2) / vec(t_stiffness.diag())(F2);
		}
	} else {
		trial_full_strain(F1) = trial_strain = t_strain;

		unsigned counter = 0;

		if(use_full_matrix)
			while(++counter < max_iteration) {
				if(base->update_trial_status(trial_full_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
				trial_full_strain(F2) -= solve(t_stiffness(F2, F2), t_stress(F2));
				const auto error = norm(t_stress(F2));
				suanpan_extra_debug("PlaneStress state determination error: %.4E.\n", error);
				if(error < tolerance) break;
			}
		else
			while(++counter < max_iteration) {
				if(base->update_trial_status(trial_full_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
				trial_full_strain(F2) -= t_stress(F2) / vec(t_stiffness.diag())(F2);
				const auto error = norm(t_stress(F2));
				suanpan_extra_debug("PlaneStress state determination error: %.4E.\n", error);
				if(error < tolerance) break;
			}

		if(counter == max_iteration) {
			suanpan_error("PlaneStress cannot converge with in %u iterations.\n", counter);
			return SUANPAN_FAIL;
		}

		suanpan_extra_debug("PlaneStress local iteration: %u.\n", counter);

		trial_stress = t_stress(F1);
	}

	trial_stiffness = form_stiffness(t_stiffness);

	return SUANPAN_SUCCESS;
}

int PlaneStress::clear_status() {
	current_full_strain.zeros(6);
	trial_full_strain.zeros(6);
	Material::clear_status();
	return base->clear_status();
}

int PlaneStress::commit_status() {
	current_full_strain = trial_full_strain;
	Material::commit_status();
	return base->commit_status();
}

int PlaneStress::reset_status() {
	trial_full_strain = current_full_strain;
	Material::reset_status();
	return base->reset_status();
}

vector<vec> PlaneStress::record(const OutputType P) { return base->record(P); }

void PlaneStress::print() {
	suanpan_info("A plane stress wrapper.\n");
	suanpan_info("current strain: ");
	current_strain.t().print();
	suanpan_info("current stress: ");
	current_stress.t().print();
	base->print();
}
