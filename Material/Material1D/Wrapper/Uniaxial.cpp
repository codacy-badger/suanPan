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

#include "Uniaxial.h"
#include <Domain/DomainBase.h>

const uvec Uniaxial::F1{0};
const uvec Uniaxial::F2{1, 2, 3, 4, 5};

mat Uniaxial::form_stiffness(const mat& full_stiffness) { return full_stiffness(F1, F1) - full_stiffness(F1, F2) * solve(full_stiffness(F2, F2), full_stiffness(F2, F1)); }

Uniaxial::Uniaxial(const unsigned T, const unsigned BT, const unsigned MI)
	: Material1D(T, 0.)
	, base_tag(BT)
	, max_iteration(MI) {}

Uniaxial::Uniaxial(const Uniaxial& old_obj)
	: Material1D(old_obj)
	, base_tag(old_obj.base_tag)
	, max_iteration(old_obj.max_iteration)
	, base(old_obj.base->get_copy())
	, trial_full_strain(old_obj.trial_full_strain)
	, current_full_strain(old_obj.current_full_strain) {}

void Uniaxial::initialize(const shared_ptr<DomainBase>& D) {
	if(D == nullptr) return;

	if(!D->find_material(base_tag) || D->get_material(base_tag)->get_material_type() != MaterialType::D3) {
		D->disable_material(get_tag());
		return;
	}

	base = get_material(D, base_tag)->get_copy();

	trial_full_strain = current_full_strain.zeros(6);

	current_stiffness = trial_stiffness = initial_stiffness = form_stiffness(base->get_initial_stiffness());
}

double Uniaxial::get_parameter(const ParameterType P) const { return base->get_parameter(P); }

unique_ptr<Material> Uniaxial::get_copy() { return make_unique<Uniaxial>(*this); }

int Uniaxial::update_trial_status(const vec& t_strain) {
	auto& t_stress = base->get_trial_stress();
	auto& t_stiffness = base->get_trial_stiffness();

	if(1 == max_iteration) {
		incre_strain = t_strain - trial_strain;

		trial_full_strain(F1) = trial_strain = t_strain;
		trial_full_strain(F2) -= solve(t_stiffness(F2, F2), t_stress(F2) + t_stiffness(F2, F1) * incre_strain);

		if(base->update_trial_status(trial_full_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		trial_stress = t_stress(F1) - t_stiffness(F1, F2) * solve(t_stiffness(F2, F2), t_stress(F2));
	} else {
		trial_full_strain(F1) = trial_strain = t_strain;

		unsigned counter = 0;

		while(++counter < max_iteration) {
			if(base->update_trial_status(trial_full_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
			trial_full_strain(F2) -= solve(t_stiffness(F2, F2), t_stress(F2));
			const auto error = norm(t_stress(F2));
			suanpan_extra_debug("Uniaxial state determination error: %.4E.\n", error);
			if(error <= tolerance) break;
		}

		if(max_iteration == counter) {
			suanpan_error("Uniaxial cannot converge with in %u iterations.\n", counter);
			return SUANPAN_FAIL;
		}

		trial_stress = t_stress(F1);
	}

	trial_stiffness = form_stiffness(t_stiffness);

	return SUANPAN_SUCCESS;
}

int Uniaxial::clear_status() {
	current_full_strain.zeros(6);
	trial_full_strain.zeros(6);
	Material::clear_status();
	return base->clear_status();
}

int Uniaxial::commit_status() {
	current_full_strain = trial_full_strain;
	Material::commit_status();
	return base->commit_status();
}

int Uniaxial::reset_status() {
	trial_full_strain = current_full_strain;
	Material::reset_status();
	return base->reset_status();
}

vector<vec> Uniaxial::record(const OutputType P) { return base->record(P); }

void Uniaxial::print() {
	suanpan_info("A uniaxial wrapper.\n");
	suanpan_info("current strain: ");
	current_strain.t().print();
	suanpan_info("current stress: ");
	current_stress.t().print();
	base->print();
}
