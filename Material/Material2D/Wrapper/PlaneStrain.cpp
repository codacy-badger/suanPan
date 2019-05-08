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

#include "PlaneStrain.h"
#include <Domain/DomainBase.h>

const uvec PlaneStrain::F{0, 1, 3};

PlaneStrain::PlaneStrain(const unsigned T, const unsigned BT)
	: Material2D(T, PlaneType::E, 0.)
	, base_tag(BT)
	, trial_full_strain(6, fill::zeros) {}

PlaneStrain::PlaneStrain(const PlaneStrain& old_obj)
	: Material2D(old_obj)
	, base_tag(old_obj.base_tag)
	, base(old_obj.base == nullptr ? nullptr : old_obj.base->get_copy())
	, trial_full_strain(old_obj.trial_full_strain) {}

void PlaneStrain::initialize(const shared_ptr<DomainBase>& D) {
	if(D == nullptr) return;

	if(!D->find_material(base_tag)) {
		D->disable_material(get_tag());
		return;
	}

	base = get_material(D, base_tag)->get_copy();

	access::rw(density) = base->get_parameter();

	current_stiffness = trial_stiffness = initial_stiffness = base->get_initial_stiffness()(F, F);
}

double PlaneStrain::get_parameter(const ParameterType P) const { return base->get_parameter(P); }

unique_ptr<Material> PlaneStrain::get_copy() { return make_unique<PlaneStrain>(*this); }

int PlaneStrain::update_trial_status(const vec& t_strain) {
	trial_full_strain(F) = trial_strain = t_strain;

	if(base->update_trial_status(trial_full_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

	trial_stress = base->get_trial_stress()(F);

	trial_stiffness = base->get_trial_stiffness()(F, F);

	return SUANPAN_SUCCESS;
}

int PlaneStrain::clear_status() {
	Material::clear_status();
	return base->clear_status();
}

int PlaneStrain::commit_status() {
	Material::commit_status();
	return base->commit_status();
}

int PlaneStrain::reset_status() {
	Material::reset_status();
	return base->reset_status();
}

vector<vec> PlaneStrain::record(const OutputType P) { return base->record(P); }

void PlaneStrain::print() {
	suanpan_info("A plane strain wrapper.\n");
	base->print();
}
