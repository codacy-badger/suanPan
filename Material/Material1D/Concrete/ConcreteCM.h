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
/**
 * @class ConcreteCM
 * @brief A ConcreteCM material class.
 *
 * The ConcreteCM class represents a concrete material model based on the Chang & Mander concrete model.
 *
 * doi: 10.1061/(ASCE)0733-9445(1988)114:8(1804)
 *
 * @author tlc
 * @date 19/06/2018
 * @version 0.2.0
 * @file ConcreteCM.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CONCRETECM_H
#define CONCRETECM_H

#include <Material/Material1D/Material1D.h>

class ConcreteCM final : public Material1D {
	enum class Status { NONE, CBACKBONE, TBACKBONE, CUNLOAD, TUNLOAD, CSUBUNLOAD, TSUBUNLOAD, CRELOAD, TRELOAD, CTRANS, TTRANS };

	Status trial_load_status = Status::NONE, current_load_status = Status::NONE;

	const double peak_stress, peak_strain, crack_stress, crack_strain;

	const double MC, NC, MT, NT;

	const bool linear_trans;

	podarray<double> compute_compression_backbone(double) const;
	podarray<double> compute_tension_backbone(double) const;
	podarray<double> compute_compression_unload(double);
	podarray<double> compute_tension_unload(double);
	podarray<double> compute_compression_sub_unload(double);
	podarray<double> compute_tension_sub_unload(double);
	podarray<double> compute_compression_reload(double);
	podarray<double> compute_tension_reload(double);
	void update_compression_reverse(double);
	void update_tension_reverse(double);
	void update_compression_new_reverse(double, double);
	void update_tension_new_reverse(double, double);

	static podarray<double> compute_transition(double, double, double, double, double, double, double, bool);
	static podarray<double> compute_linear_transition(double, double, double, double, double);
public:
	ConcreteCM(unsigned,                           // tag
	           double,                             // peak stress in negative
	           double,                             // crack stress in postive
	           double,                             // MC
	           double,                             // NC
	           double,                             // MT
	           double,                             // NT
	           double = -2E-3 - datum::eps * 1E10, // peak strain in negative
	           double = 1E-4 + datum::eps * 1E10,  // crack strain in postive
	           bool = false,                       // if to use linear transition
	           double = 0.                         // density
	);

	void initialize(const shared_ptr<DomainBase>& = nullptr) override;

	unique_ptr<Material> get_copy() override;

	double get_parameter(ParameterType) const override;

	int update_trial_status(const vec&) override;
	int update_trial_status(const vec&, const vec&) override;

	int clear_status() override;
	int commit_status() override;
	int reset_status() override;

	vector<vec> record(OutputType) override;

	void print() override;
};

#endif

//! @}
