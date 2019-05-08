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
 * @class ConcreteSimple
 * @brief A ConcreteSimple material class.
 *
 * The ConcreteSimple class defines a uniaxial concrete model with simple hysteresis behaviour.
 *
 * @author tlc
 * @date 28/01/2019
 * @version 0.1.0
 * @file ConcreteSimple.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CONCRETESIMPLE_H
#define CONCRETESIMPLE_H

#include <Material/Material1D/Material1D.h>

class ConcreteSimple : public Material1D {
	const double middle_point;

	enum class Status { NONE, CBACKBONE, TBACKBONE, CINNER, TINNER };

	Status trial_flag = Status::NONE, current_flag = Status::NONE;

	virtual podarray<double> compute_compression_backbone(double) const = 0;
	virtual podarray<double> compute_tension_backbone(double) const = 0;
	virtual double compute_compression_residual(double, double) const = 0;
	virtual double compute_tension_residual(double, double) const = 0;
	podarray<double> compute_compression_inner(double) const;
	podarray<double> compute_tension_inner(double) const;
protected:
	const double peak_stress, crack_stress;
	const double peak_strain, crack_strain;
public:
	ConcreteSimple(unsigned,   // tag
	               double,     // peak stress in negative
	               double,     // crack stress in positive
	               double,     // peak strain in negative
	               double,     // crack strain in positive
	               double,     // middle point
	               double = 0. // density
	);

	double get_parameter(ParameterType) const override;

	int update_trial_status(const vec&) override;

	int clear_status() override;
	int commit_status() override;
	int reset_status() override;

	void print() override;
};

#endif

//! @}
