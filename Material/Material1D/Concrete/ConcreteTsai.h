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
 * @class ConcreteTsai
 * @brief A ConcreteTsai material class.
 *
 * @author tlc
 * @date 08/07/2018
 * @version 0.3.0
 * @file ConcreteTsai.h
 * @addtogroup Material-1D
 * @{
 */

#ifndef CONCRETETSAI_H
#define CONCRETETSAI_H

#include "ConcreteSimple.h"

class ConcreteTsai final : public ConcreteSimple {
	const double MC, MT, NC, NT;

	podarray<double> compute_compression_backbone(double) const override;
	podarray<double> compute_tension_backbone(double) const override;
	double compute_compression_residual(double, double) const override;
	double compute_tension_residual(double, double) const override;
public:
	ConcreteTsai(unsigned,   // tag
	             double,     // peak stress in negative
	             double,     // crack stress in positive
	             double,     // MC
	             double,     // NC
	             double,     // MT
	             double,     // NT
	             double,     // middle point
	             double,     // peak strain in negative
	             double,     // crack strain in positive
	             double = 0. // density
	);
	ConcreteTsai(unsigned,   // tag
	             double,     // peak stress in negative
	             double,     // crack stress in positive
	             double,     // M
	             double,     // N
	             double,     // middle point
	             double,     // peak strain in negative
	             double,     // crack strain in positive
	             double = 0. // density
	);

	void initialize(const shared_ptr<DomainBase>& = nullptr) override;

	unique_ptr<Material> get_copy() override;
};

#endif

//! @}
