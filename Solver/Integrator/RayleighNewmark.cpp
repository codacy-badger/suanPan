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

#include "RayleighNewmark.h"
#include <Domain/DomainBase.h>
#include <Element/Utility/MatrixModifier.h>

RayleighNewmark::RayleighNewmark(const unsigned T, const double DA, const double DB, const double A, const double B)
	: Newmark(T, A, B)
	, damping_alpha(DA)
	, damping_beta(DB) {}

void RayleighNewmark::assemble_resistance() {
	const auto& element_pool = get_domain().lock()->get_element_pool();

	suanpan_for_each(element_pool.cbegin(), element_pool.cend(), [&](const shared_ptr<Element>& t_element) { suanpan::damping::rayleigh::apply(t_element, damping_alpha, damping_beta); });

	Newmark::assemble_resistance();
}
