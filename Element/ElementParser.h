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

#ifndef ELEMENTPARSER_H
#define ELEMENTPARSER_H

#include <suanPan.h>

int create_new_element(const shared_ptr<DomainBase>&, istringstream&);

void new_allman(unique_ptr<Element>&, istringstream&);
void new_b21(unique_ptr<Element>&, istringstream&);
void new_b21h(unique_ptr<Element>&, istringstream&);
void new_b31(unique_ptr<Element>&, istringstream&);
void new_c3d20(unique_ptr<Element>&, istringstream&);
void new_c3d4(unique_ptr<Element>&, istringstream&);
void new_c3d8(unique_ptr<Element>&, istringstream&);
void new_c3d8i(unique_ptr<Element>&, istringstream&);
void new_cax3(unique_ptr<Element>&, istringstream&);
void new_cax4(unique_ptr<Element>&, istringstream&);
void new_cax8(unique_ptr<Element>&, istringstream&);
void new_cin3d8(unique_ptr<Element>&, istringstream&);
void new_cinp4(unique_ptr<Element>&, istringstream&);
void new_cp3(unique_ptr<Element>&, istringstream&);
void new_cp4(unique_ptr<Element>&, istringstream&);
void new_cp4i(unique_ptr<Element>&, istringstream&);
void new_cp4r(unique_ptr<Element>&, istringstream&);
void new_cp6(unique_ptr<Element>&, istringstream&);
void new_cp8(unique_ptr<Element>&, istringstream&);
void new_damper01(unique_ptr<Element>&, istringstream&);
void new_damper02(unique_ptr<Element>&, istringstream&);
void new_dkt3(unique_ptr<Element>&, istringstream&);
void new_dkt4(unique_ptr<Element>&, istringstream&);
void new_dkts3(unique_ptr<Element>&, istringstream&);
void new_eb21(unique_ptr<Element>&, istringstream&);
void new_f21(unique_ptr<Element>&, istringstream&);
void new_f21h(unique_ptr<Element>&, istringstream&);
void new_f31(unique_ptr<Element>&, istringstream&);
void new_gcmq(unique_ptr<Element>&, istringstream&);
void new_gcmqg(unique_ptr<Element>&, istringstream&);
void new_gcmqi(unique_ptr<Element>&, istringstream&);
void new_gcmql(unique_ptr<Element>&, istringstream&);
void new_gq12(unique_ptr<Element>&, istringstream&);
void new_mass(unique_ptr<Element>&, istringstream&);
void new_mindlin(unique_ptr<Element>&, istringstream&);
void new_mvlem(unique_ptr<Element>&, istringstream&);
void new_ps(unique_ptr<Element>&, istringstream&);
void new_qe2(unique_ptr<Element>&, istringstream&);
void new_rcp4(unique_ptr<Element>&, istringstream&);
void new_rebarlayer(unique_ptr<Element>&, istringstream&);
void new_rgcmq(unique_ptr<Element>&, istringstream&);
void new_rgcmqg(unique_ptr<Element>&, istringstream&);
void new_rgcmqi(unique_ptr<Element>&, istringstream&);
void new_rgcmql(unique_ptr<Element>&, istringstream&);
void new_s4(unique_ptr<Element>&, istringstream&);
void new_sgcmqg(unique_ptr<Element>&, istringstream&);
void new_sgcmqi(unique_ptr<Element>&, istringstream&);
void new_sgcmql(unique_ptr<Element>&, istringstream&);
void new_srgcmqg(unique_ptr<Element>&, istringstream&);
void new_srgcmqi(unique_ptr<Element>&, istringstream&);
void new_srgcmql(unique_ptr<Element>&, istringstream&);
void new_singlesection2d(unique_ptr<Element>&, istringstream&);
void new_singlesection3d(unique_ptr<Element>&, istringstream&);
void new_spring01(unique_ptr<Element>&, istringstream&);
void new_spring02(unique_ptr<Element>&, istringstream&);
void new_t2d2(unique_ptr<Element>&, istringstream&);
void new_t2d2s(unique_ptr<Element>&, istringstream&);
void new_t3d2(unique_ptr<Element>&, istringstream&);
void new_t3d2s(unique_ptr<Element>&, istringstream&);
void new_tie(unique_ptr<Element>&, istringstream&);

#endif
