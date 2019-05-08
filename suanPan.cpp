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

#include <suanPan>

// ReSharper disable once CppParameterMayBeConst
int main(int argc, char** argv) {
	if(check_debugger()) return 0;

#ifdef SUANPAN_MAGMA
	magma_init();
#endif

	wall_clock T;
	T.tic();

	argument_parser(argc, argv);

	suanpan_info("Finished in %.3F seconds.\n", T.toc());

#ifdef SUANPAN_MAGMA
	magma_finalize();
#endif

	return 0;
}

void test_mode() { suanpan_info("This is the test mode that shall be used to debug internal codes.\n"); }
