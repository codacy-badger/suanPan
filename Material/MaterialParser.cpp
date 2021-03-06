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

#include <Domain/DomainBase.h>
#include <Domain/ExternalModule.h>
#include <Material/Material>
#include <Toolbox/utility.h>

int create_new_material(const shared_ptr<DomainBase>& domain, istringstream& command) {
	string material_id;
	if(!get_input(command, material_id)) {
		suanpan_info("create_new_material() needs a tag.\n");
		return 0;
	}

	unique_ptr<Material> new_material = nullptr;

	if(is_equal(material_id, "AFC01")) new_afc01(new_material, command);
	else if(is_equal(material_id, "AFC02")) new_afc02(new_material, command);
	else if(is_equal(material_id, "AFC03")) new_afc03(new_material, command);
	else if(is_equal(material_id, "Axisymmetric")) new_axisymmetric(new_material, command);
	else if(is_equal(material_id, "AxisymmetricElastic")) new_axisymmetricelastic(new_material, command);
	else if(is_equal(material_id, "Bilinear2D")) new_bilinear2d(new_material, command);
	else if(is_equal(material_id, "BilinearCC")) new_bilinearcc(new_material, command);
	else if(is_equal(material_id, "BilinearDP")) new_bilineardp(new_material, command);
	else if(is_equal(material_id, "BilinearElastic1D")) new_bilinearelastic1d(new_material, command);
	else if(is_equal(material_id, "BilinearHoffman")) new_bilinearhoffman(new_material, command);
	else if(is_equal(material_id, "BilinearJ2")) new_bilinearj2(new_material, command);
	else if(is_equal(material_id, "BilinearMises1D")) new_bilinearmises1d(new_material, command);
	else if(is_equal(material_id, "BlatzKo")) new_blatzko(new_material, command);
	else if(is_equal(material_id, "CDP")) new_cdp(new_material, command);
	else if(is_equal(material_id, "ConcreteCM")) new_concretecm(new_material, command);
	else if(is_equal(material_id, "Concrete21")) new_concrete21(new_material, command);
	else if(is_equal(material_id, "Concrete22")) new_concrete22(new_material, command);
	else if(is_equal(material_id, "ConcreteTsai")) new_concretetsai(new_material, command);
	else if(is_equal(material_id, "Elastic1D")) new_elastic1d(new_material, command);
	else if(is_equal(material_id, "Elastic2D")) new_elastic2d(new_material, command);
	else if(is_equal(material_id, "Elastic3D")) new_isotropicelastic3d(new_material, command);
	else if(is_equal(material_id, "ExpCC")) new_expcc(new_material, command);
	else if(is_equal(material_id, "ExpDP")) new_expdp(new_material, command);
	else if(is_equal(material_id, "ExpHoffman")) new_exphoffman(new_material, command);
	else if(is_equal(material_id, "ExpJ2")) new_expj2(new_material, command);
	else if(is_equal(material_id, "ExpMises1D")) new_expmises1d(new_material, command);
	else if(is_equal(material_id, "Flag01")) new_flag01(new_material, command);
	else if(is_equal(material_id, "Flag02")) new_flag02(new_material, command);
	else if(is_equal(material_id, "Gap01")) new_gap01(new_material, command);
	else if(is_equal(material_id, "IsotropicElastic3D")) new_isotropicelastic3d(new_material, command);
	else if(is_equal(material_id, "Maxwell")) new_maxwell(new_material, command);
	else if(is_equal(material_id, "MooneyRivlin")) new_mooneyrivlin(new_material, command);
	else if(is_equal(material_id, "MPF")) new_mpf(new_material, command);
	else if(is_equal(material_id, "MultilinearMises1D")) new_multilinearmises1d(new_material, command);
	else if(is_equal(material_id, "MultilinearElastic1D")) new_multilinearelastic1d(new_material, command);
	else if(is_equal(material_id, "MultilinearJ2")) new_multilinearj2(new_material, command);
	else if(is_equal(material_id, "OrthotropicElastic3D")) new_orthotropicelastic3d(new_material, command);
	else if(is_equal(material_id, "Parallel")) new_parallel(new_material, command);
	else if(is_equal(material_id, "ParabolicCC")) new_paraboliccc(new_material, command);
	else if(is_equal(material_id, "PlaneStrain")) new_planestrain(new_material, command);
	else if(is_equal(material_id, "PlaneStress")) new_planestress(new_material, command);
	else if(is_equal(material_id, "PolyJ2")) new_polyj2(new_material, command);
	else if(is_equal(material_id, "RambergOsgood")) new_rambergosgood(new_material, command);
	else if(is_equal(material_id, "Laminated")) new_laminated(new_material, command);
	else if(is_equal(material_id, "Stacked")) new_stacked(new_material, command);
	else if(is_equal(material_id, "Rebar2D")) new_rebar2d(new_material, command);
	else if(is_equal(material_id, "Rebar3D")) new_rebar3d(new_material, command);
	else if(is_equal(material_id, "Sequential")) new_sequential(new_material, command);
	else if(is_equal(material_id, "Uniaxial")) new_uniaxial(new_material, command);
	else if(is_equal(material_id, "BilinearPerzyna")) new_bilinearperzyna(new_material, command);
	else if(is_equal(material_id, "Viscosity01")) new_viscosity01(new_material, command);
	else if(is_equal(material_id, "Viscosity02")) new_viscosity02(new_material, command);
	else if(is_equal(material_id, "Bilinear1D")) new_bilinear1d(new_material, command);
	else load::object(new_material, domain, material_id, command);

	if(new_material == nullptr || !domain->insert(std::move(new_material))) suanpan_debug("create_new_material() fails to insert new material.\n");

	return 0;
}

void new_afc01(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_afc01() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_afc01() requires a valid elastic modulus.\n");
		return;
	}

	double t_yield_stress, t_hardening, t_unloading;
	double c_yield_stress, c_hardening, c_unloading;
	if(!get_input(command, t_yield_stress)) {
		suanpan_error("new_afc01() requires tension yield stress.\n");
		return;
	}
	if(!get_input(command, t_hardening)) {
		suanpan_error("new_afc01() requires tension hardening modulus.\n");
		return;
	}
	if(!get_input(command, t_unloading)) {
		suanpan_error("new_afc01() requires tension unloading modulus.\n");
		return;
	}
	if(!get_input(command, c_yield_stress)) {
		suanpan_error("new_afc01() requires compression yield stress.\n");
		return;
	}
	if(!get_input(command, c_hardening)) {
		suanpan_error("new_afc01() requires compression hardening modulus.\n");
		return;
	}
	if(!get_input(command, c_unloading)) {
		suanpan_error("new_afc01() requires compression unloading modulus.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_afc01() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_afc01() assumes zero density.\n");

	return_obj = make_unique<AFC>(tag, elastic_modulus, t_yield_stress, t_hardening, t_unloading, c_yield_stress, c_hardening, c_unloading, density);
}

void new_afc02(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_afc02() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_afc02() requires a valid elastic modulus.\n");
		return;
	}

	double t_yield_stress, t_hardening, t_unloading;
	if(!get_input(command, t_yield_stress)) {
		suanpan_error("new_afc02() requires yield stress.\n");
		return;
	}
	if(!get_input(command, t_hardening)) {
		suanpan_error("new_afc02() requires hardening modulus.\n");
		return;
	}
	if(!get_input(command, t_unloading)) {
		suanpan_error("new_afc02() requires unloading modulus.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_afc02() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_afc02() assumes zero density.\n");

	return_obj = make_unique<AFC>(tag, elastic_modulus, t_yield_stress, t_hardening, t_unloading, t_yield_stress, t_hardening, t_unloading, density);
}

void new_afc03(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_afc03() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_afc03() requires a valid elastic modulus.\n");
		return;
	}

	double t_yield_stress, t_hardening, t_unloading;
	double c_yield_stress, c_hardening, c_unloading;
	if(!get_input(command, t_yield_stress)) {
		suanpan_error("new_afc03() requires tension yield stress.\n");
		return;
	}
	if(!get_input(command, t_hardening)) {
		suanpan_error("new_afc03() requires tension hardening modulus.\n");
		return;
	}
	if(!get_input(command, t_unloading)) {
		suanpan_error("new_afc03() requires tension unloading modulus.\n");
		return;
	}
	if(!get_input(command, c_yield_stress)) {
		suanpan_error("new_afc03() requires compression yield stress.\n");
		return;
	}
	if(!get_input(command, c_hardening)) {
		suanpan_error("new_afc03() requires compression hardening modulus.\n");
		return;
	}
	if(!get_input(command, c_unloading)) {
		suanpan_error("new_afc03() requires compression unloading modulus.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_afc03() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_afc03() assumes zero density.\n");

	return_obj = make_unique<AFCN>(tag, elastic_modulus, t_yield_stress, t_hardening, t_unloading, c_yield_stress, c_hardening, c_unloading, density);
}

void new_axisymmetricelastic(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_axisymmetricelastic() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_axisymmetricelastic() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_axisymmetricelastic() requires a valid poissons ratio.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_axisymmetricelastic() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_axistmmetricelastic() assumes zero density.\n");

	return_obj = make_unique<AxisymmetricElastic>(tag, elastic_modulus, poissons_ratio, density);
}

void new_bilinear1d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_bilinear1d() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_bilinear1d() requires a valid elastic modulus.\n");
		return;
	}

	double yield_stress;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_bilinear1d() requires a valid yield stress.\n");
		return;
	}

	auto hardening_ratio = 0.;
	if(!command.eof()) {
		if(!get_input(command, hardening_ratio)) {
			suanpan_error("new_bilinear1d() requires a valid hardening ratio.\n");
			return;
		}
	} else suanpan_debug("new_bilinear1d() assumes zero hardening ratio.\n");

	auto beta = 1.;
	if(!command.eof()) {
		if(!get_input(command, beta)) {
			suanpan_error("new_bilinear1d() requires a valid beta.\n");
			return;
		}
	} else suanpan_debug("new_bilinear1d() assumes isotropic hardening.\n");
	if(beta > 1.) beta = 1.;
	else if(beta < 0.) beta = 0.;

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_bilinear1d() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_bilinear1d() assumes zero density.\n");

	return_obj = make_unique<Bilinear1D>(tag, elastic_modulus, yield_stress, hardening_ratio, beta, density);
}

void new_bilinear2d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_bilinear2d() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_bilinear2d() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_bilinear2d() requires a valid poissons ratio.\n");
		return;
	}

	double yield_stress;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_bilinear2d() requires a valid yield stress.\n");
		return;
	}

	auto hardening_ratio = 0.;
	if(!command.eof()) {
		if(!get_input(command, hardening_ratio)) {
			suanpan_error("new_bilinear2d() requires a valid hardening ratio.\n");
			return;
		}
	} else suanpan_debug("new_bilinear2d() assumes zero hardening ratio.\n");

	auto beta = 1.;
	if(!command.eof()) {
		if(!get_input(command, beta)) {
			suanpan_error("new_bilinear2d() requires a valid beta.\n");
			return;
		}
	} else suanpan_debug("new_bilinear2d() assumes isotropic hardening.\n");

	unsigned material_type = 0;
	if(!command.eof()) {
		if(!get_input(command, material_type)) {
			suanpan_error("new_bilinear2d() requires a valid material type.\n");
			return;
		}
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_bilinear2d() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_bilinear2d() assumes zero density.\n");

	return_obj = make_unique<Bilinear2D>(tag, elastic_modulus, poissons_ratio, yield_stress, hardening_ratio, beta, material_type == 0 ? PlaneType::S : PlaneType::E, density);
}

void new_bilinearj2(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_bilinearj2() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_bilinearj2() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_bilinearj2() requires a valid poissons ratio.\n");
		return;
	}

	double yield_stress;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_bilinearj2() requires a valid yield stress.\n");
		return;
	}

	auto hardening_ratio = 0.;
	if(!command.eof()) {
		if(!get_input(command, hardening_ratio)) {
			suanpan_error("new_bilinearj2() requires a valid hardening ratio.\n");
			return;
		}
	} else suanpan_debug("new_bilinearj2() assumes zero hardening ratio.\n");

	auto beta = 1.;
	if(!command.eof()) {
		if(!get_input(command, beta)) {
			suanpan_error("new_bilinearj2() requires a valid beta.\n");
			return;
		}
	} else suanpan_debug("new_bilinearj2() assumes isotropic hardening.\n");

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_bilinearj2() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_bilinearj2() assumes zero density.\n");

	return_obj = make_unique<BilinearJ2>(tag, elastic_modulus, poissons_ratio, yield_stress, hardening_ratio, beta, density);
}

void new_bilinearelastic1d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_bilinearelastic1d() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_bilinearelastic1d() requires a valid elastic modulus.\n");
		return;
	}

	double yield_stress;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_bilinearelastic1d() requires a valid yield stress.\n");
		return;
	}

	auto hardening_ratio = 0.;
	if(!command.eof()) {
		if(!get_input(command, hardening_ratio)) {
			suanpan_error("new_bilinearelastic1d() requires a valid hardening ratio.\n");
			return;
		}
	} else suanpan_debug("new_bilinearelastic1d() assumes zero hardening ratio.\n");

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_bilinearelastic1d() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_bilinear1d() assumes zero density.\n");

	return_obj = make_unique<BilinearElastic1D>(tag, elastic_modulus, yield_stress, hardening_ratio, density);
}

void new_bilinearmises1d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_bilinearmises1d() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_bilinearmises1d() requires a valid elastic modulus.\n");
		return;
	}

	double yield_stress;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_bilinearmises1d() requires a valid yield stress.\n");
		return;
	}

	auto hardening_ratio = 0.;
	if(!command.eof()) {
		if(!get_input(command, hardening_ratio)) {
			suanpan_error("new_bilinearmises1d() requires a valid hardening ratio.\n");
			return;
		}
	} else suanpan_debug("new_bilinearmises1d() assumes zero hardening ratio.\n");

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_bilinearmises1d() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_bilinearmises1d() assumes zero density.\n");

	return_obj = make_unique<BilinearMises1D>(tag, elastic_modulus, yield_stress, hardening_ratio, density);
}

void new_blatzko(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_blatzko() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_blatzko() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_blatzko() requires a valid poissons ratio.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_blatzko() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_blatzko() assumes zero density.\n");

	return_obj = make_unique<BlatzKo>(tag, elastic_modulus, poissons_ratio, density);
}

void new_cdp(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_cdp() requires a valid tag.\n");
		return;
	}
	vec para_pool(std::initializer_list<double>{3E4, .2, 3., 30., 5E-4, 5E-2, .2, 2., .5, .65, .2, 1.16, .5, 2400E-12});

	auto idx = 0;
	double para;
	while(!command.eof() && idx < 14) if(get_input(command, para)) para_pool(idx++) = para;

	return_obj = make_unique<CDP>(tag, para_pool(0), para_pool(1), para_pool(2), para_pool(3), para_pool(4), para_pool(5), para_pool(6), para_pool(7), para_pool(8), para_pool(9), para_pool(10), para_pool(11), para_pool(12), para_pool(13));
}

void new_concretecm(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_concretecm() requires a valid tag.\n");
		return;
	}

	double peak_stress;
	if(!get_input(command, peak_stress)) {
		suanpan_error("new_concretecm() requires a valid compression stress.\n");
		return;
	}

	double crack_stress;
	if(!get_input(command, crack_stress)) {
		suanpan_error("new_concretecm() requires a valid tension stress.\n");
		return;
	}

	double MC, NC, MT, NT;
	if(!get_input(command, MC)) {
		suanpan_error("new_concretecm() requires a valid parameter.\n");
		return;
	}
	if(!get_input(command, NC)) {
		suanpan_error("new_concretecm() requires a valid parameter.\n");
		return;
	}
	if(!get_input(command, MT)) {
		suanpan_error("new_concretecm() requires a valid parameter.\n");
		return;
	}
	if(!get_input(command, NT)) {
		suanpan_error("new_concretecm() requires a valid parameter.\n");
		return;
	}

	auto peak_strain = 2E-3;
	if(!command.eof() && !get_input(command, peak_strain)) {
		suanpan_error("new_concretecm() requires a valid tension stress.\n");
		return;
	}

	auto crack_strain = 1E-4;
	if(!command.eof() && !get_input(command, crack_strain)) {
		suanpan_error("new_concretecm() requires a valid tension stress.\n");
		return;
	}

	string linear_trans = "false";
	if(!command.eof() && !get_input(command, linear_trans)) {
		suanpan_error("new_concretecm() requires a valid transition switch.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_concretecm() requires a valid density.\n");
			return;
		}
	} else suanpan_extra_debug("new_concretecm() assumes zero density.\n");

	return_obj = make_unique<ConcreteCM>(tag, peak_stress, crack_stress, MC, NC, MT, NT, peak_strain + datum::eps * 1E10, crack_strain + datum::eps * 1E10, is_true(linear_trans), density);
}

void new_concrete21(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_concrete21() requires a valid tag.\n");
		return;
	}

	vector<double> para;
	double input;
	while(!command.eof() && get_input(command, input)) para.emplace_back(input);

	if(para.size() == 9) return_obj = make_unique<Concrete21>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8], 0.);
	else if(para.size() == 10) return_obj = make_unique<Concrete21>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8], para[9]);
	else suanpan_error("new_concrete21() requires 9 or 10 double inputs.\n");
}

void new_concrete22(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_concrete22() requires a valid tag.\n");
		return;
	}

	vector<double> para;
	double input;
	while(!command.eof() && get_input(command, input)) para.emplace_back(input);

	if(para.size() == 11) return_obj = make_unique<Concrete22>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8], para[9], para[10], 0.);
	else if(para.size() == 12) return_obj = make_unique<Concrete22>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8], para[9], para[10], para[11]);
	else suanpan_error("new_concrete22() requires 11 or 12 double inputs.\n");
}

void new_concretetsai(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_concretetsai() requires a valid tag.\n");
		return;
	}

	vector<double> para;
	double input;
	while(!command.eof() && get_input(command, input)) para.emplace_back(input);

	if(para.size() == 7) return_obj = make_unique<ConcreteTsai>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], 0.);
	else if(para.size() == 8) return_obj = make_unique<ConcreteTsai>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7]);
	else if(para.size() == 9) return_obj = make_unique<ConcreteTsai>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8], 0.);
	else if(para.size() == 10) return_obj = make_unique<ConcreteTsai>(tag, para[0], para[1], para[2], para[3], para[4], para[5], para[6], para[7], para[8], para[9]);
	else suanpan_error("new_concretetsai() requires 7, 8, 9 or 10 double inputs.\n");
}

void new_elastic1d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_elastic1d() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_elastic1d() requires a valid elastic modulus.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_elastic1d() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_elastic1d() assumes zero density.\n");

	return_obj = make_unique<Elastic1D>(tag, elastic_modulus, density);
}

void new_elastic2d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_elastic2d() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_elastic2d() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_elastic2d() requires a valid poissons ratio.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_elastic2d() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_elastic2d() assumes zero density.\n");

	auto material_type = 0;
	if(!command.eof()) {
		if(!get_input(command, material_type)) {
			suanpan_error("new_elastic2d() requires a valid material type.\n");
			return;
		}
	} else suanpan_debug("new_elastic2d() assumes plane stress.\n");

	return_obj = make_unique<Elastic2D>(tag, elastic_modulus, poissons_ratio, density, material_type == 0 ? PlaneType::S : PlaneType::E);
}

void new_isotropicelastic3d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_elastic3d() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_elastic3d() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_elastic3d() requires a valid poissons ratio.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_elastic3d() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_elastic3d() assumes zero density.\n");

	return_obj = make_unique<IsotropicElastic3D>(tag, elastic_modulus, poissons_ratio, density);
}

void new_expcc(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_expcc() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_expcc() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_expcc() requires a valid poisson's ratio.\n");
		return;
	}

	double beta, m, pt, a0, e0, lambda, kappa;
	if(!get_input(command, beta)) {
		suanpan_error("new_expcc() requires a valid beta.\n");
		return;
	}
	if(!get_input(command, m)) {
		suanpan_error("new_expcc() requires a valid radius ratio.\n");
		return;
	}
	if(!get_input(command, pt)) {
		suanpan_error("new_expcc() requires a valid tensile yield strength.\n");
		return;
	}
	if(!get_input(command, a0)) {
		suanpan_error("new_expcc() requires a valid initial a_0.\n");
		return;
	}
	if(!get_input(command, e0)) {
		suanpan_error("new_expcc() requires a valid initial void ratio.\n");
		return;
	}
	if(!get_input(command, lambda)) {
		suanpan_error("new_expcc() requires a valid lambda.\n");
		return;
	}
	if(!get_input(command, kappa)) {
		suanpan_error("new_expcc() requires a valid kappa.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof() && !get_input(command, density)) {
		suanpan_error("new_expcc() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<ExpCC>(tag, elastic_modulus, poissons_ratio, beta, m, pt, a0, e0, lambda, kappa, density);
}

void new_expdp(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_expdp() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_expdp() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_expdp() requires a valid poisson's ratio.\n");
		return;
	}

	double eta_yield, eta_flow, xi, cohesion, cohesion_a, cohesion_b;
	if(!get_input(command, eta_yield)) {
		suanpan_error("new_expdp() requires a valid eta for yielding criterion.\n");
		return;
	}
	if(!get_input(command, eta_flow)) {
		suanpan_error("new_expdp() requires a valid eta for plasticity flow rule.\n");
		return;
	}
	if(!get_input(command, xi)) {
		suanpan_error("new_expdp() requires a valid xi.\n");
		return;
	}
	if(!get_input(command, cohesion)) {
		suanpan_error("new_expdp() requires a valid cohesion.\n");
		return;
	}
	if(!get_input(command, cohesion_a)) {
		suanpan_error("new_expdp() requires a valid cohesion.\n");
		return;
	}
	if(!get_input(command, cohesion_b)) {
		suanpan_error("new_expdp() requires a valid cohesion.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof() && !get_input(command, density)) {
		suanpan_error("new_expdp() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<ExpDP>(tag, elastic_modulus, poissons_ratio, eta_yield, eta_flow, xi, cohesion, cohesion_a, cohesion_b, density);
}

void new_expj2(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_expj2() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_expj2() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_expj2() requires a valid poissons ratio.\n");
		return;
	}

	double yield_stress, a, b;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_expj2() requires a valid yield stress.\n");
		return;
	}
	if(!get_input(command, a)) {
		suanpan_error("new_expj2() requires a valid a.\n");
		return;
	}
	if(!get_input(command, b)) {
		suanpan_error("new_expj2() requires a valid b.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof() && !get_input(command, density)) {
		suanpan_error("new_expj2() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<ExpJ2>(tag, elastic_modulus, poissons_ratio, yield_stress, a, b, density);
}

void new_exphoffman(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_linearhoffman() requires a valid tag.\n");
		return;
	}

	vec modulus(6);
	for(unsigned I = 0; I < modulus.n_elem; ++I)
		if(!get_input(command, modulus(I))) {
			suanpan_error("new_linearhoffman() requires a valid modulus.\n");
			return;
		}

	vec poissons_ratio(3);
	for(unsigned I = 0; I < poissons_ratio.n_elem; ++I)
		if(!get_input(command, poissons_ratio(I))) {
			suanpan_error("new_linearhoffman() requires a valid poisson's ratio.\n");
			return;
		}

	vec stress(9);
	for(unsigned I = 0; I < stress.n_elem; ++I)
		if(!get_input(command, stress(I))) {
			suanpan_error("new_linearhoffman() requires a valid yield stress.\n");
			return;
		}

	double a, b;
	if(!get_input(command, a)) {
		suanpan_error("new_linearhoffman() requires a valid a.\n");
		return;
	}
	if(!get_input(command, b)) {
		suanpan_error("new_linearhoffman() requires a valid b.\n");
		return;
	}

	auto density = 0.;
	if(command.eof()) suanpan_debug("new_linearhoffman() assumes zero density.\n");
	else if(!get_input(command, density)) {
		suanpan_error("new_linearhoffman() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<ExpHoffman>(tag, std::move(modulus), std::move(poissons_ratio), std::move(stress), a, b, density);
}

void new_expmises1d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_expmises1d() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_expmises1d() requires a valid elastic modulus.\n");
		return;
	}

	double yield_stress, a, b;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_expmises1d() requires a valid yield stress.\n");
		return;
	}
	if(!get_input(command, a)) {
		suanpan_error("new_expmises1d() requires a valid a.\n");
		return;
	}
	if(!get_input(command, b)) {
		suanpan_error("new_expmises1d() requires a valid b.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof() && !get_input(command, density)) {
		suanpan_error("new_expmises1d() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<ExpMises1D>(tag, elastic_modulus, yield_stress, a, b, density);
}

void new_flag01(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_flag() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_flag() requires a valid elastic modulus.\n");
		return;
	}

	double yield_stress;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_flag() requires a valid yield stress.\n");
		return;
	}

	double residual;
	if(!get_input(command, residual)) {
		suanpan_error("new_flag() requires a valid residual stress.\n");
		return;
	}

	auto hardening_ratio = 0.;
	if(!command.eof()) {
		if(!get_input(command, hardening_ratio)) {
			suanpan_error("new_flag() requires a valid hardening ratio.\n");
			return;
		}
	} else suanpan_debug("new_flag() assumes zero hardening ratio.\n");

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_bilinear1d() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_bilinear1d() assumes zero density.\n");

	return_obj = make_unique<Flag>(tag, elastic_modulus, yield_stress, residual, hardening_ratio, density);
}

void new_flag02(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_flag() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_flag() requires a valid elastic modulus.\n");
		return;
	}

	double t_yield_stress;
	if(!get_input(command, t_yield_stress)) {
		suanpan_error("new_flag() requires a valid yield stress.\n");
		return;
	}

	double t_residual;
	if(!get_input(command, t_residual)) {
		suanpan_error("new_flag() requires a valid residual stress.\n");
		return;
	}

	double t_hardening_ratio;
	if(!get_input(command, t_hardening_ratio)) {
		suanpan_error("new_flag() requires a valid hardening ratio.\n");
		return;
	}

	double c_yield_stress;
	if(!get_input(command, c_yield_stress)) {
		suanpan_error("new_flag() requires a valid yield stress.\n");
		return;
	}

	double c_residual;
	if(!get_input(command, c_residual)) {
		suanpan_error("new_flag() requires a valid residual stress.\n");
		return;
	}

	double c_hardening_ratio;
	if(!get_input(command, c_hardening_ratio)) {
		suanpan_error("new_flag() requires a valid hardening ratio.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_bilinear1d() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_bilinear1d() assumes zero density.\n");

	return_obj = make_unique<Flag>(tag, elastic_modulus, t_yield_stress, t_residual, t_hardening_ratio, c_yield_stress, c_residual, c_hardening_ratio, density);
}

void new_gap01(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_gap01() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_gap01() requires a valid elastic modulus.\n");
		return;
	}

	double yield_stress;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_gap01() requires a valid yield stress.\n");
		return;
	}

	auto gap_strain = 0.;
	if(!command.eof() && !get_input(command, gap_strain)) {
		suanpan_error("new_gap01() requires a valid hardening ratio.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof() && !get_input(command, density)) {
		suanpan_error("new_gap01() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<Gap01>(tag, elastic_modulus, yield_stress, gap_strain, density);
}

void new_bilinearcc(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_bilinearcc() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_bilinearcc() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_bilinearcc() requires a valid poisson's ratio.\n");
		return;
	}

	double beta, m, pt, a, a_slope;
	if(!get_input(command, beta)) {
		suanpan_error("new_bilinearcc() requires a valid beta.\n");
		return;
	}
	if(!get_input(command, m)) {
		suanpan_error("new_bilinearcc() requires a valid radius ratio.\n");
		return;
	}
	if(!get_input(command, pt)) {
		suanpan_error("new_bilinearcc() requires a valid tensile yield strength.\n");
		return;
	}
	if(!get_input(command, a)) {
		suanpan_error("new_bilinearcc() requires a valid initial size.\n");
		return;
	}
	if(!get_input(command, a_slope)) {
		suanpan_error("new_bilinearcc() requires a valid hardening slope.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof() && !get_input(command, density)) {
		suanpan_error("new_bilinearcc() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<BilinearCC>(tag, elastic_modulus, poissons_ratio, beta, m, pt, a, a_slope, density);
}

void new_bilineardp(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_lineardp() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_lineardp() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_lineardp() requires a valid poisson's ratio.\n");
		return;
	}

	double eta_yield, eta_flow, xi, cohesion, cohesion_slope;
	if(!get_input(command, eta_yield)) {
		suanpan_error("new_lineardp() requires a valid eta for yielding criterion.\n");
		return;
	}
	if(!get_input(command, eta_flow)) {
		suanpan_error("new_lineardp() requires a valid eta for plasticity flow rule.\n");
		return;
	}
	if(!get_input(command, xi)) {
		suanpan_error("new_lineardp() requires a valid xi.\n");
		return;
	}
	if(!get_input(command, cohesion)) {
		suanpan_error("new_lineardp() requires a valid cohesion.\n");
		return;
	}
	if(!get_input(command, cohesion_slope)) {
		suanpan_error("new_lineardp() requires a valid cohesion.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof() && !get_input(command, density)) {
		suanpan_error("new_lineardp() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<BilinearDP>(tag, elastic_modulus, poissons_ratio, eta_yield, eta_flow, xi, cohesion, cohesion_slope, density);
}

void new_bilinearhoffman(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_linearhoffman() requires a valid tag.\n");
		return;
	}

	vec modulus(6);
	for(unsigned I = 0; I < modulus.n_elem; ++I)
		if(!get_input(command, modulus(I))) {
			suanpan_error("new_linearhoffman() requires a valid modulus.\n");
			return;
		}

	vec poissons_ratio(3);
	for(unsigned I = 0; I < poissons_ratio.n_elem; ++I)
		if(!get_input(command, poissons_ratio(I))) {
			suanpan_error("new_linearhoffman() requires a valid poisson's ratio.\n");
			return;
		}

	vec stress(9);
	for(unsigned I = 0; I < stress.n_elem; ++I)
		if(!get_input(command, stress(I))) {
			suanpan_error("new_linearhoffman() requires a valid yield stress.\n");
			return;
		}

	auto hardening = 0.;
	if(command.eof()) suanpan_debug("new_linearhoffman() assumes zero hardening.\n");
	else if(!get_input(command, hardening)) {
		suanpan_error("new_linearhoffman() requires a valid hardening ratio.\n");
		return;
	}

	auto density = 0.;
	if(command.eof()) suanpan_debug("new_linearhoffman() assumes zero density.\n");
	else if(!get_input(command, density)) {
		suanpan_error("new_linearhoffman() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<BilinearHoffman>(tag, std::move(modulus), std::move(poissons_ratio), std::move(stress), hardening, density);
}

void new_maxwell(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_maxwell() requires a valid tag.\n");
		return;
	}

	unsigned damper_tag;
	if(!get_input(command, damper_tag)) {
		suanpan_error("new_maxwell() requires a valid tag.\n");
		return;
	}

	unsigned spring_tag;
	if(!get_input(command, spring_tag)) {
		suanpan_error("new_maxwell() requires a valid tag.\n");
		return;
	}

	string matrix = "false";
	if(!command.eof() && !get_input(command, matrix)) {
		suanpan_error("new_maxwell() requires a valid algorithm switch.\n");
		return;
	}

	unsigned proceed = 0;
	if(!command.eof() && !get_input(command, proceed)) {
		suanpan_error("new_maxwell() requires a valid algorithm switch.\n");
		return;
	}

	return_obj = make_unique<Maxwell>(tag, damper_tag, spring_tag, is_true(matrix), proceed);
}

void new_mooneyrivlin(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_mooneyrivlin() requires a valid tag.\n");
		return;
	}

	double bulk_modulus;
	if(!get_input(command, bulk_modulus)) {
		suanpan_error("new_mooneyrivlin() requires a valid bulk modulus.\n");
		return;
	}

	double a10;
	if(!get_input(command, a10)) {
		suanpan_error("new_mooneyrivlin() requires a valid a10.\n");
		return;
	}

	double a01;
	if(!get_input(command, a01)) {
		suanpan_error("new_mooneyrivlin() requires a valid a01.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_mooneyrivlin() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_mooneyrivlin() assumes zero density.\n");

	return_obj = make_unique<MooneyRivlin>(tag, bulk_modulus, a10, a01, density);
}

void new_multilinearmises1d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_multilinearmises1d() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_multilinearmises1d() requires a valid elastic modulus.\n");
		return;
	}

	auto density = 0.;
	if(!get_input(command, density)) {
		suanpan_error("new_multilinearmises1d() requires a valid density.\n");
		return;
	}

	vector<double> p_strain, p_stress;
	double c_value;
	while(!command.eof()) {
		if(!get_input(command, c_value)) {
			suanpan_error("new_multilinearmises1d() requires a valid plastic strain.\n");
			return;
		}
		p_strain.emplace_back(c_value);
		if(!get_input(command, c_value)) {
			suanpan_error("new_multilinearmises1d() requires a valid plastic stress.\n");
			return;
		}
		p_stress.emplace_back(c_value);
	}

	return_obj = make_unique<MultilinearMises1D>(tag, elastic_modulus, join_rows(vec{p_strain}, vec{p_stress}), density);
}

void new_multilinearj2(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_multilinearj2() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_multilinearj2() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_multilinearj2() requires a valid poissons ratio.\n");
		return;
	}

	auto density = 0.;
	if(!get_input(command, density)) {
		suanpan_error("new_multilinearj2() requires a valid density.\n");
		return;
	}

	vector<double> p_strain, p_stress;
	double c_value;
	while(!command.eof()) {
		if(!get_input(command, c_value)) {
			suanpan_error("new_multilinearj2() requires a valid plastic strain.\n");
			return;
		}
		p_strain.emplace_back(c_value);
		if(!get_input(command, c_value)) {
			suanpan_error("new_multilinearj2() requires a valid plastic stress.\n");
			return;
		}
		p_stress.emplace_back(c_value);
	}

	return_obj = make_unique<MultilinearJ2>(tag, elastic_modulus, poissons_ratio, join_rows(vec{p_strain}, vec{p_stress}), density);
}

void new_multilinearelastic1d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_multilinearelastic1d() requires a valid tag.\n");
		return;
	}

	auto density = 0.;
	if(!get_input(command, density)) {
		suanpan_error("new_multilinearelastic1d() requires a valid density.\n");
		return;
	}

	vector<double> p_strain, p_stress;
	double c_value;
	while(!command.eof()) {
		if(!get_input(command, c_value)) {
			suanpan_error("new_multilinearelastic1d() requires a valid plastic strain.\n");
			return;
		}
		p_strain.emplace_back(c_value);
		if(!get_input(command, c_value)) {
			suanpan_error("new_multilinearelastic1d() requires a valid plastic stress.\n");
			return;
		}
		p_stress.emplace_back(c_value);
	}

	return_obj = make_unique<MultilinearElastic1D>(tag, join_rows(vec{p_strain}, vec{p_stress}), density);
}

void new_mpf(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_mpf() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_mpf() requires a valid elastic modulus.\n");
		return;
	}

	double yield_stress;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_mpf() requires a valid yield stress.\n");
		return;
	}

	auto hardening_ratio = .05;
	if(!command.eof() && !get_input(command, hardening_ratio)) {
		suanpan_error("new_mpf() requires a valid hardening ratio.\n");
		return;
	}

	auto R0 = 20.;
	if(!command.eof() && !get_input(command, R0)) {
		suanpan_error("new_mpf() requires a valid R0.\n");
		return;
	}

	auto A1 = 18.5;
	if(!command.eof() && !get_input(command, A1)) {
		suanpan_error("new_mpf() requires a valid A1.\n");
		return;
	}

	auto A2 = .15;
	if(!command.eof() && !get_input(command, A2)) {
		suanpan_error("new_mpf() requires a valid A2.\n");
		return;
	}

	auto A3 = .01;
	if(!command.eof() && !get_input(command, A3)) {
		suanpan_error("new_mpf() requires a valid A3.\n");
		return;
	}

	auto A4 = 7.;
	if(!command.eof() && !get_input(command, A4)) {
		suanpan_error("new_mpf() requires a valid A4.\n");
		return;
	}

	string iso = "false";
	if(!command.eof() && !get_input(command, iso)) {
		suanpan_error("new_mpf() requires a valid isotropic hardening switch.\n");
		return;
	}

	string con = "false";
	if(!command.eof() && !get_input(command, con)) {
		suanpan_error("new_mpf() requires a valid constant radius switch.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof() && !get_input(command, density)) {
		suanpan_error("new_mpf() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<MPF>(tag, elastic_modulus, yield_stress, hardening_ratio, R0, A1, A2, A3, A4, is_true(iso), is_true(con), density);
}

void new_orthotropicelastic3d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_orthotropicelastic3d() requires a valid tag.\n");
		return;
	}

	vec modulus(6);
	for(unsigned I = 0; I < modulus.n_elem; ++I)
		if(!get_input(command, modulus(I))) {
			suanpan_error("new_orthotropicelastic3d() requires a valid modulus.\n");
			return;
		}

	vec poissons_ratio(3);
	for(unsigned I = 0; I < poissons_ratio.n_elem; ++I)
		if(!get_input(command, poissons_ratio(I))) {
			suanpan_error("new_orthotropicelastic3d() requires a valid poisson's ratio.\n");
			return;
		}

	auto density = 0.;
	if(command.eof()) suanpan_debug("new_orthotropicelastic3d() assumes zero density.\n");
	else if(!get_input(command, density)) {
		suanpan_error("new_orthotropicelastic3d() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<OrthotropicElastic3D>(tag, std::move(modulus), std::move(poissons_ratio), density);
}

void new_polyj2(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_polyj2() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_polyj2() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_polyj2() requires a valid poissons ratio.\n");
		return;
	}

	auto density = 0.;
	if(!get_input(command, density)) {
		suanpan_error("new_polyj2() requires a valid density.\n");
		return;
	}

	vector<double> p_para;
	double c_value;
	while(!command.eof()) {
		if(!get_input(command, c_value)) {
			suanpan_error("new_polyj2() requires a valid plastic strain.\n");
			return;
		}
		p_para.emplace_back(c_value);
	}

	return_obj = make_unique<PolyJ2>(tag, elastic_modulus, poissons_ratio, vec{p_para}, density);
}

void new_rambergosgood(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_rambergosgood() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_rambergosgood() requires a valid elastic modulus.\n");
		return;
	}

	double yield_stress;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_rambergosgood() requires a valid yield stress.\n");
		return;
	}

	auto offset = 1.;
	if(!command.eof() && !get_input(command, offset)) {
		suanpan_error("new_rambergosgood() requires a valid offset.\n");
		return;
	}

	auto n = 10.;
	if(!command.eof() && !get_input(command, n)) {
		suanpan_error("new_rambergosgood() requires a valid n.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof() && !get_input(command, density)) {
		suanpan_error("new_rambergosgood() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<RambergOsgood>(tag, elastic_modulus, yield_stress, offset, n, density);
}

void new_laminated(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_laminated() requires a valid tag.\n");
		return;
	}

	uword c_value;
	vector<uword> mat_tag;
	while(!command.eof() && get_input(command, c_value)) mat_tag.emplace_back(c_value);

	return_obj = make_unique<Laminated>(tag, uvec(mat_tag));
}

void new_stacked(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_stacked() requires a valid tag.\n");
		return;
	}

	uword c_value;
	vector<uword> mat_tag;
	while(!command.eof() && get_input(command, c_value)) mat_tag.emplace_back(c_value);

	return_obj = make_unique<Stacked>(tag, uvec(mat_tag));
}

void new_rebar2d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_rebarlayer() requires a valid tag.\n");
		return;
	}

	unsigned major_tag, minor_tag;
	if(!get_input(command, major_tag)) {
		suanpan_error("new_rebarlayer() requires a valid material tag.\n");
		return;
	}
	if(!get_input(command, minor_tag)) {
		suanpan_error("new_rebarlayer() requires a valid material tag.\n");
		return;
	}

	double major_ratio, minor_ratio;
	if(!get_input(command, major_ratio)) {
		suanpan_error("new_rebarlayer() requires a valid reinforcement ratio.\n");
		return;
	}
	if(!get_input(command, minor_ratio)) {
		suanpan_error("new_rebarlayer() requires a valid reinforcement ratio.\n");
		return;
	}

	return_obj = make_unique<Rebar2D>(tag, major_tag, minor_tag, major_ratio, minor_ratio);
}

void new_rebar3d(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_rebar3d() requires a valid tag.\n");
		return;
	}

	unsigned tag_x, tag_y, tag_z;
	if(!get_input(command, tag_x)) {
		suanpan_error("new_rebar3d() requires a valid material tag.\n");
		return;
	}
	if(!get_input(command, tag_y)) {
		suanpan_error("new_rebar3d() requires a valid material tag.\n");
		return;
	}
	if(!get_input(command, tag_z)) {
		suanpan_error("new_rebar3d() requires a valid material tag.\n");
		return;
	}

	double ratio_x, ratio_y, ratio_z;
	if(!get_input(command, ratio_x)) {
		suanpan_error("new_rebar3d() requires a valid reinforcement ratio.\n");
		return;
	}
	if(!get_input(command, ratio_y)) {
		suanpan_error("new_rebar3d() requires a valid reinforcement ratio.\n");
		return;
	}
	if(!get_input(command, ratio_z)) {
		suanpan_error("new_rebar3d() requires a valid reinforcement ratio.\n");
		return;
	}

	return_obj = make_unique<Rebar3D>(tag, tag_x, tag_y, tag_z, ratio_x, ratio_y, ratio_z);
}

void new_viscosity01(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_viscosity01() requires a valid tag.\n");
		return;
	}

	double alpha;
	if(!get_input(command, alpha)) {
		suanpan_error("new_viscosity01() requires a valid alpha.\n");
		return;
	}

	double damping;
	if(!get_input(command, damping)) {
		suanpan_error("new_viscosity01() requires a valid damping coefficient.\n");
		return;
	}

	auto limit = 1.;
	if(!command.eof() && !get_input(command, limit)) {
		suanpan_error("new_viscosity01() requires a valid limit.\n");
		return;
	}

	return_obj = make_unique<Viscosity01>(tag, alpha, damping, limit);
}

void new_viscosity02(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_viscosity02() requires a valid tag.\n");
		return;
	}

	double alpha;
	if(!get_input(command, alpha)) {
		suanpan_error("new_viscosity02() requires a valid alpha.\n");
		return;
	}

	double damping_a;
	if(!get_input(command, damping_a)) {
		suanpan_error("new_viscosity02() requires a valid damping coefficient for the first quadrant.\n");
		return;
	}

	auto damping_b = damping_a;
	if(!command.eof() && !get_input(command, damping_b)) {
		suanpan_error("new_viscosity02() requires a valid damping coefficient for the second quadrant.\n");
		return;
	}

	auto damping_c = damping_a;
	if(!command.eof() && !get_input(command, damping_c)) {
		suanpan_error("new_viscosity02() requires a valid damping coefficient for the third quadrant.\n");
		return;
	}

	auto damping_d = damping_a;
	if(!command.eof() && !get_input(command, damping_d)) {
		suanpan_error("new_viscosity02() requires a valid damping coefficient for the fourth quadrant.\n");
		return;
	}

	auto gap_a = 1E3;
	if(!command.eof() && !get_input(command, gap_a)) {
		suanpan_error("new_viscosity02() requires a valid gap size for strain axis.\n");
		return;
	}

	auto gap_b = 1E3;
	if(!command.eof() && !get_input(command, gap_b)) {
		suanpan_error("new_viscosity02() requires a valid gap size for strain rate axis.\n");
		return;
	}

	auto limit = 1.;
	if(!command.eof() && !get_input(command, limit)) {
		suanpan_error("new_viscosity02() requires a valid limit.\n");
		return;
	}

	return_obj = make_unique<Viscosity02>(tag, alpha, damping_a, damping_b, damping_c, damping_d, gap_a, gap_b, limit);
}

void new_bilinearperzyna(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_bilinearperzyna() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_bilinearperzyna() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_bilinearperzyna() requires a valid poissons ratio.\n");
		return;
	}

	double yield_stress;
	if(!get_input(command, yield_stress)) {
		suanpan_error("new_bilinearperzyna() requires a valid yield stress.\n");
		return;
	}

	double mu, epsilion;
	if(!get_input(command, mu)) {
		suanpan_error("new_bilinearperzyna() requires a valid mu.\n");
		return;
	}
	if(!get_input(command, epsilion)) {
		suanpan_error("new_bilinearperzyna() requires a valid epsilion.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof()) {
		if(!get_input(command, density)) {
			suanpan_error("new_bilinearperzyna() requires a valid density.\n");
			return;
		}
	} else suanpan_debug("new_bilinearperzyna() assumes zero density.\n");

	return_obj = make_unique<BilinearPerzyna>(tag, elastic_modulus, poissons_ratio, yield_stress, mu, epsilion, density);
}

void new_uniaxial(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_uniaxial() requires a valid tag.\n");
		return;
	}

	unsigned full_tag;
	if(!get_input(command, full_tag)) {
		suanpan_error("new_uniaxial() requires a valid reference material tag.\n");
		return;
	}

	auto max_iteration = 1;
	if(!command.eof() && !get_input(command, max_iteration)) {
		suanpan_error("new_uniaxial() requires a number for maximum iteration.\n");
		return;
	}

	return_obj = make_unique<Uniaxial>(tag, full_tag, max_iteration);
}

void new_axisymmetric(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_axisymmetric() requires a valid tag.\n");
		return;
	}

	unsigned full_tag;
	if(!get_input(command, full_tag)) {
		suanpan_error("new_axisymmetric() requires a valid reference material tag.\n");
		return;
	}

	return_obj = make_unique<Axisymmetric>(tag, full_tag);
}

void new_paraboliccc(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_paraboliccc() requires a valid tag.\n");
		return;
	}

	double elastic_modulus;
	if(!get_input(command, elastic_modulus)) {
		suanpan_error("new_paraboliccc() requires a valid elastic modulus.\n");
		return;
	}

	double poissons_ratio;
	if(!get_input(command, poissons_ratio)) {
		suanpan_error("new_paraboliccc() requires a valid poisson's ratio.\n");
		return;
	}

	double beta, m, pt, a, a_slope;
	if(!get_input(command, beta)) {
		suanpan_error("new_paraboliccc() requires a valid beta.\n");
		return;
	}
	if(!get_input(command, m)) {
		suanpan_error("new_paraboliccc() requires a valid radius ratio.\n");
		return;
	}
	if(!get_input(command, pt)) {
		suanpan_error("new_paraboliccc() requires a valid tensile yield strength.\n");
		return;
	}
	if(!get_input(command, a)) {
		suanpan_error("new_paraboliccc() requires a valid initial size.\n");
		return;
	}
	if(!get_input(command, a_slope)) {
		suanpan_error("new_paraboliccc() requires a valid hardening slope.\n");
		return;
	}

	auto density = 0.;
	if(!command.eof() && !get_input(command, density)) {
		suanpan_error("new_paraboliccc() requires a valid density.\n");
		return;
	}

	return_obj = make_unique<ParabolicCC>(tag, elastic_modulus, poissons_ratio, beta, m, pt, a, a_slope, density);
}

void new_parallel(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_parallel() requires a valid tag.\n");
		return;
	}

	uword m_tag;
	vector<uword> m_pool;
	while(!command.eof()) if(get_input(command, m_tag)) m_pool.emplace_back(m_tag);

	return_obj = make_unique<Parallel>(tag, uvec(m_pool));
}

void new_sequential(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_sequential() requires a valid tag.\n");
		return;
	}

	uword m_tag;
	vector<uword> m_pool;
	while(!command.eof()) if(get_input(command, m_tag)) m_pool.emplace_back(m_tag);

	if(m_pool.size() == 1) {
		suanpan_error("new_sequential() requires at least two material models.\n");
		return;
	}

	return_obj = make_unique<Sequential>(tag, uvec(m_pool));
}

void new_planestrain(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_planestrain() requires a valid tag.\n");
		return;
	}

	unsigned full_tag;
	if(!get_input(command, full_tag)) {
		suanpan_error("new_planestrain() requires a valid reference material tag.\n");
		return;
	}

	return_obj = make_unique<PlaneStrain>(tag, full_tag);
}

void new_planestress(unique_ptr<Material>& return_obj, istringstream& command) {
	unsigned tag;
	if(!get_input(command, tag)) {
		suanpan_error("new_planestress() requires a valid tag.\n");
		return;
	}

	unsigned full_tag;
	if(!get_input(command, full_tag)) {
		suanpan_error("new_planestress() requires a valid reference material tag.\n");
		return;
	}

	auto max_iteration = 1;
	if(!command.eof() && !get_input(command, max_iteration)) {
		suanpan_error("new_planestress() requires a number for maximum iteration.\n");
		return;
	}

	string use_matrix = "true";
	if(!command.eof() && !get_input(command, use_matrix)) {
		suanpan_error("new_planestress() requires a valid flag to indicate if to use the matrix in iteration.\n");
		return;
	}

	return_obj = make_unique<PlaneStress>(tag, full_tag, max_iteration, is_true(use_matrix));
}
