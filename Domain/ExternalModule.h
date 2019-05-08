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
 * @class ExternalModule
 * @brief A ExternalModule class handles communication between the main program
 * and external library.
 *
 * @author tlc
 * @date 28/09/2017
 * @version 0.1.1
 * @file ExternalModule.h
 * @addtogroup Utility
 * @{
 */

#ifndef EXTERNALMODULE_H
#define EXTERNALMODULE_H

#include <suanPan.h>
#include <Domain/DomainBase.h>
#include <Toolbox/utility.h>

class Element;
class Load;
class Material;
class Section;
class Solver;
class Amplitude;
class Modifier;

class ExternalModule {
	void* ext_library = nullptr;
	void* ext_creator = nullptr;
public:
	const string library_name;

	explicit ExternalModule(string);
	ExternalModule(const ExternalModule&) = delete;
	ExternalModule(ExternalModule&&) = delete;
	ExternalModule& operator=(const ExternalModule&) = delete;
	ExternalModule& operator=(ExternalModule&&) = delete;
	~ExternalModule();

	bool locate_module(string);

	void new_object(unique_ptr<Element>&, istringstream&) const;
	void new_object(unique_ptr<Load>&, istringstream&) const;
	void new_object(unique_ptr<Material>&, istringstream&) const;
	void new_object(unique_ptr<Section>&, istringstream&) const;
	void new_object(unique_ptr<Solver>&, istringstream&) const;
	void new_object(unique_ptr<Amplitude>&, istringstream&) const;
	void new_object(unique_ptr<Modifier>&, istringstream&) const;
};

class load {
public:
	template<typename T> static void object(unique_ptr<T>&, const shared_ptr<DomainBase>&, const string&, istringstream&);
};

template<typename T> void load::object(unique_ptr<T>& new_object, const shared_ptr<DomainBase>& domain, const string& id, istringstream& command) {
	// check if the library is already loaded
	auto code = false;
	for(const auto& I : domain->get_external_module_pool())
		if(is_equal(I->library_name, id) || I->locate_module(id)) {
			code = true;
			break;
		}

	// not loaded then try load it
	if(!code && domain->insert(make_shared<ExternalModule>(id))) code = true;

	// if loaded find corresponding function
	if(code)
		for(const auto& I : domain->get_external_module_pool()) {
			if(I->locate_module(id)) I->new_object(new_object, command);
			if(new_object != nullptr) break;
		}
}

#endif
