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
 * @class Orientation
 * @brief A Orientation class.
 * 
 * The Orientation class handles transformations between coordinate systems.
 * Be sure to call set_element_ptr() before updating anything.
 * 
 * @author tlc
 * @date 27/06/2018
 * @version 0.1.0
 * @file Orientation.h
 * @addtogroup Utility
 * @ingroup Element
 * @{
 */

#ifndef ORIENTATION_H
#define ORIENTATION_H

#include <Domain/Tag.h>

class Element;

class Orientation : public Tag {
protected:
	const Element* element_ptr = nullptr;

	vec z_axis;

	double length = 0., inclination = 0.;

	mat direction_cosine;

	void check_element_ptr() const;

	virtual void update_transformation() = 0;
public:
	explicit Orientation(unsigned = 0, vec&& = {});
	Orientation(const Orientation&) = default;           // copy allowed
	Orientation(Orientation&&) = delete;                 // move forbidden
	Orientation& operator=(const Orientation&) = delete; // copy assign forbidden
	Orientation& operator=(Orientation&&) = delete;      // move assign forbidden
	virtual ~Orientation() = default;

	void update_axis(const vec&);

	void set_element_ptr(const Element*);

	double get_length() const;
	double get_inclination() const;
	const mat& get_tranformation() const;

	virtual unique_ptr<Orientation> get_copy() = 0;

	void update_status();

	virtual vec to_local_vec(double) const;
	virtual vec to_global_vec(double) const;
	virtual mat to_global_mass_mat(double) const;
	virtual mat to_global_geometry_mat(double) const;
	virtual mat to_global_stiffness_mat(double) const;

	/**
	 * \brief transform anything from global to local system
	 *        e.g., disp -> disp, vel -> vel, acc -> acc,
	 *        not appilcable to conversion such as disp -> strain
	 * \return variable in local system
	 */
	virtual vec to_local_vec(const vec&) const = 0;
	/**
	 * \brief transform anything from local to global system
	 *        e.g., disp -> disp, vel -> vel, acc -> acc,
	 *        not appilcable to conversion such as disp -> strain
	 * \return variable in global system
	 */
	virtual vec to_global_vec(const vec&) const = 0;
	/**
	 * \brief transform anything from local to global system
	 *        e.g., stiffness -> stiffness.
	 * \return variable in global system
	 */
	virtual mat to_global_mass_mat(const mat&) const;
	virtual mat to_global_geometry_mat(const mat&) const;
	virtual mat to_global_stiffness_mat(const mat&) const = 0;
};

#endif

//! @}
