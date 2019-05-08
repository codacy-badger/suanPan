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
 * @class Element
 * @brief A Element class.
 * @author tlc
 * @date 24/05/2018
 * @version 0.2.0
 * @file Element.h
 * @addtogroup Element
 * @{
 */

#ifndef ELEMENT_H
#define ELEMENT_H

#include <Domain/Tag.h>
#include <Element/Visualisation/vtkBase.h>

class Node;
class DomainBase;
class Material;
class Section;
enum class OutputType;

using std::array;
using std::vector;

struct ElementData {
	const uvec node_encoding; // node encoding
	const uvec material_tag;  // material tags
	const uvec section_tag;   // section tags

	const bool nlgeom = false; // nonlinear geometry switch

	bool update_mass = true;      // flag to indicate if update matrix
	bool update_damping = true;   // flag to indicate if update matrix
	bool update_stiffness = true; // flag to indicate if update matrix
	bool update_geometry = true;  // flag to indicate if update matrix

	uvec dof_encoding; // DoF encoding vector

	mat initial_mass;      // mass matrix
	mat initial_damping;   // damping matrix
	mat initial_stiffness; // stiffness matrix
	mat initial_geometry;  // geometry matrix

	mat trial_mass;      // mass matrix
	mat trial_damping;   // damping matrix
	mat trial_stiffness; // stiffness matrix
	mat trial_geometry;  // geometry matrix

	mat current_mass;      // mass matrix
	mat current_damping;   // damping matrix
	mat current_stiffness; // stiffness matrix
	mat current_geometry;  // geometry matrix

	vec trial_resistance;   // resistance vector
	vec current_resistance; // resistance vector
	vec trial_body_force;
	vec current_body_force;
	vec trial_traction;
	vec current_traction;
};

class Element : protected ElementData, public Tag, public vtkBase {
	const unsigned num_node;                      // number of nodes
	const unsigned num_dof;                       // number of DoFs
	const unsigned num_size = num_dof * num_node; // number of size

	const bool initialized = false;
	const bool symmetric = false;
protected:
	vector<weak_ptr<Node>> node_ptr; // node pointers

	virtual mat get_coordinate(unsigned) const;

	virtual vec get_incre_displacement() const;
	virtual vec get_incre_velocity() const;
	virtual vec get_incre_acceleration() const;
	virtual vec get_trial_displacement() const;
	virtual vec get_trial_velocity() const;
	virtual vec get_trial_acceleration() const;
	virtual vec get_current_displacement() const;
	virtual vec get_current_velocity() const;
	virtual vec get_current_acceleration() const;

	virtual vec get_node_incre_resistance() const;
	virtual vec get_node_trial_resistance() const;
	virtual vec get_node_current_resistance() const;

	friend mat get_coordinate(const Element*, unsigned);

	friend vec get_incre_displacement(const Element*);
	friend vec get_incre_velocity(const Element*);
	friend vec get_incre_acceleration(const Element*);
	friend vec get_trial_displacement(const Element*);
	friend vec get_trial_velocity(const Element*);
	friend vec get_trial_acceleration(const Element*);

	friend vec get_node_incre_resistance(const Element*);
	friend vec get_node_trial_resistance(const Element*);

	vector<shared_ptr<Material>> get_material(const shared_ptr<DomainBase>&) const;
	vector<shared_ptr<Section>> get_section(const shared_ptr<DomainBase>&) const;
public:
	explicit Element(unsigned = 0, // tag
	                 unsigned = 0, // number of nodes
	                 unsigned = 0, // number of dofs
	                 uvec&& = {},  // node encoding
	                 uvec&& = {},  // material tags
	                 uvec&& = {},  // section tags
	                 bool = false  // nonlinear geometry switch
	);
	Element(const Element&) = delete;            // copy forbidden
	Element(Element&&) = delete;                 // move forbidden
	Element& operator=(const Element&) = delete; // assign forbidden
	Element& operator=(Element&&) = delete;      // assign forbidden

	~Element() override;

	virtual void initialize(const shared_ptr<DomainBase>&) = 0;

	void set_initialized(bool) const;
	void set_symmetric(bool) const;
	bool is_initialized() const;
	bool is_symmetric() const;
	bool is_nlgeom() const;

	void update_dof_encoding();

	const uvec& get_dof_encoding() const;
	const uvec& get_node_encoding() const;

	unsigned get_dof_number() const;
	unsigned get_node_number() const;
	unsigned get_total_number() const;

	const vector<weak_ptr<Node>>& get_node_ptr() const;

	virtual const vec& get_trial_resistance() const;
	virtual const vec& get_current_resistance() const;

	virtual const mat& get_trial_mass() const;
	virtual const mat& get_trial_damping() const;
	virtual const mat& get_trial_stiffness() const;
	virtual const mat& get_trial_geometry() const;
	virtual const mat& get_trial_secant() const;

	virtual const mat& get_current_mass() const;
	virtual const mat& get_current_damping() const;
	virtual const mat& get_current_stiffness() const;
	virtual const mat& get_current_geometry() const;
	virtual const mat& get_current_secant() const;

	virtual const mat& get_initial_mass() const;
	virtual const mat& get_initial_damping() const;
	virtual const mat& get_initial_stiffness() const;
	virtual const mat& get_initial_geometry() const;
	virtual const mat& get_initial_secant() const;

	virtual int update_status() = 0;
	virtual int clear_status() = 0;
	virtual int commit_status() = 0;
	virtual int reset_status() = 0;

	virtual vector<vec> record(OutputType);
};

#endif

//! @}
