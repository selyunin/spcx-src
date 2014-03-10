#ifndef GUARD_variable_id_map_h
#define GUARD_variable_id_map_h
/***************************************************************************
 *   Copyright (C) 2005 by Goran Frehse   *
 *   goran.frehse@imag.fr   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

//#include <vector>
//#include <set>
//#include <map>
//#include "bidirectional_map.h"
#include "core/continuous/polyhedra/ppl_polyhedron/bidirectional_map.h"
#include <string>
#include <set>
#include <iostream>
#include "boost/shared_ptr.hpp"
#include "core/continuous/polyhedra/ppl_polyhedron/general.h"
#include "utility/stl_helper_functions.h"
#include "math/vdom/index_to_index_bimap.h"

namespace ppl_polyhedron {

// Create and maintain a bijective function from variable names to variable ids (integers)

//using namespace std;
//using namespace codeproject;

typedef std::string variable_name;
typedef std::set <variable_name> variable_name_set;
typedef unsigned int dimension_t;
typedef size_t variable_ref;
typedef std::set<variable_ref> variable_ref_set;

typedef dimension_t clock_ref;

class variable_id_map;
typedef boost::shared_ptr<variable_id_map> variable_id_map_ptr;

// Bidirectional Map from id to name
// Property 1: The map provides a name for every id from 0 to n, where n is referred to as the dimension.
// If a name is removed, the corresponding id(name) is also removed and the remaining ids>id(name) must
// be "shifted" down by 1, so property 1 is preserved.
// As a consequence of property 1, the dimension always equals the size of mymap.

class variable_id_map {
	//   typedef boost::shared_ptr<variable_id_map> variable_id_map_ptr;
	typedef bidirectional_map<dimension_t,variable_name>::type
			variable_id_map_impl_type;

public:
	typedef variable_id_map_impl_type::const_iterator const_iterator;
	variable_id_map() {
	}
	;
	variable_id_map(dimension_t s); // create default map of size s

	dimension_t size() const;
	dimension_t get_dimension() const;
	void clear() { //mymap.clear();
		mymap=variable_id_map_impl_type();
	}
	;
	bool empty() const;
	variable_ref get_id(const variable_name& s) const;
	const variable_name& get_name(const variable_ref& vr) const;
	bool contains_name(const variable_name& s) const;
	bool contains_id(const variable_ref& vr) const;
	bool contains_names(const variable_name_set& ss) const;
	bool contains_ids(const variable_ref_set& vrs) const;
	bool contains(const variable_id_map& v) const;
	bool operator==(const variable_id_map& v) const;
	bool operator!=(const variable_id_map& v) const {
		return !operator==(v);
	}
	;
	void insert(const variable_ref& vr, const variable_name& s);
	void erase_id(const variable_ref& vr);
	void erase_name(const variable_name& vn);
	void erase_ids(const variable_ref_set& vrs);
	void erase_names(const variable_name_set& vns);
	void union_assign(const variable_id_map& vim);
	void intersection_assign(const variable_id_map& vim);
	void append_to_names(const variable_name& s); // add s to all names
	void
			shifted_union_assign(const variable_id_map& vim,
					const variable_name& s);
	void rename_variable(const variable_name& var1, const variable_name& var2);
	variable_name get_new_name(variable_ref vr);
	void add_space_dimensions(dimension_t dim);
	void map_space_dimensions(const index_to_index_bimap& pfunc);
//	void apply_transform(const variable_transform& transf,
//			bool is_relation=false);
	const_iterator begin() const {
		return mymap.begin();
	}
	;
	const_iterator end() const {
		return mymap.end();
	}
	;

private:
	variable_id_map_impl_type mymap;
};

std::ostream& operator<<(std::ostream& os, const variable_id_map &c);

void get_common_var_names(variable_id_map a1, variable_id_map a2,
		variable_id_map& c_map, index_to_index_bimap& pfunc1, index_to_index_bimap& pfunc2,
		dimension_t& newdim);

void get_from_to_map(variable_id_map from_map, variable_id_map to_map,
		index_to_index_bimap& pfunc);

class clock_ref_set : public std::set<clock_ref> {
	friend std::ostream& operator<<(std::ostream& os, const clock_ref_set &c);

public:
	clock_ref_set() {
	}
	;
	clock_ref_set(clock_ref dim);
	clock_ref_set(clock_ref sdim, clock_ref edim);

	bool contains(const clock_ref& lab) const;

	bool contains(const clock_ref_set& set) const;

	bool equals(const clock_ref_set& set) const;

	void intersection_assign(const clock_ref_set& cset);
	void difference_assign(const clock_ref_set& cset);

	void union_assign(const clock_ref_set& cset);

	clock_ref_set complement(const clock_ref& m, const clock_ref& M);
};

clock_ref_set remap(clock_ref_set crs, const index_to_index_bimap& pfunc);

}

#endif
