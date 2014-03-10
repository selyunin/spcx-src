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

#include "core/continuous/polyhedra/ppl_polyhedron/variable_id_map.h"

namespace ppl_polyhedron {

using namespace std;

variable_id_map::variable_id_map(dimension_t s) // create default map of size s
{
	for (dimension_t i=0; i<s; ++i)
		//         mymap[i]="x"+int2string(i);
		mymap.insert(variable_id_map_impl_type::value_type(i, "x"+int2string(i)));
}

dimension_t variable_id_map::size() const {
	return mymap.size();
}

dimension_t variable_id_map::get_dimension() const {
	// to do: this is not consistent, needs fixing
	return mymap.size();
}

bool variable_id_map::empty() const {
	return mymap.empty();
}

variable_ref variable_id_map::get_id(const variable_name& s) const {
	//      return mymap.to[s];
	//const_iterator it=get<to>(mymap).find(s);
	bidirectional_map<dimension_t,variable_name>::type::index<to>::type::const_iterator
			it = mymap.get<to>().find(s);
	if (it!=mymap.get<to>().end()) {
		return it->first;
	} else {
		throw_error("variable_id_map::get_id: Variable " + s + " not found");
		return variable_ref();
	}
}

const variable_name& variable_id_map::get_name(const variable_ref& vr) const {
	//return mymap.from[vr];
	/*cout << "hello" << flush;
	 variable_id_map_impl_type tmp=mymap;
	 cout << tmp.from[vr] << flush;      */
	//variable_id_map_impl_type::index<from>::type::const_iterator it = mymap.get<from>().find(vr);
	bidirectional_map<dimension_t,variable_name>::type::index<from>::type::const_iterator
			it = mymap.get<from>().find(vr);
	if (it!=mymap.get<from>().end()) {
		return it->second;
	} else {
		throw_error("variable_id_map::get_name: Variable " + int2string(vr)
				+ " not found");
		return mymap.get<from>().begin()->second;
	}
}

bool variable_id_map::contains_name(const variable_name& s) const {
	return mymap.get<to>().find(s) != mymap.get<to>().end();
}

bool variable_id_map::contains_id(const variable_ref& vr) const {
	return mymap.get<from>().find(vr) != mymap.get<from>().end();
}

bool variable_id_map::contains_names(const variable_name_set& ss) const {
	for (variable_name_set::const_iterator it=ss.begin(); it!=ss.end(); ++it) {
		if (!contains_name(*it))
			return false;
	}
	return true;
}

bool variable_id_map::contains_ids(const variable_ref_set& vrs) const {
	for (variable_ref_set::const_iterator it=vrs.begin(); it!=vrs.end(); ++it) {
		if (!contains_id(*it))
			return false;
	}
	return true;
}

bool variable_id_map::contains(const variable_id_map& v) const {
	variable_id_map_impl_type::index<from>::type::const_iterator k=
			mymap.get<from>().begin();
	for (variable_id_map_impl_type::index<from>::type::const_iterator it=
			v.mymap.get<from>().begin(); it!=v.mymap.get<from>().end(); ++it) {
		k=mymap.get<from>().find(it->first);
		if (k!=mymap.get<from>().end()) {
			if (k->second!=it->second)
				return false;
		} else
			return false;
		//         if ((mymap.from[it->first].get())!=(it->second))
		//            return false;
	}
	return true;
}

bool variable_id_map::operator==(const variable_id_map& v) const
{
	return contains(v) && v.contains(*this);
}

void variable_id_map::insert(const variable_ref& vr, const variable_name& s) {
	// check for consistency is done in bidirectional_map
	//      mymap[vr]=s;
	mymap.insert(variable_id_map_impl_type::value_type(vr, s));
}

void variable_id_map::erase_id(const variable_ref& vr) {
	if (contains_id(vr)) {
		/*         mymap.erase(vr);
		 // decrease the ids of the following ids by 1
		 for (bidirectional_map<dimension_t,variable_name>::reverse_iterator it=mymap.rbegin();it!=mymap.rend();++it)
		 {
		 if (it->first > vr)
		 mymap[it->second]=it->first-1;
		 }*/
		variable_id_map_impl_type newmap;
		for (variable_id_map_impl_type::iterator it=mymap.begin(); it
				!=mymap.end(); ++it) {
			if (it->first> vr) {
				//               newmap[it->second]=it->first-1;
				newmap.insert(variable_id_map_impl_type::value_type(
								it->first-1, it->second));
			} else if (it->first < vr) {
				//               newmap[it->second]=it->first;
				newmap.insert(variable_id_map_impl_type::value_type(it->first,
								it->second));
			}
		}
		mymap.swap(newmap);
	}
}

void variable_id_map::erase_name(const variable_name& vn) {
	if (contains_name(vn))
	erase_id(get_id(vn));
}

void variable_id_map::erase_ids(const variable_ref_set& vrs) {
	// must start with the highest id, otherwise id shifting leads to false results
	for (variable_ref_set::const_reverse_iterator it=vrs.rbegin(); it
			!=vrs.rend(); ++it) {
		erase_id(*it);
	}
}

void variable_id_map::erase_names(const variable_name_set& vns) {
	for (variable_name_set::const_iterator it=vns.begin(); it!=vns.end(); ++it) {
		erase_name(*it);
	}
}

void variable_id_map::union_assign(const variable_id_map& vim) {
	// to do: check whether identical entries are accepted or not
	mymap.insert(vim.begin(), vim.end());
}

void variable_id_map::intersection_assign(const variable_id_map& vim) {
	// to do: this is not a pretty implementation
	// erase all those of mymap that are not in vim
	variable_id_map_impl_type::iterator it=mymap.begin();
	while (it!=mymap.end()) {
		if (vim.contains_id(it->first) && vim.contains_name(it->second)
				&& vim.get_name(it->first)==it->second) {
			++it;
		} else {
			mymap.erase(it);
		}
	}

	// to do: shift everything down to zero
}

void variable_id_map::append_to_names(const variable_name& s) // add s to all names
{
	// to do: for now, don't mess with changing values directly, since this might temporarily violate uniqueness of values
	// instead, recontruct the map
	variable_id_map_impl_type newmap;
	for (variable_id_map_impl_type::const_iterator it=mymap.begin(); it
			!=mymap.end(); ++it) {
		//         newmap[it->first]=it->second+s;
		newmap.insert(variable_id_map_impl_type::value_type(it->first,
						it->second+s));
	}
	mymap.swap(newmap);
}

void variable_id_map::shifted_union_assign(const variable_id_map& vim,
		const variable_name& s) {
	// add the shifted entries of vim, extended with string(s)
	dimension_t mysize=size();
	for (variable_id_map_impl_type::const_iterator it=vim.mymap.begin(); it
			!=vim.mymap.end(); ++it) {
		//         mymap[it->first+mysize]=it->second+s;
		mymap.insert(variable_id_map_impl_type::value_type(it->first+mysize,
						it->second+s));
	}
}

void variable_id_map::rename_variable(const variable_name& var1,
		const variable_name& var2) {
	if (var1!=var2) {
		if (contains_name(var1)) {
			variable_ref ref=get_id(var1);
			// don't use variable_id_map::erase because the ids don't need to be shifted -- the same id is used to define (ref,var2)
			mymap.get<to>().erase(var1);
			//            mymap[ ref ]=var2;
			mymap.insert(variable_id_map_impl_type::value_type(ref, var2));
		} else
		throw_error("variable_id_map::rename_variable : variable " +var1
				+" not found");
	}
}

variable_name variable_id_map::get_new_name(variable_ref vr) {
	variable_name s="x";
	while (mymap.get<to>().find(s+int2string(vr))!=mymap.get<to>().end()) {
		++vr;
	}
	return s+int2string(vr);
}

void variable_id_map::add_space_dimensions(dimension_t dim) {
	dimension_t old_dim=get_dimension();
	for (dimension_t i=old_dim; i<dim; ++i) {
		mymap.insert(variable_id_map_impl_type::value_type(i, get_new_name(i)));
	}
}

void variable_id_map::map_space_dimensions(const index_to_index_bimap& pfunc) {
	variable_id_map_impl_type newmap;
	variable_ref y;
	for (variable_id_map_impl_type::const_iterator it=mymap.begin(); it
			!=mymap.end(); ++it) {
		if (pfunc.maps(it->first, y))
		newmap.insert(variable_id_map_impl_type::value_type(y, it->second));
	}
	mymap.swap(newmap);
}

//void variable_id_map::apply_transform(const variable_transform& transf,
//		bool is_relation) {
//	variable_id_map_impl_type newmap;
//	variable_ref y;
//	for (variable_id_map_impl_type::const_iterator it=mymap.begin(); it
//			!=mymap.end(); ++it) {
//		if (transf.transform_variable(it->first, y, is_relation))
//		newmap.insert(variable_id_map_impl_type::value_type(y, it->second));
//	}
//	// add names for the new variables
//	clock_ref_set crs=transf.get_new_variables(size(), is_relation);
//	for (clock_ref_set::const_iterator it=crs.begin(); it!=crs.end(); ++it) {
//		newmap.insert(variable_id_map_impl_type::value_type(*it,
//						get_new_name(*it)));
//	}
//	mymap.swap(newmap);
//}

std::ostream& operator<<(std::ostream& os, const variable_id_map &c) {
	os << "[";
	for (dimension_t i = 0; i<c.size(); ++i) {
		if (i!=0)
		os << ",";
		os << i << ":" << c.get_name(i);
	};
	os << "]";
	return os;
}

void get_common_var_names(variable_id_map a1, variable_id_map a2,
		variable_id_map& c_map, index_to_index_bimap& pfunc1,
		index_to_index_bimap& pfunc2, dimension_t& newdim) {
	// Builds a common map of variable names for the composition C=A1||A2 of two automata or state sets
	// Input: a1_map, a2_map
	// Output: c_map, pfunc1, pfunc2, newdim
	// The variables of A1 remain in the same place, while the variables of A2 must
	// be remapped according to pfunc2. pfunc1 is supplied for completeness.
	// The new dimension (number of variables) is given in newdim.
	//
	// Note: pfunc2 refers to a variable set of dimension newdim, because the mapping (e.g. in the PPL)
	// may not be capable of augmenting the dimension. Thus we assume that pfunc is a mapping from newdim to newdim,
	// i.e., that the set is expanded to dimension newdim before being mapped.

	// Note: spec_aut = A2, caut = C

	// Todo: Assert that every number in [0,a1_map.size()-1] occurs in the codomain of a1_map,
	//       and analogously for a2map
	// Initialize

	c_map=a1; // the variables of A1 stay the same
	newdim=0;
	for (variable_id_map::const_iterator i=a1.begin(); i!=a1.end(); ++i)
	newdim=max(newdim, i->first);

	index_to_index_bimap newfunc;
	pfunc2=newfunc;
	newdim=a1.size();
	string var_name;

	// Add the variables of A2
	// proceed by order of the variables -> look for each name and add its reference to the new map

	variable_id_map::const_iterator kc;
	for (kc=a2.begin(); kc!=a2.end(); ++kc) {
		if (c_map.contains_name(kc->second)) // name known, insert mapping in pfunc
		{
			pfunc2.insert(kc->first, c_map.get_id(kc->second));
		} else // new variable, add to dimension of a1
		{
			c_map.insert(newdim, kc->second);
			pfunc2.insert(kc->first, newdim);
			++newdim;
		}
	}
	// all variables from 0 to newdim need to be in pfunc, so add
	pfunc2.fill_up_to(newdim);

	pfunc1=newfunc;
	pfunc1.fill_up_to(newdim); // variables stay in their place
	//cout << pfunc;
}

ostream& operator<<(ostream& os, const clock_ref_set &c) {
	os << "[";
	for (clock_ref_set::const_iterator j=c.begin(); j!=c.end(); ++j) {
		if (j!=c.begin())
		os << "," << *j;
		else
		os << *j;
	};
	os << "]";
	return os;
}

clock_ref_set::clock_ref_set(clock_ref dim) {
	for (clock_ref i=0; i<dim; ++i)
	insert(i);
}

clock_ref_set::clock_ref_set(clock_ref sdim, clock_ref edim) {
	// create dimensions sdim ... edim-1
	for (clock_ref i=sdim; i<edim; ++i)
	insert(i);
}

bool clock_ref_set::contains(const clock_ref& lab) const {
	return (this->find(lab)!=end());
}

bool clock_ref_set::contains(const clock_ref_set& set) const {
	bool cont=true;

	clock_ref_set::const_iterator i;
	for (i=set.begin(); i!=set.end(); ++i) {
		if (!contains(*i))
		cont=false;
	};
	return cont;
}

bool clock_ref_set::equals(const clock_ref_set& set) const {
	return contains(set) && set.contains(*this);
}

void clock_ref_set::difference_assign(const clock_ref_set& cset) {
	clock_ref_set::iterator i=begin();
	while (i!=end()) {
		if (cset.contains(*i)) {
			erase(i);
			i=begin();
		} else
		++i;
	};
}

void clock_ref_set::intersection_assign(const clock_ref_set& cset) {
	clock_ref_set::iterator i=begin();
	while (i!=end()) {
		if (cset.contains(*i)) {
			++i;
		} else {
			erase(i);
			i=begin();
		}
	}
}

void clock_ref_set::union_assign(const clock_ref_set& cset) {
	insert(cset.begin(), cset.end());
}

clock_ref_set clock_ref_set::complement(const clock_ref& m, const clock_ref& M) {
	clock_ref_set crs;
	for (clock_ref i=m; i!=M; ++i)
	crs.insert(i);
	crs.difference_assign(*this);
	return crs;
}

clock_ref_set remap(clock_ref_set crs, const index_to_index_bimap& pfunc) {
	// remap the elements of crs according to pfunc
	clock_ref_set newcrs;
	for (clock_ref_set::const_iterator it=crs.begin(); it!=crs.end(); ++it) {
		if (pfunc.in_domain(*it))
		newcrs.insert(pfunc.get_map(*it));
	};

	return newcrs;
}
;



void get_from_to_map(variable_id_map from_map, variable_id_map to_map,
		index_to_index_bimap& pfunc) {
	// map the variables of from_map to to_map
	// the variables that are not found are projected away
	// note that remapping according to pfunc, the dimension can be smaller than that of to_map,
	// so that adding space dimensions is warranted.
	pfunc=index_to_index_bimap();
	variable_id_map::const_iterator it=from_map.begin();
	while (it!=from_map.end()) {
		if (to_map.contains_name(it->second)) {
			pfunc.insert(it->first, to_map.get_id(it->second));
		}
		++it;
	}
}

// ------------------------------------------------------



}

