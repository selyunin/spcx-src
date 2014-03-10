/*
 * positional_vdomain.cpp
 *
 *  Created on: Mar 25, 2010
 *      Author: frehse
 */

#include "positional_vdomain.h"

#include "index_to_index_bimap.h"

positional_vdomain::positional_vdomain(const std::vector<
		variable>& vars) {
	std::vector<variable_id> new_v(vars.size());
	for (std::vector<variable_id>::size_type i = 0; i < new_v.size(); ++i) {
		new_v[i] = vars[i].get_id();
	}
	my_iimap
			= index_to_variable_id_map::get_index_to_variable_id_map_ptr(new_v);
}

positional_vdomain::positional_vdomain(const variable_id_list& vars) {
	std::vector<variable_id> new_v(vars.size());
	variable_id_list::const_iterator it=vars.begin();
	for (std::vector<variable_id>::size_type i = 0; i < new_v.size(); ++i) {
		new_v[i] = *it;
		++it;
	}
	my_iimap
			= index_to_variable_id_map::get_index_to_variable_id_map_ptr(new_v);
}

void positional_vdomain::add_variables(const std::set<variable>& x) {
	variable_id_set vis;
	for (std::set<variable>::const_iterator it = x.begin(); it != x.end(); ++it) {
		vis.insert(it->get_id());
	}
	my_iimap = my_iimap->index_to_variable_id_map::get_map_with_ids_added(vis);
}

void positional_vdomain::remove_variable(const variable& x) {
	variable_id_set vis;
	vis.insert(x.get_id());
	my_iimap
			= my_iimap->index_to_variable_id_map::get_map_with_ids_removed(vis);
}

void positional_vdomain::remove_variables(
		const std::set<variable>& x) {
	variable_id_set vis;
	for (std::set<variable>::const_iterator it = x.begin(); it != x.end(); ++it) {
		vis.insert(it->get_id());
	}
	my_iimap
			= my_iimap->index_to_variable_id_map::get_map_with_ids_removed(vis);
}

std::set<variable> positional_vdomain::get_variables() const {
	std::set<variable> res;
	const variable_id_set& vis = my_iimap->get_ids();
	for (variable_id_set::const_iterator it = vis.begin(); it != vis.end(); ++it) {
		res.insert(variable(*it));
	}
	return res;
}

std::vector<variable> positional_vdomain::get_variable_vector() const {
	std::vector<variable> res;
	const std::vector<variable_id>& ids = my_iimap->get_id_vector();
	for (std::vector<variable_id>::size_type i=0; i < ids.size(); ++i) {
		res[i]=variable(ids[i]);
	}
	return res;
}

void positional_vdomain::print(std::ostream& os) const {
	os << "[";
	for (size_type i = 0; i < size(); ++i) {
		if (i > 0)
			os << ",";
		os << variable(my_iimap->get_id(i));
	}
	os << "]";
}

position_map compute_reordering(const positional_vdomain& D1,
		const positional_vdomain& D2, bool use_invalid_pos) {
	position_map f;
	positional_vdomain::size_type j;
	for (size_t i = 0; i < D1.size(); ++i) {
		bool found = D2.in_domain(D1.get_variable(i),j);
		if (!found) {
			j = positional_vdomain::invalid_pos();
		}
		f.insert(i,j);
	}
	return f;
}

positional_vdomain compose(const positional_vdomain& D1,
		const positional_vdomain& D2, position_map& f1, position_map& f2) {
	positional_vdomain::size_type newdim;
	index_to_variable_id_map_ptr M;
	get_common_map(D1.my_iimap, D2.my_iimap, M, newdim, f2);
	positional_vdomain D(M);
	f1 = position_map::identity(D1.size());
	return D;
}

positional_vdomain compose(const positional_vdomain& D1,
		const positional_vdomain& D2) {
	position_map f1;
	position_map f2;
	return compose(D1,D2,f1,f2);
}


std::ostream& operator<<(std::ostream& os, const positional_vdomain& D) {
	D.print(os);
	return os;
}

std::set<variable> create_variable_set(const variable_id_set& vis) {
	std::set<variable> res;
	for (variable_id_set::const_iterator it = vis.begin(); it != vis.end(); ++it) {
		res.insert(variable(*it));
	}
	return res;
}

variable_id_set create_variable_id_set(const std::set<variable>& vars) {
	variable_id_set vis;
	for (std::set<variable>::const_iterator it = vars.begin(); it != vars.end(); ++it) {
		vis.insert(it->get_id());
	}
	return vis;
}


