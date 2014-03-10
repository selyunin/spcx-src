/*
 * index_to_variable_id_map.hpp
 *
 *  Created on: May 12, 2010
 *      Author: frehse
 */

#ifndef INDEX_TO_VARIABLE_ID_MAP_HPP_
#define INDEX_TO_VARIABLE_ID_MAP_HPP_

#include <stdexcept>
#include <sstream>
#include "utility/stl_helper_functions.h"
#include "utility/basic_exception.h"

inline
index_to_variable_id_map::size_type index_to_variable_id_map::dimensions() const {
	return my_id_vector.size();
}

inline
const variable_id& index_to_variable_id_map::get_id(const index_type& i) const {
	if (i >= dimensions()) {
		std::stringstream ss;
		print(ss);
		throw std::out_of_range("index " + int2string(i)
				+ " not in index_to_variable_id_map "+ss.str());
	}
	return my_id_vector[i];
}

inline
const variable_id_set& index_to_variable_id_map::get_ids() const {
	return my_ids;
}

inline
const std::vector<variable_id>& index_to_variable_id_map::get_id_vector() const {
	return my_id_vector;
}

inline
const index_type& index_to_variable_id_map::get_index(const variable_id& id) const {
	std::map<variable_id, index_type>::const_iterator pos =
			my_id_to_index_map.find(id);
	if (pos == my_id_to_index_map.end()) {
		throw basic_exception("id " + int2string(id)
				+ " not in index_to_variable_id_map, name is "+variable::get_name(id)+".");
	}
	return pos->second;
}

inline
bool index_to_variable_id_map::has_id(const variable_id& id) const {
	std::map<variable_id, index_type>::const_iterator pos =
			my_id_to_index_map.find(id);
	return !(pos == my_id_to_index_map.end());
}

inline
index_type index_to_variable_id_map::check_for_index(const variable_id& id,
		bool& has_id) const {
	std::map<variable_id, index_type>::const_iterator pos =
			my_id_to_index_map.find(id);
	has_id = !(pos == my_id_to_index_map.end());
	if (has_id)
		return pos->second;
	else
		return 0;
}

inline
const index_to_variable_id_map_ptr& index_to_variable_id_map::empty_map() {
	static index_to_variable_id_map_ptr empty_map_ptr =
			get_index_to_variable_id_map_ptr(std::vector<variable_id>());
	return empty_map_ptr;
}

#endif /* INDEX_TO_VARIABLE_ID_MAP_HPP_ */
