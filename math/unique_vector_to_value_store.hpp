/*
 * unique_vector_to_value_store.hpp
 *
 *  Created on: Oct 19, 2012
 *      Author: notroot
 */

#ifndef UNIQUE_VECTOR_TO_VALUE_STORE_HPP_
#define UNIQUE_VECTOR_TO_VALUE_STORE_HPP_

#include "unique_vector_to_value_store.h"

#include "math/vector.h"
#include "math/numeric/container_comp.h"
#include "utility/stl_helper_functions.h"

namespace math {

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
unsigned int unique_vector_to_value_store<scalar_type,vector_type,data_type>::size() const {
	return my_map.size();
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
typename unique_vector_to_value_store<scalar_type, vector_type, data_type>::iterator unique_vector_to_value_store<scalar_type,vector_type,data_type>::begin() {
	return my_map.begin();
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
typename unique_vector_to_value_store<scalar_type, vector_type, data_type>::const_iterator unique_vector_to_value_store<scalar_type,vector_type,data_type>::begin() const {
	return my_map.begin();
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
typename unique_vector_to_value_store<scalar_type, vector_type, data_type>::iterator unique_vector_to_value_store<scalar_type,vector_type,data_type>::end() {
	return my_map.end();
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
typename unique_vector_to_value_store<scalar_type, vector_type, data_type>::const_iterator unique_vector_to_value_store<scalar_type,vector_type,data_type>::end() const {
	return my_map.end();
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
const data_type& unique_vector_to_value_store<scalar_type,vector_type,data_type>::get_value(
		const vector_type<scalar_type>& vec) const {
	const_iterator it = my_map.lower_bound(vec); // first element GE
	const_iterator jt = my_map.upper_bound(vec); // first element GT

	if (it == my_map.end() || it == jt) {
		std::stringstream ss;
		ss << vec;
		throw math_exception(
				"unique_vector_to_value_store: Requested value for vector that is not in store: "
						+ ss.str());
	}
	return it->second;
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
typename unique_vector_to_value_store<scalar_type, vector_type, data_type>::iterator unique_vector_to_value_store<scalar_type,vector_type,data_type>::find(
		const vector_type<scalar_type>& vec) {
	iterator it = my_map.lower_bound(vec); // first element GE
	iterator jt = my_map.upper_bound(vec); // first element GT
	if (it == my_map.end() || it == jt) {
		return my_map.end();
	}
//	numeric::lex_comp_less<scalar_type, vector_type> comparator;
//	std::cout << comparator(vec,it->first) << comparator(it->first,vec) << (it == jt) << (jt == my_map.end()) << std::endl;
//	if (jt != my_map.end()) {
//		std::cout << comparator(vec,jt->first) << comparator(jt->first,vec) << std::endl;
//	}
//	std::pair<typename vec_to_val_map::iterator,typename vec_to_val_map::iterator> pr = my_map.equal_range(vec);
//	std::cout << (pr.first != pr.second) << (pr.first == it) << comparator(vec,pr.first->first) << comparator(pr.first->first,vec) << std::endl;

	// the following is a safeguard against a bug (?) which returns it GT but not equal to jt
	if (numeric::lex_comp_less<scalar_type, vector_type>()(vec,it->first)) {
		return my_map.end();
	}
	return it;
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
const data_type& unique_vector_to_value_store<scalar_type,vector_type,data_type>::insert(
		const vector_type<scalar_type>& vec, const data_type& val) {
	iterator it = find(vec);
	if (it == my_map.end()) {
		it = my_map.insert(it, std::make_pair(vec, val));
		//std::cout << std::setprecision(16) << "inserting new " << v << " with val " << val << std::endl;
	}
	return it->second;
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
typename unique_vector_to_value_store<scalar_type, vector_type, data_type>::iterator unique_vector_to_value_store<scalar_type,vector_type,data_type>::insert_missing(
		const vector_type<scalar_type>& vec, const data_type& val) {
	iterator it = my_map.insert(std::make_pair(vec, val)).first;
	return it;
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
template<class comp_op>
const data_type& unique_vector_to_value_store<scalar_type,vector_type,data_type>::update(
		const vector_type<scalar_type>& vec, const data_type& val) {
	iterator it = find(vec);
	if (it == my_map.end()) {
		it = my_map.insert(it, std::make_pair(vec, val));
		//std::cout << std::setprecision(16) << "inserting new " << v << " with val " << val << std::endl;
	} else {
		// replace the stored value with val if it is smaller.
		comp_op comp;
		if (comp(it->second, val)) {
			it->second = val;
		}
	}
	return it->second;
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
void unique_vector_to_value_store<scalar_type,vector_type,data_type>::print(std::ostream& os) const {
	//os << my_map;
	os << "[";
	for (const_iterator it = begin(); it != end(); ++it) {
		if (it != begin())
			os << ",";
		os << it->first << "->" << it->second;
	}
	os << "]";
}

template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
void unique_vector_to_value_store<scalar_type,vector_type,data_type>::clear() {
	my_map.clear();
}

}

#endif /* UNIQUE_VECTOR_TO_VALUE_STORE_HPP_ */
