/*
 * unique_scalar_to_value_store.h
 *
 *  Created on: Jan 20, 2011
 *      Author: frehse
 */

#ifndef UNIQUE_SCALAR_TO_VALUE_STORE_H_
#define UNIQUE_SCALAR_TO_VALUE_STORE_H_

#include <map>
#include "math/numeric/comp.h"
#include "utility/stl_helper_functions.h"

namespace math {

/** A class for storing unique scalars together with a value, with fast search access.
 *
 * Stored scalars are unique up to numerical approximation. The math::numeric comparator
 * is used to tell when two vectors are equivalent and when not.
 */
template<typename scalar_type, typename value_type>
class unique_scalar_to_value_store {
public:
	typedef std::map<scalar_type, value_type,
			math::numeric::comp_less<scalar_type> >
			vec_to_val_map;

	typedef typename vec_to_val_map::iterator iterator;
	typedef typename vec_to_val_map::const_iterator const_iterator;

	/** The iterator pointing to the first element. */
	iterator begin() {
		return my_map.begin();
	}
	;

	/** The const iterator pointing to the first element. */
	const_iterator begin() const {
		return my_map.begin();
	}
	;

	/** The iterator pointing to the last element. */
	iterator end() {
		return my_map.end();
	}
	;

	/** The iterator pointing to the last element. */
	const_iterator end() const {
		return my_map.end();
	}
	;

	/** Return the value stored for scalar vec.
	 *
	 * Throws if the vector is not in the store. */
	const value_type& get_value(const scalar_type& vec) const {
		const_iterator it = my_map.lower_bound(vec); // first element GE
		const_iterator jt = my_map.upper_bound(vec); // first element GT
		if (it == my_map.end() || it == jt) {
			std::stringstream ss;
			ss << vec;
			throw math_exception(
					"unique_scalar_to_value_store: Requested value for scalar that is not in store: "
							+ ss.str());
		}
		return it->second;
	}
	;

	/** Find the equivalent scalar if already in the store.
	 *
	 * Returns an iterator to the equivalent scalar if there
	 * is one, and end() otherwise. */
	const_iterator find(const scalar_type& vec) const {
		const_iterator it = my_map.lower_bound(vec); // first element GE
		const_iterator jt = my_map.upper_bound(vec); // first element GT
		if (it == my_map.end() || it == jt) {
			return my_map.end();
		}
		return it;
	}
	;

	/** Find the equivalent scalar if already in the store.
	 *
	 * Returns an iterator to the equivalent scalar if there
	 * is one, and end() otherwise. */
	iterator find(const scalar_type& vec) {
		iterator it = my_map.lower_bound(vec); // first element GE
		iterator jt = my_map.upper_bound(vec); // first element GT
		if (it == my_map.end() || it == jt) {
			return my_map.end();
		}
		return it;
	}
	;

	/** Insert a scalar only if there is no equivalent already in the store.
	 * If there is, return stored the value. */
	const value_type& insert(const scalar_type& vec,
			const value_type& val = value_type()) {
		iterator it = find(vec);
		if (it == my_map.end()) {
			it = my_map.insert(it, std::make_pair(vec, val));
			//std::cout << std::setprecision(16) << "inserting new " << v << " with val " << val << std::endl;
		}
		return it->second;
	}
	;

	/** Insert a scalar and value that are known not to be in the store.
	 *
	 * This variant is faster, since it skips the check whether vec is
	 * already in the store.
	 *
	 * Returns an iterator to the inserted element.
	 * */
	iterator insert_missing(const scalar_type& vec,
			const value_type& val) {
		iterator it = my_map.insert(std::make_pair(vec, val)).first;
		return it;
	}
	;

	/** Replace any stored value for vec (or equivalent) by val if it is larger.
	 *
	 * If vec is not yet stored, insert it together with val.
	 * Uses the operator comp_op to compare values.
	 *
	 * For example, if less is used, the stored value for each vector
	 * is the maximum. Using a different operator, other bounds can be obtained (min, etc.).
	 *
	 * The comparison operator can be defined just as STL comparison
	 * operators for sets, maps, etc. */
	template<class comp_op>
	const value_type& update(const scalar_type& vec,
			const value_type& val) {
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
	;

	/** Output to stream in human-readable form. */
	void print(std::ostream& os) const {
		//os << my_map;
		os << "[";
		for (const_iterator it = begin(); it != end(); ++it) {
			if (it != begin())
				os << ",";
			os << it->first << "->" << it->second;
		}
		os << "]";
	}
	;

	/** Clear store. */
	void clear() {
		my_map.clear();
	}
	;

private:
	vec_to_val_map my_map;
};

}


#endif /* UNIQUE_SCALAR_TO_VALUE_STORE_H_ */
