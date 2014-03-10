/*
 * unique_vector_to_value_store.h
 *
 *  Created on: Oct 25, 2010
 *      Author: frehse
 */

#ifndef UNIQUE_VECTOR_TO_VALUE_STORE_H_
#define UNIQUE_VECTOR_TO_VALUE_STORE_H_

#include <map>

/** Forward declarations */
namespace math {
namespace numeric {
template<typename scalar_type, template<typename > class container_type>
class lex_comp_less;
}
}

namespace math {

/** A class for storing unique vectors together with a value, with fast search access.
 *
 * Vectors are unique up to numerical approximation. The math::numeric comparator
 * is used to tell when two vectors are equivalent and when not.
 */
template<typename scalar_type, template<typename > class vector_type,
		typename data_type>
class unique_vector_to_value_store {
public:
	typedef std::map<vector_type<scalar_type> , data_type,
			math::numeric::lex_comp_less<scalar_type, vector_type> >
			vec_to_val_map;

	typedef typename vec_to_val_map::iterator iterator;
	typedef typename vec_to_val_map::const_iterator const_iterator;
	typedef typename vec_to_val_map::value_type value_type;

	/** Returns the size of the store. */
	unsigned int size() const;

	/** The iterator pointing to the first element. */
	iterator begin();

	/** The const iterator pointing to the first element. */
	const_iterator begin() const;

	/** The iterator pointing to the last element. */
	iterator end();

	/** The iterator pointing to the last element. */
	const_iterator end() const;

	/** Return the value stored for vector vec.
	 *
	 * Throws if the vector is not in the store. */
	const data_type& get_value(const vector_type<scalar_type>& vec) const;

	/** Find the equivalent vector if already in the store.
	 *
	 * Returns an iterator to the equivalent vector if there
	 * is one, and end() otherwise. */
	iterator find(const vector_type<scalar_type>& vec);

	/** Insert a vector only if there is no equivalent already in the store.
	 *
	 * If there is, return the stored value. */
	const data_type& insert(const vector_type<scalar_type>& vec,
			const data_type& val = data_type());

	/** Insert a vector and value that are known not to be in the store.
	 *
	 * This variant is faster, since it skips the check whether vec is
	 * already in the store.
	 *
	 * Returns an iterator to the inserted element.
	 * */
	iterator insert_missing(const vector_type<scalar_type>& vec,
			const data_type& val);

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
	const data_type& update(const vector_type<scalar_type>& vec,
			const data_type& val);

	/** Output to stream in humain-readable form. */
	void print(std::ostream& os) const;

	/** Clear store. */
	void clear();

private:
	vec_to_val_map my_map;
};

}

#include "unique_vector_to_value_store.hpp"

#endif /* UNIQUE_VECTOR_TO_VALUE_STORE_H_ */
