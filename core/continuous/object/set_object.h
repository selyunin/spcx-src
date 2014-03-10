/*
 * set_object.h
 *
 *  Created on: Mar 25, 2010
 *      Author: frehse
 */

#ifndef SET_OBJECT_H_
#define SET_OBJECT_H_

#include <iostream>
#include "core/continuous/continuous_set.h"
#include "math/vdom/positional_vdomain.h"

typedef positional_vdomain positional_domain;

namespace continuous {
namespace object {

typedef global_types::rational_type rational_type;
typedef global_types::float_type float_type;

/** This class is intended to provide a user-friendly interface
 * for creating and manipulating sets and their support function
 * representations.
 * Sets can be defined over rational or floating point numbers.
 *
 * It wraps corresponding calls to basic and derived classes of
 * the continuous_set hierarchy.
 *
 * @attention The class uses shallow copy (like Java).
 * To make a deep copy of an object, use clone().
 */

class set_object {
public:
	enum number_type {
		RATIONAL = global_types::default_rational_type_id,
		FLOAT = global_types::default_float_type_id
		//		RATIONAL = global_types::type_identifier<rational_type>::coeff,
		//		FLOAT = global_types::type_identifier<float_type>::coeff
	};

	set_object();

	/** Create a continuous set from a string using.
	 *
	 * Accepts as optional argument the scalar type, default being
	 * fast floating point .*/
	set_object(std::string s, number_type t = FLOAT);

	/** Create a continuous set from an implementation object.
	 *
	 * @note The implementation object is adopted, i.e., there
	 * should be no modification from the outside.
	 * It will be destroyed automatically when no longer
	 * used. */
	set_object(continuous_set_ptr p, number_type t);

	/** Create a continuous set from an implementation object.
	 *
	 * @note The implementation object is adopted, i.e., there
	 * should be no modification from the outside.
	 * It will be destroyed automatically when no longer
	 * used. */
	//set_object(continuous_set_ptr p);

	/** Shallow copy constructor */
	set_object(const set_object& X);

	/** Deep copy. */
	set_object clone() const;

	/** Shallow assignment. */
	set_object& operator=(const set_object& X);


	/** Returns whether the set is empty. */
	bool is_empty() const;

	/** Returns whether the set contains the entire space. */
	bool is_universe() const;

	/** Returns whether the intersection with X is empty. */
	bool is_disjoint_from(const set_object& X) const;

	/** Conservative containment check.
	 *
	 * If true, the set contains X. If false, X might not
	 * be contained, or the precision might not suffice
	 * to give a definite result.
	 */
	bool contains(const set_object& X) const;

	/** Remove redundant elements */
	void simplify();

	/** Get the domain. */
	positional_domain get_positional_domain() const;

	/** Get a shared pointer to the implementation object. */
	continuous_set_const_ptr get_impl() const;

	/** Get a shared pointer to the implementation object. */
	continuous_set_ptr get_impl();

	/** Get the number type. */
	number_type get_number_type() const;

	/** Set the number type. */
	void set_number_type(number_type t);

private:
	continuous_set_ptr my_cset;
	number_type my_number_type;
};

}
}

std::ostream& operator<<(std::ostream& os, const continuous::object::set_object& p);

#endif /* SET_OBJECT_H_ */
