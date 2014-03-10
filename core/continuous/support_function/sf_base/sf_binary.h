/*
 * sf_binary.h
 *
 *  Created on: Apr 21, 2010
 *      Author: frehse
 */

#ifndef SF_BINARY_H_
#define SF_BINARY_H_

#include "sf_set.h"

namespace continuous {
namespace support_function {

/** A base class for sf_set implementations that use two support function
 * providers.
 */

template<typename scalar_type> class sf_binary: public sf_set<scalar_type> {
public:
	typedef typename sf_set<scalar_type>::affine_map affine_map;

	/** Construct a support function set representation of the result of a binary operation on
	 * two continuous sets. */
	sf_binary(const support_function_provider::const_ptr& s1,
			const support_function_provider::const_ptr& s2) :
		my_set1(s1), my_set2(s2) {
	}
	;

	/** Construct a support function set representation of the result of a binary operation on
	 * two continuous sets, transformed by an affine map M. */
	sf_binary(const support_function_provider::const_ptr& s1,
			const support_function_provider::const_ptr& s2, const affine_map& M) :
		sf_set<scalar_type> (M), my_set1(s1), my_set2(s2) {
	}
	;

	virtual ~sf_binary() {
	}
	;

	virtual sf_binary<scalar_type>* clone() const = 0;

	/** Returns the memory occupied by *this. */
	virtual int get_memory() const {
		throw std::runtime_error(
				"sf_binary : missing implementation get_memory");
	}

	/** Returns the set in predicate form. */
	virtual continuous_set_predicate::ptr get_predicate() const {
		throw std::runtime_error(
				"sf_binary : missing implementation decrease_primedness");
	}

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const {
		throw std::runtime_error(
				"sf_binary : missing implementation print");
	}

	//-------------------------------------------------------
	// primed_variable_provider functions
	//-------------------------------------------------------
	/** Returns the ids of all variables over which the set is defined. */
	virtual const variable_id_set& get_variable_ids() const {
		if (!this->get_map() || this->get_map()->is_empty()) {
			variable_id_set vis1 = my_set1->get_variable_ids();
			variable_id_set vis2 = my_set2->get_variable_ids();
			// @todo make this a member variable that is updated on the fly
			static variable_id_set vis;
			vis = variable_id_set();
			std::set_union(vis1.begin(), vis1.end(), vis2.begin(), vis2.end(),
					std::inserter(vis,vis.begin()));
			return vis;
		} else {
			return this->get_map()->codomain().get_variable_ids();
		}
	}
	;

	/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
	 */
	virtual void reassign_primedness(unsigned int, unsigned int = 0) {
		throw std::runtime_error(
				"sf_binary : missing implementation reassign_primedness");
	}

	/** Increase the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, increase all. */
	virtual void increase_primedness(unsigned int = 0) {
		throw std::runtime_error(
				"sf_binary : missing implementation increase_primedness");
	}
	/** Decrease the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, decrease all. */
	virtual void decrease_primedness(unsigned int = 0) {
		throw std::runtime_error(
				"sf_binary : missing implementation decrease_primedness");
	}
	//-------------------------------------------------------
	// support_function_provider functions
	//-------------------------------------------------------
	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<Rational>& l,
			Rational& max_value, math::vdom_vector<Rational>& support_vec,
			bool& is_empty, bool& is_bounded) const = 0;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<double>& l,
			double& max_value, math::vdom_vector<double>& support_vec,
			bool& is_empty, bool& is_bounded) const = 0;

protected:
	support_function_provider::const_ptr my_set1;
	support_function_provider::const_ptr my_set2;
};

}
}

#endif /* SF_BINARY_H_ */
