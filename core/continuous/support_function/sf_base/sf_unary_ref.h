/*
 * sf_unary_ref_ref.h
 *
 *  Created on: Mar 31, 2011
 *      Author: frehse
 */

#if !defined(SF_UNARY_REF_H_) || !defined(SF_UNARY_REF_H__BODY)
//#ifndef SF_UNARY_REF_H_
#define SF_UNARY_REF_H_

#include "sf_set.h"

#ifndef SF_UNARY_REF_H__BODY
#define SF_UNARY_REF_H__BODY

namespace continuous {
namespace support_function {

/** An implicit representation of a continuous set by means of its support function,
 *  given as a const reference.
 *
 *  This differs from sf_unary in that the support_function_provider is passed
 *  as a const reference, and therefore does not need to be owned by a
 *  shared_ptr.
 */

template<typename scalar_type, class implementor = support_function_provider> class sf_unary_ref: public sf_set<scalar_type>
{
public:
	typedef typename sf_set<scalar_type>::affine_map affine_map;

	/** Construct a support function set representation from a continuous set
	 * that provides a support function. */
	explicit sf_unary_ref(const implementor& s);

	/** Construct a support function set representation from a continuous set
	 * that provides a support function, transformed by an affine map M. */
	sf_unary_ref(const implementor& s, const affine_map& M);

	virtual ~sf_unary_ref();

	/** Shallow copy.
	 *
	 * The clone shares the support_function_provider with *this. */
	virtual sf_unary_ref<scalar_type,implementor>* clone() const;

	/** Returns the memory occupied by *this. */
	virtual int get_memory() const;

	/** Returns the set in predicate form. */
	virtual continuous_set_predicate::ptr get_predicate() const;

	/** Returns whether *this is universe. */
	virtual math::tribool is_universe() const;

	/** Returns whether *this is empty. */
	virtual math::tribool is_empty() const;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const;

	//-------------------------------------------------------
	// primed_variable_provider functions
	//-------------------------------------------------------
	/** Returns the ids of all variables over which the set is defined. */
	virtual const variable_id_set& get_variable_ids() const;

	/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
	 */
	virtual void reassign_primedness(unsigned int, unsigned int = 0);

	/** Increase the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, increase all. */
	virtual void increase_primedness(unsigned int = 0);

	/** Decrease the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, decrease all. */
	virtual void decrease_primedness(unsigned int = 0);

	//-------------------------------------------------------
	// support_function_provider functions
	//-------------------------------------------------------
	/** Returns true if compute_support returns a support vector
		 * and false otherwise.
		 */
	virtual bool computes_support_vector() const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(
			const math::vdom_vector<Rational>& l,
			Rational& max_value,
			math::vdom_vector<Rational>& support_vec, bool& is_empty,
			bool& is_bounded) const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<double>& l,
			double& max_value, math::vdom_vector<double>& support_vec,
			bool& is_empty, bool& is_bounded) const;

protected:
	/** Compute the support function with a given map.
	 *
	 * Note that max v.x, x \in AS+b = max v.Ay+v.b, y \in S.
	 * */
	template<typename fun_type> static void
	compute_support_mapped(const implementor& a_set,
			const typename sf_set<scalar_type>::affine_map_const_ptr& a_map, const math::vdom_vector<
					fun_type>& l, fun_type& max_value,
			math::vdom_vector<fun_type>& support_vec, bool& is_empty,
			bool& is_bounded);

	const implementor& my_set;
};

}
}

#include "sf_unary_ref.hpp"

#endif /* SF_UNARY_REF_H__BODY */
#endif /* SF_UNARY_REF_H_ */
