/*
 * sf_ode_solv.h
 *
 *  Created on: Oct 11, 2010
 *      Author: frehse
 */

#ifndef SF_ODE_SOLV_H_
#define SF_ODE_SOLV_H_

#include "core/continuous/support_function/sf_base/sf_unary.h"

namespace continuous {
namespace support_function {

/** An implicit representation of a continuous set by means of its support function.
 */

template<typename scalar_type> class sf_ode_solv: public sf_unary<scalar_type> {
public:
	typedef typename sf_unary<scalar_type>::affine_map affine_map;

	/** Construct a support function set representation of
	 * the flowpipe from X0 according to dynamics dyn
	 * up to time delta.
	 * */
	sf_ode_solv(const support_function_provider::const_ptr& X0,
			const typed_dynamics<scalar_type>& dyn, scalar_type delta);

	/** Construct a support function set representation of
	 * the flowpipe from X0 according to dynamics dyn
	 * up to time delta.
	 * */
	sf_ode_solv(const support_function_provider::const_ptr& X0,
			const typed_dynamics<scalar_type>& dyn, scalar_type delta, const affine_map& M);

	virtual ~sf_ode_solv();

	/** Shallow copy. */
	virtual sf_ode_solv<scalar_type>* clone() const;

	/** Returns the memory occupied by *this. */
	virtual int get_memory() const;

	/** Returns the set in predicate form. */
	virtual continuous_set_predicate::ptr get_predicate() const;

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
	virtual bool computes_support_vector() const {
		return true;
	}

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<Rational>& l,
			Rational& max_value, math::vdom_vector<Rational>& support_vec,
			bool& is_empty, bool& is_bounded) const;

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
	const typed_dynamics<scalar_type>& my_dyn;
	scalar_type my_delta;
};

}
}

#include "sf_ode_solv.hpp"

#endif /* SF_ODE_SOLV_H_ */
