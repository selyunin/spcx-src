#ifndef SUPPORT_FUNCTION_PROVIDER_H_
#define SUPPORT_FUNCTION_PROVIDER_H_

#include "core/continuous/continuous_set.h"
#include "math/vdom/lin_expression.h"

namespace continuous {

/* A class defining the interface for computation of support functions
 * on a set. */
class support_function_provider : public continuous_set {
public:
	typedef boost::shared_ptr<support_function_provider> ptr;
	typedef boost::shared_ptr<const support_function_provider> const_ptr;

	virtual ~support_function_provider() {
	}
	;

	/** Creates an identical copy of *this. */
	virtual support_function_provider* clone() const = 0;

	/** Accept a visitor. */
	virtual void accept(dispatching::dispatcher<continuous_set_typelist>& d) const {
		d.dispatch(this);
	}
	;

	/** Returns true if compute_support returns a support vector
	 * and false otherwise.
	 *
	 * The statement must hold for the current state of the object,
	 * and remain until the object is modified.
	 */
	virtual bool computes_support_vector() const = 0;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<Rational>& l,
			Rational& max_value,
			math::vdom_vector<Rational>& support_vec, bool& is_empty,
			bool& is_bounded) const = 0;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<double>& l,
			double& max_value, math::vdom_vector<double>& support_vec,
			bool& is_empty, bool& is_bounded) const = 0;

};


/** A converter for compute_support that provides the type closest to the passed type
 *
 * By default use double.
 */
template <typename scalar_type>
struct support_type {
	// by default use double
	typedef double type;
};

/** Use Rational for Rational. */
template <>
struct support_type<Rational> {
	typedef Rational type;
};

}

#endif /*SUPPORT_FUNCTION_PROVIDER_H_*/
