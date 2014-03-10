/*
 * sf_lin_con_intersection.h
 *
 *  Created on: Apr 27, 2010
 *      Author: frehse
 */

#ifndef SF_LIN_CON_INTERSECTION_H_
#define SF_LIN_CON_INTERSECTION_H_

#include "core/continuous/support_function/sf_base/sf_unary.h"

namespace continuous {
namespace support_function {

/** An implicit representation of a continuous set by means of its support function.
 */

template<typename scalar_type> class sf_lin_con_intersection: public sf_unary<
		scalar_type> {
public:
	typedef boost::shared_ptr<sf_lin_con_intersection<scalar_type> > ptr;
	typedef boost::shared_ptr<const sf_lin_con_intersection<scalar_type> > const_ptr;

	typedef typename sf_unary<scalar_type>::affine_map affine_map;

	/** Construct a support function set representation of the intersection
	 * of a continuous set with the constraint con.
	 *
	 * The variables of con must be a subset of the variables of s. */
	sf_lin_con_intersection(const support_function_provider::const_ptr& s,
			const math::lin_constraint<scalar_type>& con, const std::string& minbrak_type,
			const double& intersection_error);

	/** Construct a support function set representation of the intersection
	 * of a continuous set with the constraint con, transformed by an affine map M.
	 *
	 * The variables of con must be a subset of the variables of s. */
	sf_lin_con_intersection(const support_function_provider::const_ptr& s,
			const math::lin_constraint<scalar_type>& con, const affine_map& M,
			const std::string& minbrak_type, const double& intersection_error);

	virtual ~sf_lin_con_intersection();

	/** Shallow copy. */
	virtual sf_lin_con_intersection<scalar_type>* clone() const;

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
		return false;
	}

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<Rational>& l,
			Rational& max_value, math::vdom_vector<Rational>& support_vec,
			bool& is_empty, bool& is_bounded) const ;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<double>& l,
			double& max_value, math::vdom_vector<double>& support_vec,
			bool& is_empty, bool& is_bounded) const;

	/**
	 * Computes the support function of the intersection set with a bounded number of
	 * function samplings in the underlying lower bound search algorithm.
	 *
	 * @param sample_bound The bound on the number of function sampling in the
	 * lower bound search algorithm.
	 */
	virtual void bounded_compute_support(const math::vdom_vector<double>& l,
				double& max_value, math::vdom_vector<double>& support_vec,
				const unsigned int sample_bound,
				bool& is_empty, bool& is_bounded) const;

	//-------------------------------------------------------
	// special functions
	//-------------------------------------------------------
	virtual const math::lin_constraint<scalar_type>& get_constraint() const;

protected:
	math::lin_constraint<scalar_type> my_con;
	scalar_type my_intersection_error;
	std::string my_minbrak_type;
};

}
}

#include "sf_lin_con_intersection.hpp"

#endif /* SF_LIN_CON_INTERSECTION_H_ */
