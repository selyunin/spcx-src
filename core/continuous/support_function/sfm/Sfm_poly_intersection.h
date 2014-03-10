/*
 * Sfm_poly_intersection.h
 *
 *  Created on: Dec 30, 2010
 *      Author: ray
 */

#ifndef SFM_POLY_INTERSECTION_H_
#define SFM_POLY_INTERSECTION_H_

#include "sf_unary.h"

namespace continuous {
	namespace support_function {

template <typename scalar_type> class Sfm_poly_intersection : public sf_unary<scalar_type>{
public:
	typedef typename sf_unary<scalar_type>::affine_map affine_map;

	Sfm_poly_intersection(const typename sfm_cont_set<scalar_type>::const_ptr& s,
			const math::numeric::interval<unsigned int>& my_intv,
			const continuous::polyhedron<scalar_type>& cons_poly);
	/** Construct a support function set representation of the intersection
	 * of a continuous set with the constraint con, transformed by an affine map M.
	 *
	 * The variables of con must be a subset of the variables of s. */
	Sfm_poly_intersection(const typename sfm_cont_set<scalar_type>::const_ptr& s,
			const math::numeric::interval<unsigned int>& my_intv,
			const continuous::polyhedron<scalar_type>& cons_poly, const affine_map& M);


	/** Shallow copy. */
	virtual Sfm_poly_intersection<scalar_type>* clone() const;

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

	//-------------------------------------------------------
	// special functions
	//-------------------------------------------------------
	virtual const continuous::polyhedron<scalar_type>& get_poly() const;

	virtual ~Sfm_poly_intersection();

protected:
	continuous::polyhedron<scalar_type> my_poly;
	math::numeric::interval<unsigned int> my_interval;

	virtual ~Sfm_poly_intersection();
};

} // end of support_function namespace
} // end of continuous namespace

#endif /* SFM_POLY_INTERSECTION_H_ */

