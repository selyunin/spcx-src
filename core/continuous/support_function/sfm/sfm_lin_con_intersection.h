/*
 * sfm_lin_cons_intersection.h
 *
 *  Created on: May 29, 2010
 *      Author: ray
 */

#ifndef SFM_LIN_CON_INTERSECTION_H_
#define SFM_LIN_CON_INTERSECTION_H_

#include "core/continuous/support_function/sf_base/sf_unary.h"
#include "core/continuous/support_function/sfm/sfm_guard_intersection.h"

namespace continuous {
namespace support_function {

/** An implicit representation of a continuous set by means of its support function.
 * The continuous set being represented is the intersection of an SFM flowpipe section
 * with a linear constraint (Guard). The section of the SFM flowpipe is determined with
 * the interval of the SFM indices passed to the constructor.
 *
 */

template<typename scalar_type> class sfm_lin_con_intersection: public sf_unary<
		scalar_type> {
public:
	typedef typename std::set<math::vdom_vector<scalar_type>,
							math::numeric::lex_comp_less<scalar_type,
									math::vdom_vector> > vector_set;
	typedef boost::shared_ptr<sfm_lin_con_intersection<scalar_type> > ptr;

	typedef typename sf_unary<scalar_type>::affine_map affine_map;

	/** Construct a support function set representation of the intersection
	 * of a continuous set with the constraint con.
	 *
	 * The variables of con must be a subset of the variables of s. */
	sfm_lin_con_intersection(const typename sfm_cont_set<scalar_type>::const_ptr& s,
			const math::numeric::interval<unsigned int>& my_intv,
			const math::lin_constraint<scalar_type>& con,
			const std::string& minbrak_type,
			const double& inters_error);

	/** Construct a support function set representation of the intersection
	 * of a continuous set with the constraint con, transformed by an affine map M.
	 *
	 * The variables of con must be a subset of the variables of s. */
	sfm_lin_con_intersection(const typename sfm_cont_set<scalar_type>::const_ptr& s,
			const math::numeric::interval<unsigned int>& my_intv,
			const math::lin_constraint<scalar_type>& con, const affine_map& M,
			const std::string& minbrak_type,
			const double& inters_error);

	virtual ~sfm_lin_con_intersection();

	/** Shallow copy. */
	virtual sfm_lin_con_intersection<scalar_type>* clone() const;

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

	/**
	 * Computes the outer polytopes of the sfm section member
	 * sets and the linear constraint intersection, on the passed directions.
	 * The result returned is a polyhedron collection object.
	 * This call the simultaneous lb_search routine to compute the result optimally.
	 *
	 * @param directions
	 * @param split_size
	 * @return
	 */
	std::vector<typename polyhedron<scalar_type>::ptr> get_outer_polys(vector_set directions) const;

	/**
	 * Computes the outer polytopes of the sfm section member
	 * sets and the linear constraint intersection, on the passed directions.
	 * The result returned is a polyhedron collection object.
	 * The split size parameter denotes the sfm interval size considered to compute
	 * the intersection with the linear constraint. By default, split size is set
	 * to 1 meaning that all the individual sfm section members will be intersected
	 * with the guard constraint and the number of returned outer polys will be
	 * same as the sfm section size.
	 */

	std::vector<typename polyhedron<scalar_type>::ptr> get_chull_outer_polys(vector_set directions, size_t split_size=1) const;
	//-------------------------------------------------------
	// special functions
	//-------------------------------------------------------
	virtual const math::lin_constraint<scalar_type>& get_constraint() const;

protected:
	math::lin_constraint<scalar_type> my_con;
	math::numeric::interval<unsigned int> my_interval;
	std::string my_minbrak_type;
	double my_intersection_error;
};

} // end of namespace support_function
} // end of namespace continuous

#include "sfm_lin_con_intersection.hpp"

#endif /* SFM_LIN_CON_INTERSECTION_H_ */
