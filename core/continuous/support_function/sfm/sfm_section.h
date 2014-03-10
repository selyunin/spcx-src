/*
 * sfm_section.h
 *
 *  Created on: May 3, 2010
 *      Author: frehse
 */

#ifndef SFM_SECTION_H_
#define SFM_SECTION_H_

#include "core/continuous/support_function/sf_base/sf_set.h"
#include "sfm_cont_set.h"

namespace continuous {
namespace support_function {

/** A representation of a connected subset of a sfm_cont_set.
 */

template<typename scalar_type> class sfm_section: public sf_set<scalar_type> {
public:
	typedef typename boost::shared_ptr<sfm_section<scalar_type> > ptr;
	typedef typename boost::shared_ptr<const sfm_section<scalar_type> > const_ptr;
	typedef typename sfm_cont_set<scalar_type>::ptr sfm_ptr;
	typedef typename sf_set<scalar_type>::affine_map_const_ptr affine_map_const_ptr;
	typedef typename sf_set<scalar_type>::affine_map affine_map;
	typedef typename sfm_cont_set<scalar_type>::index_interval index_interval;


	/** Construct a support function set representation from a continuous set
	 * that provides a support function.`
	 *
	 * If compact_this is true, the sfm_section is compacted, unless
	 * forbidden by the static member compact_sfms. */
	sfm_section(const sfm_ptr& s, index_interval intv, bool compact_this=true);

	/** Construct a support function set representation from a continuous set
	 * that provides a support function, transformed by an affine map M. */
	sfm_section(const sfm_ptr& s, index_interval intv, const affine_map& M);

	virtual ~sfm_section();

	/** Shallow copy. */
	virtual sfm_section<scalar_type>* clone() const;

	/** Returns the memory occupied by *this. */
	virtual int get_memory() const;

	/** Returns the set in predicate form. */
	virtual continuous_set_predicate::ptr get_predicate() const;

	/** Returns whether *this is empty. */
	virtual math::tribool is_empty() const;

	/** Returns whether *this is universe. */
	virtual math::tribool is_universe() const;

	/*
	 * returns the sfm index interval.
	 */
	virtual index_interval get_interval() const;
	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const;

	/** Reduce the sfm to the interval of *this */
	virtual void compact();

	/** Obtain the underlying sfm set */
	virtual sfm_ptr get_sfm_set();

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
	/** Compute the support function with a given map. */
	template<typename fun_type> static void
	compute_support_mapped(const sfm_ptr& a_set, index_interval intv,
			const affine_map_const_ptr& a_map, const math::vdom_vector<
					fun_type>& l, fun_type& max_value, math::vdom_vector<
					fun_type>& support_vec, bool& is_empty, bool& is_bounded);

	/** Compute the support function on a given interval. */
	template<typename fun_type> static void
	compute_support_interval(const sfm_ptr& a_set, index_interval intv,
			const math::vdom_vector<fun_type>& l, fun_type& max_value,
			math::vdom_vector<fun_type>& support_vec, bool& is_empty,
			bool& is_bounded);

	sfm_ptr my_set;
	index_interval my_intv;

public:
	static bool compact_sfms;

};
}
}

#include "sfm_section.hpp"

#endif /* SFM_SECTION_H_ */
