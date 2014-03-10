/*
 * sf_set.h
 *
 *  Created on: Apr 1, 2010
 *      Author: frehse
 */

#ifndef SF_SET_H_
#define SF_SET_H_

#include "core/continuous/support_function_provider.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "math/vdom/affine_map.h"
#include "math/vdom/vdom_vector.h"
#include "math/vdom/lin_constraint.h"
#include "math/numeric/container_comp.h"

namespace continuous {
namespace support_function {

/** An interface class for implicit representations of continuous sets by means of
 * their support function.
 *
 * The sf_set class defines the interface, and sf_unary is an implementation
 * for a single continuous set that has the capacity to compute the support function,
 * i.e., solve a linear optimization problem over its elements.
 *
 * The sf_set class also provides a member my_map, which is an affine map that is
 * used to transform sets. Transforming the set (including existential quantification)
 * will simply change my_map, and will propagate any changes to the underlying set.
 * The functions map and existentially_quantify_variables are implemented by
 * changing my_map.
 *
 * There are also default implementations for outer_poly and is_empty, which
 * call the compute_support function.
 *
 * An auxiliary function compute_support_mapped is provided to help derived classes
 * implement their own compute_support functions, in which they must take my_map into
 * account.
 */

template<typename scalar_type> class sf_set: public support_function_provider
//	public boost::enable_shared_from_this<continuous_set>
{
public:
	typedef sf_set<scalar_type> self_type;
	typedef boost::shared_ptr<self_type> ptr;
	typedef boost::shared_ptr<const self_type> const_ptr;
	typedef math::affine_map<scalar_type> affine_map;
	typedef boost::shared_ptr<affine_map> affine_map_ptr;
	typedef boost::shared_ptr<const affine_map> affine_map_const_ptr;
	typedef math::vdom_vector<scalar_type> vector_type;
	typedef math::lin_constraint<scalar_type> lin_constraint;
	typedef constr_polyhedron<scalar_type> poly;
	typedef typename poly::ptr poly_ptr;
	//typedef std::set<vector_type,math::numeric::comp::lex_less_comp> vector_set;
	typedef std::set<vector_type, math::numeric::lex_comp_less<scalar_type,
			math::vdom_vector> > vector_set;

	//	/** Return a shared_ptr to *this. */
	//	virtual continuous_set::ptr get_ptr();
	//
	//	/** Return a shared_ptr to const *this. */
	//	virtual continuous_set::const_ptr get_const_ptr() const;

	/** Default constructor without affine map. */
	sf_set();

	/** Constructor for sets transformed by an affine map M. */
	sf_set(const affine_map& M);

	virtual ~sf_set();

	/** Returns the outer polyhedral approximation for a given set of directions.
	 */
	virtual poly_ptr outer_poly(const vector_set& directions) const;

	/** Apply an affine map> */
	virtual void map(const affine_map& M);

	/** Shallow copy. */
	virtual self_type* clone() const = 0;

	/** Returns the universe.
	 *
	 * Constructs an sf_unary object. */
	virtual self_type* create_universe() const;

	/** Returns the empty set.
	 *
	 * Constructs an sf_unary object. */
	virtual self_type* create_empty() const;

	/** Returns the memory occupied by *this. */
	virtual int get_memory() const = 0;

	/** Returns the dimension of *this. */
	virtual unsigned int get_dim() const;

	/** Adds the variables to the domain of *this. */
	virtual void embed_variables(const variable_id_set& id_set);

	/** Existentially quantifies over the domain of *this. */
	virtual void
	existentially_quantify_variables(const variable_id_set& id_set);

	/** Returns the set in predicate form. */
	virtual continuous_set_predicate::ptr get_predicate() const = 0;

//	/** Accept a visitor. */
//	virtual void
//	accept(dispatching::dispatcher<continuous_set_typelist>& d) const;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const = 0;

	//-------------------------------------------------------
	// primed_variable_provider functions
	//-------------------------------------------------------
	/** Returns the ids of all variables over which the set is defined. */
	virtual const variable_id_set& get_variable_ids() const = 0;

	/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
	 */
	virtual void reassign_primedness(unsigned int, unsigned int = 0) = 0;

	/** Increase the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, increase all. */
	virtual void increase_primedness(unsigned int = 0) = 0;

	/** Decrease the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, decrease all. */
	virtual void decrease_primedness(unsigned int = 0) = 0;

	//-------------------------------------------------------
	// support_function_provider functions
	//-------------------------------------------------------
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

	/** Returns the affine map associated with the object */
	affine_map_const_ptr get_map() const;

protected:
	/** Compute the support function with a given map.
	 *
	 * Note that max v.x, x \in AS+b = max v.Ay+v.b, y \in S.
	 * */
	template<typename fun_type> static void
	compute_support_mapped(const support_function_provider& a_set,
			const affine_map_const_ptr& a_map, const math::vdom_vector<
					fun_type>& l, fun_type& max_value,
			math::vdom_vector<fun_type>& support_vec, bool& is_empty,
			bool& is_bounded);

	/** Returns the affine map associated with the object */
	affine_map_ptr get_map();

	/** Set the affine map associated with the object */
	void set_map(const affine_map& M);

private:
	affine_map_ptr my_map;
};

}
}

#include "sf_set.hpp"

#endif /* SF_SET_H_ */
