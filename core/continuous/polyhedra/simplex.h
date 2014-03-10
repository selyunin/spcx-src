/*
 * simplex.h
 *
 *  Created on: Jul 8, 2011
 *      Author: frehse
 */

#ifndef SPACEEX_SIMPLEX_H_
#define SPACEEX_SIMPLEX_H_

#include "math/vdom/index_to_variable_id_map_provider.h"
//#include "math/vdom/affine_map.h"
#include "math/vector.h"
#include "math/vector_operators.h"
#include "math/unique_vector_to_value_store.h"
#include "polyhedron.h"
#include "core/continuous/continuous_dynamics/continuous_dynamics.h"

namespace continuous {

/** A simplex represented by its center and a generator.
 *
 * A simplex (c,g,dom) defines the set x given component-wise as
 *    x_i = c_i + u_i*g_i, -1 <= u_i <= 1
 * for all x_i in dom.
 *
 * A simplex can not be universe, nor zero-dimensional.
 * An exception is throws if such objects are requested
 * through calls of the base class interface.
 *
 * A simplex is empty iff there is a g_i<0. */
template<typename scalar_type>
class simplex: public polyhedron<scalar_type> ,
		public index_to_variable_id_map_provider,
		public boost::enable_shared_from_this<simplex<scalar_type> > {
public:
	typedef simplex<scalar_type> my_type;
	typedef scalar_type value_type;
	typedef math::vector<scalar_type> point_type;
	typedef math::vdom_vector<scalar_type> vdom_vector_type;
	typedef typename point_type::size_type size_type;

	typedef math::unique_vector_to_value_store<scalar_type, math::vector,
			bool> point_set_type;

	/** Inconsistent constructor */
	simplex();

	/** A simplex defined over dom.
	 */
	simplex(positional_vdomain dom);

	/** A regular simplex defined over dom and constained inside
	 * a ball of radius given
	*/
	simplex(positional_vdomain dom, scalar_type radius);

//	/** Copy constructor */
//	simplex(const simplex<scalar_type>& box);

//	/** Assignment constructor */
//	simplex<scalar_type>& operator=(
//			const simplex<scalar_type>& box);

	// ----------------------------------------
	// Continuous_set interface
	// ----------------------------------------

	/** Return a shared_ptr to *this. */
	continuous_set::ptr get_ptr();

	/** Return a shared_ptr to const *this. */
	continuous_set::const_ptr get_const_ptr() const;

	my_type* clone() const;

	my_type* create_universe() const;
	my_type* create_empty() const;

	int get_memory() const;

	unsigned int get_dim() const;

	math::tribool is_empty() const;
	math::tribool is_universe() const;

	void embed_variables(const variable_id_set& id_set);
	void
	existentially_quantify_variables(const variable_id_set& id_set);

	// Visitors are handled by polyhedron base class
	//	/** Accept a visitor. */
	//	 void
	//			accept(dispatching::dispatcher<continuous_set_typelist>& d) const;

	// ----------------------------------------
	// Support Function functions
	// ----------------------------------------

	/** Returns true if compute_support returns a support vector
	 * and false otherwise.
	 */
	bool computes_support_vector() const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	void compute_support(const math::vdom_vector<Rational>& l,
			Rational& max_value, math::vdom_vector<Rational>& support_vec,
			bool& is_empty, bool& is_bounded) const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	void compute_support(const math::vdom_vector<double>& l, double& max_value,
			math::vdom_vector<double>& support_vec, bool& is_empty,
			bool& is_bounded) const;

	// ----------------------------------------
	// Polyhedron functions
	// ----------------------------------------

	/** Get the constraints of *this. */
	typename math::lin_constraint_system<scalar_type>::const_ptr
	get_constraints() const;

	/** Get the vertices of *this. */
	const point_set_type& get_vertices() const;

	/*! Adds the constraint \p c to \p *this.
	 *
	 * Throws if the constraint is not a box constraint.
	 */
	void add_constraint(const math::lin_constraint<scalar_type> &c,
			bool check_redundancy = false);

	/*! Adds the vertice \p v to \p *this.
	 */
	void add_vertice(const point_type &c);

	/** Remove all redundant constraints.
	 *
	 * Does nothing, since no redundant constraints.
	 * */
	virtual void remove_redundant_constraints();

private:

	point_set_type my_vertices;
};

}

#include "simplex.hpp"

#endif /* SIMPLEX_H_ */
