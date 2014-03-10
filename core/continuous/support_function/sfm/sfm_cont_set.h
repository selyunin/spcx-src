#ifndef SFM_CONT_SET
#define SFM_CONT_SET

/***************************************************************************
 *   Copyright (C) 2008 by Goran Frehse   *
 *   goran.frehse@imag.fr   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "core/continuous/support_function_provider.h"
#include "core/continuous/continuous_set_collection.h"
#include "core/continuous/polyhedra/polyhedron_collection.h"
#include "math/vdom/index_to_variable_id_map_provider.h"
#include "math/matrix.h"
#include "math/scalar_types/scalar_with_infinity.h"
#include "math/numeric/interval.h"
#include "math/vdom/lin_constraint.h"
#include "math/unique_vector_to_value_store.h"

/** Forward declarations */
namespace continuous {
namespace support_function {

template<class T>
class sf_evaluator;

}
}

namespace continuous {
namespace support_function {

/** Forward declaration of friend class
 */
template<typename scalar_type> class sfm_section;

/**
 * This class encapsulates the representation of continuous sets as Support Function Matrix.
 * A support function matrix is defined for a choosen set of directions, say D. A direction is
 * a vector \f$l \in R^n\f$.
 * A support function matrix is a k X N matrix where k is the number of directions in D and N is the
 * time horizon. The \f$ (i,j)^th \f$ entry denotes the support function of the reachable set at the
 * jth iteration of post_c computation using support functions in the direction \f$ l_i \f$.
 *
 * The root of an SFM is a set such that post_c(root) \subseteq SFM.
 *
 * @author Rajarshi Ray
 */

/**
 * This structure to hold all the information that uniquely defines a post_c problem
 * @author Rajarshi
 */
template <typename scalar_type>
struct postc_params{
	double delta; // time step
	double time_horizon; // The time upper bound on the post_c iterations
	math::matrix<scalar_type> dynamics_A; // The matrix e^{delta*A} saying the dynamics.
	math::vector<scalar_type> dynamics_b;
	boost::shared_ptr<const continuous_set> invariant_set_ptr;
	boost::shared_ptr<const continuous_set> input_set_ptr;
	boost::shared_ptr<const continuous_set> initial_set_ptr;
	std::vector<scalar_type> delta_vec;
};

template <typename scalar_type> class sfm_cont_set : public index_to_variable_id_map_provider, public support_function_provider
//	public boost::enable_shared_from_this<continuous_set>
{
public:
	typedef boost::shared_ptr<sfm_cont_set<scalar_type> > ptr;
	typedef boost::shared_ptr<const sfm_cont_set<scalar_type> > const_ptr;
	typedef math::matrix<scalar_type> matrix_type;
	typedef math::vector<scalar_type> vector_type;
	typedef std::list<vector_type> vector_list;
	typedef vector_type direction;
	typedef math::numeric::interval<unsigned int> index_interval;
	typedef typename math::matrix<scalar_type>::size_type index_type; // take the same type as is used by the matrix
	typedef math::unique_vector_to_value_store<scalar_type, math::vector,
			index_type> direction_store; // stores the row index for each vector

	typedef enum {TEXTUAL, DOUBLE_CONSTRAINTS, DOUBLE_GENERATORS, JVX} output_format;

//	/** Return a shared_ptr to *this. */
//	virtual continuous_set::ptr get_ptr();
//
//	/** Return a shared_ptr to const *this. */
//	virtual continuous_set::const_ptr get_const_ptr() const;

	/** Initialize with an element. */
	sfm_cont_set(
					const postc_params<scalar_type>& postc_prb,
					const matrix_type& sfm,
					const direction_store& dirs,
					const index_to_variable_id_map_ptr& iimap,
					const math::lin_constraint_system<scalar_type>& nonstate_contraints);

	/** Initialize with an element. */
	sfm_cont_set(const postc_params<scalar_type>& postc_prb, const matrix_type& sfm, const direction_store& dirs, const index_to_variable_id_map_ptr&);

	/** Initialize without postc_prb. */
	sfm_cont_set(const matrix_type& sfm, const direction_store& dirs, const index_to_variable_id_map_ptr&);

	virtual ~sfm_cont_set();

	virtual const matrix_type& get_sfm() const;
	virtual void set_sfm(matrix_type new_sfm);

	virtual const direction_store& get_directions() const;
	virtual void set_directions(const direction_store& new_dirs);

	/** Returns the number of objects (columns of the SFM) defined. */
	virtual size_t get_size() const;

	virtual const postc_params<scalar_type>& get_postc_prb() const;

//	/**
//	 * Computes an SFM representation of the transformation AX + v.
//	 *
//	 * X is a continuous set represented by *this. If A is non-singular, then an exact image of the map is computed.
//	 * If A is singular then we compute a polyhedral tight over-approximation of the exact image.
//	 *
//	 * @param M an affine map
//	 * @return a SFM which represents the overapproximation of the transformed set.
//	 */
//	virtual sfm_cont_set<scalar_type> affine_transform(const math::affine_map<scalar_type>& M) const;

	/** Computes \f$ \Omega_j \f$ as facets represented polytope for a given parameter j.
	 *
	 * If j is more than the time horizon N, an exception is thrown
	 * (run_time_error("Invalid Parameter j - should be less than N")).
	 *
	 * @param j Iteration index of the \f$post_c \f$ computation using support functions.
	 * @return The reach set at jth iteration of \f$ post_c \f with support functions$, as polytope.
	 */
	virtual constr_polyhedron<scalar_type> get_polytope(unsigned int j) const;

	/** Returns the outer polyhedral approximation of the flowpipe from step i to j including. */
	virtual polyhedron_collection<scalar_type> get_outer_polytope_collection(unsigned int i, unsigned int j) const;

	/** Returns the outer polyhedral approximation of the flowpipe. */
	virtual polyhedron_collection<scalar_type> get_outer_polytope_collection() const;

	/** Return true if Omega_j is empty. */
	virtual bool is_column_empty(unsigned int j) const;

	/** Remove all Omega_j that are empty. */
	virtual void remove_empty_columns();

//	/**
//	 * @param The normal vector to the guard half space d
//	 * @param The distance of the guard from the origin b
//	 * @return The first column index j such that \omega_j intersects with the guard set which is a half space
//	 */
//	virtual unsigned int get_first_j(direction d, scalar_type b) const;

	/**
	 * Adds a new row to the SFM with the support function evaluations of the new direction.
	 * If the new direction is already present in the list of directions, then simply return the
	 * row index corresponding to the new direction
	 *
	 *@param New direction d
	 *@return Row index i and scaling factor f for which the direction D_i = f*d
	 */
	virtual std::pair<unsigned int,scalar_type> extend_sfm(const direction& d);

	/** Returns the polyhedron that is the outer polyhedral approximation in
	 * the directions of *this of the Omega_i with i in intv.
	 */
	virtual constr_polyhedron<scalar_type> compute_template_hull(index_interval intv=index_interval()) const;

	/** Extend the sfm to this direction and replace with con wherever con is tighter. */
	void intersection_with_constraint(const math::lin_constraint<scalar_type>& con, bool add_to_ini = true);

	/**
	 * Computes the SFM thats represents the intersection of *this with
	 * a constraint represented polytope,i.e., polytope as facets.
	 *
	 * @param A polytope represented with facets.
	 * @return continuous set pointer to a SFM_cont_set.
	 *
	 * @todo Map *this and con_poly to the same variables. Current version can
	 * throw if there is a variable in con_poly that is not in *this.
	 * (reordering of coefficients can throw).
	 */
	void intersection_with_poly(const polyhedron<scalar_type>& poly);

	/**
	 * With some optimizations
	 * @param poly
	 * @return intersection set
	 */
	virtual ptr intersection_with_poly_improved(const polyhedron<scalar_type>& poly) const;

	/** Returns true if the initial set of *this contains the initial set of s.
	 *
	 * If either initial set is null, return false. */
	virtual bool contains_initial_set(const sfm_cont_set& s) const;

	/** Returns true if the initial set of *this contains a set.
	 *
	 * If initial set is null, return false. */
	virtual bool contains_initially(const continuous_set::const_ptr& s) const;

	/** Copy
	 *
	 * Uses normal copy constructor. This makes a copy of the sfm, but only a shallow
	 * copy of the problem instance.
	 */
	virtual sfm_cont_set<scalar_type>* clone() const;

	/**
	 * Create empty sfm with no directions, empty matrix and empty index_to_id_map
	 */
	virtual sfm_cont_set<scalar_type>* create_universe() const;
	virtual sfm_cont_set<scalar_type>* create_empty() const;

	/** Create an empty sfm. */
	static sfm_cont_set<scalar_type>* empty_set();

	virtual int get_memory() const;

	virtual unsigned int get_dim() const;

	/** The collection is by construction empty iff all the elements is empty
	 * or the root is empty. */
	virtual math::tribool is_empty() const;
	virtual math::tribool is_universe() const;

	/** Return whether the outer polytope approximation is empty */
	virtual math::tribool is_empty_outer() const;

	virtual void embed_variables(const variable_id_set& id_set);
	virtual void existentially_quantify_variables(const variable_id_set& id_set);

	/** Returns the set in predicate form. */
	virtual continuous_set_predicate::ptr get_predicate() const {
		throw std::runtime_error("missing implementation get_predicate");
		return continuous_set_predicate::ptr();
	};

	/** Accept a visitor. */
	virtual void accept(dispatching::dispatcher<continuous_set_typelist>& d) const;
	/** Set the output format of the set.
	 */
	static const output_format& get_output_format() {
		return my_output_format;
	}
	;
	static void set_output_format(output_format new_format) {
		my_output_format=new_format;
	}
	;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const;

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

private:
	matrix_type _my_sfm;
	//vector_list directions;
	direction_store my_directions;
	postc_params<scalar_type> my_postc_prb;
protected:
	math::lin_constraint_system<scalar_type> my_nonstate_contraints;
private:
	boost::shared_ptr<sf_evaluator<scalar_type> > my_evaluator;

	static output_format my_output_format;
	//positional_vdomain my_pos_dom;

	friend class sfm_section<scalar_type>;
};

/** Utility function to convert the direction list into a direction_store
 *
 * The indices are assigned in the order of the list.
 * */
template <typename scalar_type> typename sfm_cont_set<scalar_type>::direction_store convert_to_direction_store(const typename sfm_cont_set<scalar_type>::vector_list& dirs ) {
	typename sfm_cont_set<scalar_type>::direction_store store;
	typename sfm_cont_set<scalar_type>::index_type i = 0;
	for (typename sfm_cont_set<scalar_type>::vector_list::const_iterator it=dirs.begin();it!=dirs.end();++it) {
		store.insert_missing(*it,i);
		++i;
	}
	return store;
}

}
}

#include "sfm_cont_set.hpp"

#endif /*SFM_CONT_SET*/

