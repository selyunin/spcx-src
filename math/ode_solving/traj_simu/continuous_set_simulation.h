#ifndef continuous_set_SIMULATION_h
#define continuous_set_SIMULATION_h

#include "core/continuous/continuous_set.h"
//included in continuous_set.h
#include <iostream>
#include "boost/shared_ptr.hpp"
#include <boost/enable_shared_from_this.hpp>
#include "utility/shared_ptr_output.h"
#include "utility/printable.h"
#include "utility/tree_node.h"
#include "math/unique_vector_to_value_store.h"
#include "math/vdom/primed_variable_provider.h"
#include "math/vdom/index_to_variable_id_map_provider.h"
#include "core/continuous/continuous_set_declarations.h"
#include "utility/dispatching/double_dispatch.h"

//our includes
#include "math/vdom/positional_vdomain.h"
#include "core/continuous/polyhedra/hyperbox/finite_hyperbox.h"
#include "math/vdom/trajectory.h"

namespace continuous {
/**
 * A class defining continuous sets for simulation purpose.
 * The continuous set is defined by points and a zero-centered
 * hyperbox, which defines if a point is close enough of another one
 * in order to decide not to develop it
 *
 * @note The set is defined over a given variable domain.
 */

template<typename scalar_type>
class continuous_set_simulation: public support_function_provider, public index_to_variable_id_map_provider
/*public boost::enable_shared_from_this<continuous_set_simulation<scalar_type> >*/ {

public:

	typedef boost::shared_ptr<continuous_set_simulation> ptr;
	typedef boost::shared_ptr<const continuous_set_simulation> const_ptr;

	typedef math::trajectory<scalar_type> trajectory;
	typedef math::vdom_vector<scalar_type> state;
	typedef math::unique_vector_to_value_store<scalar_type,math::vdom_vector,trajectory> state_to_trajectory_map;
	typedef typename state_to_trajectory_map::iterator iterator;
	typedef typename state_to_trajectory_map::const_iterator const_iterator;
	typedef typename state_to_trajectory_map::value_type value_type;

	typedef math::vdom_vector<scalar_type> vdom_vector;
	typedef finite_hyperbox<scalar_type> fhyperbox;

	// --------------------------------------------
	// Constructors
	// The given hyperbox should be zero centered !!!
	/** Create a continuous_set_simulation object.
	 *
	 * The given hyperbox modelling the uncertainty
	 * around each trajectory.
	 */
	continuous_set_simulation(const fhyperbox  & hbox) : index_to_variable_id_map_provider(hbox.domain()),
		myhbox(hbox) {
	}


	// --------------------------------------------

	/** Creates an identical copy of *this. */
	virtual continuous_set_simulation* clone() const;


	/** Return a shared_ptr to *this. */
	//virtual continuous_set::ptr get_ptr();

	/** Return a shared_ptr to const *this. */
	//virtual continuous_set::const_ptr get_const_ptr() const;

	/** Creates a set of dimension zero, containing the entire state space
	 * (equivalent to true). */
	virtual continuous_set_simulation* create_universe() const;

	/** Creates an empty set of dimension zero (equivalent to false). */
	virtual continuous_set_simulation* create_empty() const;

	/** Virtual Destructor. */
	virtual ~continuous_set_simulation() {
	}
	;

	/* \} */
	// --------------------------------------------
	/** \name Non-modifying non-semantic methods
	 *  \{ */
	// --------------------------------------------

	/** Return a shared_ptr to *this. */
	//ptr get_ptr();

	/** Return a shared_ptr to const *this. */
	//const_ptr get_const_ptr() const;

	/** Returns the memory consumed by \p *this in bytes. */
	virtual int get_memory() const;

	/* \} */
	// --------------------------------------------
	/** \name Non-modifying semantic methods
	 *  \{ */
	// --------------------------------------------

	/** Returns the dimension of the state space.
	 *  A utility function.
	 */
	virtual dimension_t get_dim() const;

	/** Returns <CODE>true</CODE> if and only if \p *this is empty. */
	virtual math::tribool is_empty() const;

	/** Returns <CODE>true</CODE> if and only if \p *this contains
	 the entire state space.
	 The default implementation constructs the universe, subtracts \p *this
	 and tests for emptiness.
	 */
	virtual math::tribool is_universe() const;

	/** Returns <CODE>true</CODE> if and only if the intersection of \p *this and \p *ps is
	 empty.
	 The default implementation corresponds to the mathematical definition, but it may be
	 overridden by a more efficient implementation.
	 */
	//virtual bool is_disjoint_from(const continuous_set_const_ptr& ps) const;

	/** Returns <CODE>true</CODE> if and only if \p *this contains \p *ps.
	 The default implementation corresponds to difference assign followed by
	 checking emptiness, but it may be overridden by a more efficient implementation.
	 */
	//virtual bool contains(const continuous_set_const_ptr& ps) const;

	/** Iterface to a (possibly time consuming) simplication operator.
	 The default implementation is to do nothing.
	 */
	/*
	virtual void simplify() {
	}
	;*/

	/* \} */
	// --------------------------------------------
	/** \name Domain manipulation
	 *  \{ */
	// --------------------------------------------

	/** Expand the state space to incorporate the variables in id_set in the sense of embedding.
	 If an id in id_set is already in the state space, perform the operation using id_set without id.
	 */
	virtual void embed_variables(const variable_id_set& id_set);

	/** \brief Existential quantification over the variables in id_set.
	 *
	 If an id in id_set is not in the state space, perform the operation using id_set without id.
	 */
	virtual void existentially_quantify_variables(const variable_id_set& id_set);

	/* \} */
	// --------------------------------------------
	/** \name Other base class methods
	 *  \{ */
	// --------------------------------------------

	/** Output set to stream */
	virtual void print(std::ostream&) const{
		throw std::runtime_error("feature not implemented yet");
	}

	/* \} */

	/** Returns the set in predicate form.*/
	virtual continuous_set_predicate::ptr get_predicate() const;

	/**  Accept a const_visitor. A const_visitor must provide the function
	 * <code> void dispatch(const T* c) </code>
	 * for all derived classes of continuous_set that are listed in
	 * continuous_set_typelist (see continuous_set_declarations). */
	virtual void accept(const_visitor& d) const;

	/* \} */
	// --------------------------------------------
	/** \name Support function methods
	 *  \{ */
	// --------------------------------------------

	/** Returns true if compute_support returns a support vector
	 * and false otherwise.
	 *
	 * The statement must hold for the current state of the object,
	 * and remain until the object is modified.
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

	/* \} */
	// --------------------------------------------
	/** \name Trajectory set methods
	 *  \{ */
	// --------------------------------------------

	/**
	 *	Containment check implementation
	 *
	 * @param d a continuous set
	 * @return if the roots of *this contains the roots of d
	 * @todo implement an efficient containment check
	 * this one is only for dev/debug
	 **/
	bool contains_roots_of(const continuous_set_simulation& d) const;

	/** Check if any of the roots contain the state v */
	bool root_contains(const state& v) const;

	/** Intersect with another continuous_set_simulation object
	 *
	 * @todo To be implemented
	 * */
	continuous_set_ptr intersection_with(const continuous_set_simulation<scalar_type> & d) const;

	/** Intersect with a set of linear constraints
	 *
	 * Portions of trajectories are returned whose sampling points satisfy cons.
	 * No interpolation points are added. Non-contiguous satisfying intervals are stored in
	 * separate trajectories.
	 *
	 * @attention This is a very unefficient implementation.
	 * */
	continuous_set_ptr intersection_with(const math::lin_constraint_system<scalar_type> & cons) const;

//	const trajectory & get_roots() const;
//	void add_roots(const trajectory & t);
//	void set_roots(const trajectory & t);
//
//	const trajectory & get_traj() const;
//	void add_traj(const trajectory & t);
//	void set_traj(const trajectory & t);

	/** Returns the number of roots/trajectories in the set */
	unsigned int size() const;

	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;

	/** Insert a root with its trajectory
	 */
	void insert(const state& x, const trajectory& traj);

	const fhyperbox & get_hbox() const;
	void set_hbox(const fhyperbox & b);

	/* \} */

private:
	state_to_trajectory_map my_map;
	fhyperbox myhbox;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	void
	compute_support_impl(const math::vdom_vector<scalar_type>& l,
			scalar_type& max_value, math::vdom_vector<scalar_type>& support_vec,
			bool& is_empty, bool& is_bounded) const;
};

}

#include "continuous_set_simulation.hpp"

#endif
