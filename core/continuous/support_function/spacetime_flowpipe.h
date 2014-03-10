/*
 * spacetime_flowpipe.h
 *
 *  Created on: Nov 3, 2012
 *      Author: notroot
 */

#ifndef SPACETIME_FLOWPIPE_H_
#define SPACETIME_FLOWPIPE_H_

#include <boost/enable_shared_from_this.hpp>
#include "core/continuous/continuous_set.h"
#include "math/vdom/index_to_variable_id_map_provider.h"
#include "core/post_operators/spacetime_post/spacetime_plif.h"

/** Forward declarations */
namespace spacetime {
template<typename scalar_type>
struct affine_support_flowpipe_description;
}

namespace continuous {

/** Representation of a flowpipe (subset) mapped and bloated.
 *
 * @attention The map and bloating are not included in evolutions. Directions are pre-map.
 *
 * @todo Fix the semantics - either move map and bloating elsewhere or take them into account when restricting domains.
 *  */
template<typename scalar_type>
class spacetime_flowpipe: public continuous_set,
		public index_to_variable_id_map_provider {
public:
	typedef boost::shared_ptr<spacetime_flowpipe<scalar_type> > ptr;
	typedef boost::shared_ptr<const spacetime_flowpipe<scalar_type> > const_ptr;
	typedef spacetime::affine_support_flowpipe_description<scalar_type> flowpipe_description;
	typedef spacetime::spacetime_plif spacetime_plif;
	typedef spacetime_plif::cut_point_method cut_point_method;
	typedef spacetime_plif::time_interval time_interval;
	typedef spacetime_plif::error_type error_type;
	typedef spacetime_plif::direction direction;
	typedef spacetime_plif::annotated_plif annotated_plif;
	typedef math::affine_map<scalar_type> affine_map;
	typedef boost::shared_ptr<affine_map> affine_map_ptr;

	/** Construct from an affine_support_flowpipe_description over a given time domain */
	spacetime_flowpipe(const flowpipe_description& asfd,
			const time_interval& tdom);

	/** Deep copy */
	spacetime_flowpipe(const spacetime_flowpipe<scalar_type>& flp);

	/** Virtual destructor */
	virtual ~spacetime_flowpipe();

	/** Creates a shallow copy of *this. */
	virtual spacetime_flowpipe<scalar_type>* clone() const;

	virtual spacetime_flowpipe<scalar_type>* create_universe() const;
	virtual spacetime_flowpipe<scalar_type>* create_empty() const;

	virtual int get_memory() const;

	virtual unsigned int get_dim() const;

	virtual math::tribool is_empty() const;

	virtual void embed_variables(const variable_id_set& id_set);
	virtual void existentially_quantify_variables(
			const variable_id_set& id_set);

	/** Define the map */
	void set_map(const affine_map& M);

	/** Define the bloating */
	void set_bloating(const support_function_provider::const_ptr& W);

	/** Get the map */
	const affine_map_ptr& get_map() const;

	/** Get the initial states */
	support_function_provider::const_ptr get_initial_states() const;

	/** Get the bloating */
	const support_function_provider::const_ptr& get_bloating() const;

	/** Get the time domain over which all evolutions are defined */
	const time_interval& get_time_domain() const;

	/** Add a constraint that will be added to the outer approximation */
	virtual void add_outer_constraint(
			const math::lin_constraint<scalar_type>& con, bool restrict_evolutions = false);

	/** Add a constraint that will be added to the outer approximation */
	virtual void add_outer_constraints(
			const math::lin_constraint_system<scalar_type>& cons, bool restrict_evolutions = false);

	/** Get or compute support evolution for given direction and error bound
	 *
	 * The evolution is cached for later retrieval. The direction is normed according to the norm defined by
	 * get_norm(). The evolution for an arbitrary (nonzero) direction d is therefore obtained with
	 * evo=get_or_compute_evolution(d,err)*get_norm(d);
	 * */
	const spacetime_plif::annotated_plif& get_or_compute_evolution(
			const direction& d, const error_type& err);

	/** Get a support evolution for given direction if already in cache, and null otherwise
	 *
	 * */
	const spacetime_plif::annotated_plif* get_evolution(const direction& d);

	/** Restrict the upper bound of the domain so that for all times
	 * the constraints are maybe satisfied.
	 *
	 * A constraint is "maybe satisfied" at a time instant if there is a point in the set
	 * that satisfies the constraint.
	 *
	 * This is intended for intersection with an invariant.
	 * The result can be empty. */
	virtual void restrict_to_maybe_sat_prefix(
			const math::lin_constraint_system<scalar_type>& cons,
			const error_type& err);

	/** Restrict the upper bound of the domain so that for all times
	 * the constraints are definitely satisfied.
	 *
	 * A constraint is "definitely satisfied" at a time instant if all points in the set
	 * satisfy the constraint.
	 *
	 * This is intended for intersection with an invariant.
	 * The result can be empty. */
	virtual void restrict_to_definitely_sat_prefix(
			const math::lin_constraint_system<scalar_type>& cons,
			const error_type& err);

	/** Get the subdomains in which the constraint is maybe satisfied
	 *
	 * A constraint is "maybe satisfied" by an interval if there is at least one point in the
	 * interval that satisfies the constraint.
	 * The subdomains are empty iff an empty vector is returned.
	 * */
	std::vector<time_interval> get_maybe_satisfying_subdomains(
			const math::lin_constraint<scalar_type>& con,
			const error_type& err, bool widen_by_error);

	/** Get the subdomains in which all the constraints in the system are maybe satisfied at the same time
	 *
	 * A constraint is "maybe satisfied" by an interval if there is at least one point in the
	 * interval that satisfies the constraint.
	 * The subdomains are empty iff an empty vector is returned.
	 * */
	std::vector<time_interval> get_maybe_satisfying_subdomains(
			const math::lin_constraint_system<scalar_type>& cons,
			const error_type& err, bool widen_by_error);

	/** Restrict the domain of *this (including cached evolutions) to a subdomain
	 * */
	void restrict_to_subdomain(const time_interval& intv);

	/** Obtain a set of convex flowpipe that cover *this */
	std::vector<spacetime_flowpipe<scalar_type> > convexify(const cut_point_method& m) const;

	/** Get outer polyhedral approximation with directional error less than err
	 *
	 * If no evolutions have yet been computed, a universe polyhedron is returned. */
	polyhedron_collection<scalar_type> compute_outer_polyhedra(
			const cut_point_method& m) const;

	/** Get outer polyhedral approximation as a single convex polyhedron
	 *
	 * If no evolutions have yet been computed, a universe polyhedron is returned. */
	typename polyhedron<scalar_type>::ptr compute_convex_outer_polyhedron(error_type err) const;

	/** Returns whether the outer approximation contains the flowpipe f
	 *
	 * If refine is true, then f is refined until the property can be decided.
	 */
	math::tribool decide_outer_contains(spacetime_flowpipe<scalar_type>& f,
			bool refine) const;

	/** Returns whether the outer approximation contains the outer approximation of the flowpipe f
	 */
	math::tribool decide_outer_contains(const spacetime_flowpipe<scalar_type>& f) const;

	/** Get the states at the end of the time domain
	 *
	 * Last states are neither mapped nor bloated. */
	support_function_provider::const_ptr get_last_states() const;

	/** Get the time variable of the set */
	const variable& get_time_variable() const;

	/** Returns the set in predicate form.*/
	virtual continuous_set_predicate::ptr get_predicate() const;

	/** Returns the flowpipe description */
	const flowpipe_description& get_flowpipe_description() const;

	/** Accept a visitor. */
	virtual void accept(
			dispatching::dispatcher<continuous_set_typelist>& d) const;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const;

	/** Filename for temp output
	 *
	 * @todo this is a temporary hack */
	static std::string temp_filename;

	/** Get flag for blocking refinement.
	 *
	 * If true, no further refinement of *this.
	 * This is useful for checking containment.
	 */
	bool is_refinement_locked() const;

	/** Block any further refinement
	 */
	void set_refinement_locked();

protected:
	/** Forbidden */
	spacetime_flowpipe();


	/** Construct for given splif, map, W, lock, and outer constraints cons */
	spacetime_flowpipe(spacetime_plif splif, affine_map_ptr map,
			support_function_provider::const_ptr W, bool refinement_locked, math::lin_constraint_system<scalar_type> cons);

private:
	spacetime_plif my_splif;
	affine_map_ptr my_map;
	support_function_provider::const_ptr my_W;
	bool my_refinement_locked;

public:
	/** Activates or deactivates simplification in evolution computations */
	static bool simplify_concave;
	static bool simplify_convex;

protected:
	math::lin_constraint_system<scalar_type> my_outer_constraints; // constraints that are added to the outer polytopes
};

template<typename s> std::string spacetime_flowpipe<s>::temp_filename;
template<typename s> bool spacetime_flowpipe<s>::simplify_concave(true);
template<typename s> bool spacetime_flowpipe<s>::simplify_convex(true);
}

#include "spacetime_flowpipe.hpp"

#endif /* SPACETIME_FLOWPIPE_H_ */
