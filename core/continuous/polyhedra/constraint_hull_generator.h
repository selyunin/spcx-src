/*
 * constraint_hull_generator.h
 *
 *  Created on: Oct 25, 2010
 *      Author: frehse
 */

#ifndef CONSTRAINT_HULL_GENERATOR_H_
#define CONSTRAINT_HULL_GENERATOR_H_

#include <iterator>

#include "polyhedron_collection.h"
#include "constr_polyhedron/constr_polyhedron.h"
#include "constr_polyhedron/constr_polyhedron_operators.h"

#include "math/numeric/comp.h"
#include "math/unique_vector_to_value_store.h"
#include "math/vdom/vdom_vector_operators.h"
#include "utility/logger_stopwatch.h"

namespace continuous {

/** Computes the constraint hull of a collection of polyhedra.
 *
 * The constraint hull is the set of constraints that are
 * satisfied by all polyhedra in the collection.
 *
 * If all polyhedra are template polyhedra over the same
 * set of directions, the constraint hull is identical
 * to the template hull.
 *
 * The normal vectors are set to infinity normal form:
 * They are divided by their infinity norm.
 *
 * @attention All polyhedra P_i must have non-redundant constraints,
 * i.e., for all constraints a.x<=b it must hold
 * b == max_{x \in P_i} a.x.
 * @attention Converts all constraints to non-strict inequalities.
 * */

template<typename scalar_type>
class constraint_hull_generator {
public:
	typedef math::unique_vector_to_value_store<scalar_type, math::vector,
			scalar_type> normal_vec_store;

	constraint_hull_generator(const polyhedron_collection<scalar_type>& p) :
		polys(p) {
		// get a domain over all poly (otherwise we couldn't compare vectors)
		dom = positional_vdomain(polys.get_variable_ids());
	}
	;

private:
	/** Create a store of all directions with their max. b. */
	void construct_stores() {
		assert(polys.size()>1);
		assert(!polys.is_empty());

		// create the store
		max_store = normal_vec_store();
		min_stores = std::vector<normal_vec_store>(polys.size());

		unsigned int poly_count = 0;
		// process the polyhedra one by one,
		// constructing a global store with all directions
		// and local ones with non-redundant directions.
		// the local ones will later also be used for fast look-up.
		for (typename polyhedron_collection<scalar_type>::const_iterator pt =
				polys.begin(); pt != polys.end(); ++pt) {

			// get the constraints in the poly in canonic form on the same domain
			typename math::lin_constraint_system<scalar_type>::const_ptr cons =
					(*pt)->get_constraints();

			for (typename math::lin_constraint_system<scalar_type>::const_iterator
					con_it = cons->begin(); con_it != cons->end(); ++con_it) {
				// process each constraint one by one
				math::vdom_vector<scalar_type> adom = con_it->get_normal();
				// put a into the global domain
				adom.reorder(dom);
				math::vector<scalar_type> a = adom.get_vector();

				scalar_type b = -con_it->get_canonic_inh_coeff();
				// make a normed
				scalar_type f = a.infinity_norm();
				a /= f;
				b /= f;
				// update the global store with the max
				max_store.template update<std::less<scalar_type> > (a, b);
				//max_store.print(std::cout); std::cout << std::endl;
				// update the local stores with the min (to eliminate redundant constraints)
				min_stores[poly_count].template update<math::numeric::more<
						scalar_type> > (a, b);
				//min_stores[poly_count].print(std::cout); std::cout << std::endl;
				// if it's an equality, update also with the neg constraint
				if (con_it->is_equality()) {
					max_store.template update<std::less<scalar_type> > (-a, -b);
					min_stores[poly_count].template update<math::numeric::more<
							scalar_type> > (-a, -b);
				}
			}
			++poly_count;
		}
	}
	;

	/** Find the support in a given direction in poly k.
	 */
	scalar_type get_or_compute_support(const math::vector<scalar_type>& vec,
			unsigned int k) {
		scalar_type local_b;
		// find the normal vector
		typename normal_vec_store::const_iterator jt = min_stores[k].find(vec);
		// if it's not there, compute the support
		if (jt == min_stores[k].end()) {
			// use double, since I'm lazy
			math::vdom_vector<double>
					l(dom, vec.template convert_to<double> ());
			double max_value;
			math::vdom_vector<double> support_vec;
			bool is_empty;
			bool is_bounded;
			// get an iterator to the k-th poly
			typename polyhedron_collection<scalar_type>::const_iterator pt =
					polys.begin();
			std::advance(pt, k);
			(*pt)->compute_support(l, max_value, support_vec, is_empty,
					is_bounded);
			if (is_empty) {
				throw basic_exception(
						"constraint_hull_generator: can't deal with polys");
			}
			if	(!is_bounded){
				std::stringstream ss;
				ss << l;
				throw basic_exception("constraint_hull_generator: can't deal with unbounded polys in the direction:"+ss.str());
			}
			local_b = convert_element<scalar_type> (max_value);

			// add the new constraint also the local store
			min_stores[k].insert(vec, local_b);
		} else {
			// otherwise, we're already ok
			local_b = jt->second;
		}
		return local_b;
	}
	;

	/** If a poly doesn't direction vec in the store, compute it using the support
	 * and return the max. */
	scalar_type complete_local_stores(
			const math::vector<scalar_type>& global_vec) {
		assert(polys.size()>0);

		scalar_type global_b;

		unsigned int poly_count = 0;
		for (typename polyhedron_collection<scalar_type>::const_iterator pt =
				polys.begin(); pt != polys.end(); ++pt) {
			// find the normal vector
			// if it's not there, compute the support
			scalar_type local_b =
					get_or_compute_support(global_vec, poly_count);

			// take the max
			if (poly_count == 0 || global_b < local_b) {
				global_b = local_b;
			}

			++poly_count;
		}

		return global_b;
	}
	;

	/** If a poly doesn't have a direction in the store, compute it using the support. */
	void complete_local_stores() {
		// for each poly, check if all global constraints are present.
		// if not, compute the support and set the global constraint to the max
		for (typename normal_vec_store::iterator it = max_store.begin(); it
				!= max_store.end(); ++it) {
			const math::vector<scalar_type>& global_vec = it->first;
			scalar_type& global_b = it->second;

			scalar_type local_b = complete_local_stores(global_vec);
			if (global_b < local_b) {
				global_b = local_b;
			}
		}
	}
	;

	/** Construct an initial (universe) hull and treat trivial
	 * cases (empty or singular set).
	 *
	 * Set done to true if the result needs no further processing,
	 * i.e., can be returned as final result.
	 */
	typename polyhedron<scalar_type>::ptr get_initial_hull(bool& done) {
		// start with a universe poly
		typename polyhedron<scalar_type>::ptr res =
				typename polyhedron<scalar_type>::ptr(
						polys.default_element()->create_universe());
		done = false;
		if (polys.size() == 1) {
			res = typename polyhedron<scalar_type>::ptr(
					(*polys.begin())->clone());
			done = true;
		} else if (polys.is_empty()) {
			res = typename polyhedron<scalar_type>::ptr(
					polys.default_element()->create_empty());
			done = true;
		}
		return res;
	}
	;

	/** Convert a store to a polyhedron */
	typename polyhedron<scalar_type>::ptr convert_to_poly(
			const normal_vec_store& store) {
		// start with a universe poly
		typename polyhedron<scalar_type>::ptr res =
				typename polyhedron<scalar_type>::ptr(
						polys.default_element()->create_universe());
		// reconvert the global store to a polyhedron
		for (typename normal_vec_store::const_iterator it = store.begin(); it
				!= store.end(); ++it) {
			// because the con is defined as a.x+b <= 0, we have
			// to switch the sign of b
			math::lin_constraint<scalar_type> con(
					math::vdom_vector<scalar_type>(dom, it->first),
					-it->second, LE);
			res->add_constraint(con);
		}
		return res;
	}
	;

public:
	/** Compute the constraint hull of all polyhedra.
	 */
	typename polyhedron<scalar_type>::ptr compute_constraint_hull() {
		bool done;
		typename polyhedron<scalar_type>::ptr res = get_initial_hull(done);

		if (!done) {
			// construct a store for each polyhedron and
			// add the max in every direction to the max_store
			construct_stores();
//						std::cout << "first max:";
//						max_store.print(std::cout);
//						std::cout << std::endl;

			// add the missing directions to each polyhedron,
			// updating the global max if necessary
			complete_local_stores();
//						std::cout << "final max:";
//						max_store.print(std::cout);
//						std::cout << std::endl;

			// reconvert the global store to a polyhedron
			res = convert_to_poly(max_store);
//			std::cout << "final poly:";
//			std::cout << res << std::endl;
		}
		return res;
	}
	;

	/** Compute the constraint hull of all polyhedra up to a relative error err.
	 *
	 * Starts a new polytope each time the relative error is hit.
	 * @note An absolute error of 0 is nonsense since this results in no clustering
	 * as soon as the polyhedra involve the same equality.
	 */
	polyhedron_collection<scalar_type> compute_constraint_hull_up_to_error(
			double err_rel, double err_abs = 0) {
		using namespace math;
		using namespace numeric;

		bool done;
		polyhedron_collection<scalar_type> res;
		typename polyhedron<scalar_type>::ptr res_poly = get_initial_hull(done);

		if (done) {
			res.insert(res_poly);
			return res;
		} else {
			construct_stores();

			// add opposite directions
			for (typename normal_vec_store::iterator it = max_store.begin(); it
					!= max_store.end(); ++it) {
				math::vector<scalar_type> neg_vec = -it->first;
				if (max_store.find(neg_vec) == max_store.end()) {
					scalar_type neg_b = complete_local_stores(neg_vec);
					max_store.insert(neg_vec, neg_b);
				}
			}

			complete_local_stores();

			// include one poly at a time until error bound is
			// exceeded. Then start a new store
			unsigned int poly_count = 0;
			normal_vec_store store = min_stores[poly_count];
			typename polyhedron_collection<scalar_type>::const_iterator pt =
					polys.begin();
			++pt; // start with the second poly
			poly_count = 1;
			for (; pt != polys.end(); ++pt) {
				// check if the error is greater than the bound
				bool too_nonconvex = false;
				for (typename normal_vec_store::iterator it = store.begin(); it
						!= store.end() && !too_nonconvex; ++it) {
					math::vector<scalar_type> vec = it->first;
					scalar_type b_current = it->second;
					scalar_type b_max = max_store.get_value(vec);
					scalar_type b_min = -max_store.get_value(-vec);
					scalar_type b_poly =
							get_or_compute_support(vec, poly_count);
//					std::cout << "vec:" << vec << " i:" << poly_count << " b:"
//							<< b_poly << std::endl;

					// test if abs(b_max-b_poly) >= (1+err_rel)*(b_max-b_min)+err_abs
					scalar_type err = b_current - b_poly;
					scalar_type max_err = scalar_type(err_rel)
							* (b_max - b_min) + scalar_type(err_abs);
					if (definitely(is_GT(abs(err),max_err))) {
						too_nonconvex = true;
					}
					//					std::cout << "b_err: " << err << " vs " << convert_element<
					//							double> (max_err) << std::endl;
				}
				// if not, expand the current store with all the constraints in the poly
				if (!too_nonconvex) {
					for (typename normal_vec_store::iterator it = store.begin(); it
							!= store.end(); ++it) {
						math::vector<scalar_type> vec = it->first;
						scalar_type b_poly = get_or_compute_support(vec,
								poly_count);
						store.template update<std::less<scalar_type> > (vec,
								b_poly);
					}
				} else {
					// add the current store as a poly to the result set
					// and start a new store with the current poly
					res.insert(convert_to_poly(store));
					store = min_stores[poly_count];
				}
				++poly_count;
			}
			res.insert(convert_to_poly(store));

			return res;
		}
	}
	;

private:
	const polyhedron_collection<scalar_type>& polys;
	normal_vec_store max_store;
	std::vector<normal_vec_store> min_stores;
	// domain over all poly
	positional_vdomain dom;
};

}

#endif /* CONSTRAINT_HULL_GENERATOR_H_ */
