/*
 * choose_directions.h
 *
 *  Created on: Nov 27, 2009
 *      Author: frehse
 */

#ifndef CHOOSE_DIRECTIONS_H_
#define CHOOSE_DIRECTIONS_H_

//#include "../math/matrix.h"
#include "math/vector.h"
#include "math/vdom/vdom_vector.h"
#include "core/continuous/polyhedra/polyhedron.h"
#include "math/uni_sphere.h"
#include "math/unique_vector_to_value_store.h"
#include "io/common_input/parse_policy.h"

//#include "../abstract_framework/continuous/polyhedron_utility.h"
//#include "../abstract_framework/continuous/continuous_dynamics/continuous_dynamics.h"
//#include "../abstract_framework/symbolic_states/symbolic_state_collection.h"
//#include "../set_implementations/constr_polyhedron/constr_polyhedron.h"
//#include "../set_implementations/ppl_polyhedron/convert_to_ppl.h"
//#include "../set_implementations/stl_set_implementations/discrete_set_stl_set.h"
//#include "../set_implementations/sfm_poly/sfm_cont_set.h"

namespace continuous {
namespace support_function {

/** Construct R direction vectors and insert them in the list.
 *
 * We construct vectors as the normal to the plane given by the constraints of polytope cset. This should be replaced
 * with more sensible choice of directions later in future.
 */
template<class T> void add_poly_directions(
		const continuous_set::const_ptr& cset,
		std::list<math::vector<T> >& directions) {
	math::vector<T> l;
	typename polyhedron<T>::const_ptr p = boost::static_pointer_cast<
			const polyhedron<T> >(cset);

	const typename math::lin_constraint_system<T>::const_ptr& cs =
			p->get_constraints();
	int i = 0;
	for (typename math::lin_constraint_system<T>::const_iterator it =
			cs->begin(); it != cs->end(); ++it, ++i) {
		l = it->get_l().get_vector();
		directions.push_back(l);
	}
}
;

/** Return the normal vectors of a hyperbox.
 *
 * The vectors are normed with the infinity norm. */
template<class T> void add_box_directions(unsigned int dim, std::list<
		math::vector<T> >& directions) {
	/**
	 * Chooses directions so as to compute a box .
	 */
	for (unsigned int i = 0; i < dim; i++) {
		math::vector<T> l = math::vector<T>(dim, T(0));
		l[i] = T(1);
		directions.push_back(l);
		l = math::vector<T>(dim, T(0));
		l[i] = T(-1);
		directions.push_back(l);
	}
}
;

/** Return the normal vectors of a bounding box.
 *
 * The vectors are normed with the infinity norm. */
template<class T> void add_octagonal_directions(unsigned int dim, std::list<
		math::vector<T> >& directions) {
	/**
	 * Chooses directions so as to compute the octagon of the set.
	 * An octagon has two nonzero coefficients, which leads to
	 * directions of the form:
	 * x_i \in {-1,+1}
	 * x_j \in {-1,+1}
	 * for i,j in 0,...,dim.
	 * Note that bounding box directions are contained as the
	 * case i=j.
	 */
	math::vector<T> l = math::vector<T>(dim, T(0));
	for (unsigned int i = 0; i < dim; i++) {
		for (unsigned int j = i; j < dim; j++) {
			l[i] = T(1);
			l[j] = T(1);
			directions.push_back(l); // 1 1
			l[i] = T(-1);
			directions.push_back(l); // -1 1
			if (i != j) { // otherwise, these are redundant
				l[j] = T(-1);
				directions.push_back(l); // -1 -1
				l[i] = T(1);
				directions.push_back(l); // 1 -1
			}
			l[j] = T(0);
		}
		l[i] = T(0);
	}
}
;

/** Return the normal vectors of a circle in the i,j subspace of a
 * dim-dimensional space at n sampling points. */
template<class T> void add_uniform_2D_directions(unsigned int dim,
		unsigned int i, unsigned int j, unsigned int n, std::list<math::vector<
				T> >& directions) {
	double theta = 2 * M_PI / n;
	for (unsigned int k = 0; k < n; k++) {
		math::vector<T> l = math::vector<T>(dim, T(0));
		l[i] = T(std::cos(k * theta));
		l[j] = T(std::sin(k * theta));
		directions.push_back(l);
	}
}
;

/** Return nb normal vectors that are uniformly distributed.
 *
 * The vectors are normed with the infinity norm. */
template<class T> void add_uniform_directions(unsigned int dim,
		unsigned int nb, std::list<math::vector<T> >& directions) {
	/*	if (nb == 2 * dim) {
	 add_box_directions<T>(cset, directions);
	 } else if (nb == 2 * dim * dim) {
	 add_octagonal_directions<T>(cset, directions);
	 } else if (dim == 2) {
	 add_uniform_2D_directions(2, 0, 1, nb, directions);
	 } else { */
	/** Define static cache */
	typedef std::vector<math::vector<double> > vector_vector;
	typedef std::pair<unsigned int, unsigned int> key_type;
	static std::map<key_type, vector_vector> vecs_cache;

	key_type key = std::make_pair(dim, nb);

	/** Check if directions are in cache */
	std::vector<math::vector<double> > vecs;
	std::map<key_type, vector_vector>::const_iterator cache_it =
			vecs_cache.find(key);
	if (cache_it == vecs_cache.end()) {
		/* If not in cache, compute directions */
		vecs = math::uni_sphere(nb, dim, 10000 * nb, 1e-3);

		/* Norm the vector. */
		for (std::vector<math::vector<double> >::iterator it =
				vecs.begin(); it != vecs.end(); ++it) {
			(*it) = (*it)/(it->infinity_norm());
		}
		vecs_cache[key] = vecs;
	} else {
		vecs = cache_it->second;
	}

	/** Add directions to return set */
	for (std::vector<math::vector<double> >::const_iterator it = vecs.begin(); it
			!= vecs.end(); ++it) {
		directions.push_back((*it).template convert_to<T> ());
	}
	//	}
}

class direction_chooser {
public:
	typedef enum {
		BOX, OCTAGONAL, POLYH, UNIFORM
	} direction_type;
	typedef std::set<direction_type> type_store;
	typedef global_types::float_type scalar_type;
	typedef math::vdom_vector<scalar_type> direction_vector;
	typedef math::unique_vector_to_value_store<scalar_type, math::vdom_vector,
			bool> direction_store; // bool is just used as a dummy

	direction_chooser() {
		my_number = 0;
	}

	static void reset() {
		my_types.clear();
		my_store.clear();
	}

	static const type_store& get_types() {
		return my_types;
	}
	static void set_types(const type_store& types, const direction_store& store = direction_store()) {
		my_types = types;
		my_store = store;
	}

	/** Adds the directions given by the string opt.
	 *
	 * If uniform directions are added several times,
	 * their number is added up.
	 *
	 * Look up variable symbols within context.
	 */
	static void add_directions(const std::string& opt, const std::string& context = "", const parser::parse_policy& ppol = parser::parse_policy());

	/** Get directions for a space of dimension dim.
	 *
	 * @attention Does not add polyhedral directions, since
	 * the domain can not be reproduced.
	 * Throws if polyehdral directions are activated. */
	template<typename T>
	static std::list<math::vector<T> > get_directions(unsigned int dim) {
		typedef math::unique_vector_to_value_store<T, math::vector, bool>
				result_direction_store;

		result_direction_store res_store; // for eliminating duplicates
		std::list<math::vector<T> > dirs; // the other routines
		for (type_store::const_iterator it = my_types.begin(); it
				!= my_types.end(); ++it) {
			unsigned int nb;
			switch (*it) {
			case BOX:
				add_box_directions(dim, dirs);
				break;
			case OCTAGONAL:
				add_octagonal_directions(dim, dirs);
				break;
			case UNIFORM:
				nb = my_number;
				if (nb == 0) {
					nb = (unsigned int) std::pow(double(dim), 3.0) * 2;
				}
				add_uniform_directions(dim, nb, dirs);
				break;
			case POLYH:
				throw basic_exception(
						"Polyhedral directions requested without domain.");
				break;
			default:
				throw basic_exception("unknown direction type "
						+ to_string(*it));
				break;
			}
		}
		// add dirs to res_store to eliminate doubles
		// use a default domain for the store
		positional_vdomain dom;
		for (typename std::list<math::vector<T> >::const_iterator it =
				dirs.begin(); it != dirs.end(); ++it) {
			res_store.insert(*it);
		}
		// produce the final list
		dirs.clear();
		for (typename result_direction_store::const_iterator it = res_store.begin(); it
				!= res_store.end(); ++it) {
			dirs.push_back(it->first);
		}
		return dirs;
	}
	;

	/** Get directions for the domain dom. */
	template<typename T>
	static std::list<math::vector<T> > get_directions(
			const positional_vdomain& dom) {
		typedef math::unique_vector_to_value_store<T, math::vector, bool>
				result_direction_store;

		unsigned int dim = dom.size();
		result_direction_store res_store; // for eliminating duplicates
		std::list<math::vector<T> > dirs; // the other routines
		for (type_store::const_iterator it = my_types.begin(); it
				!= my_types.end(); ++it) {
			unsigned int nb;
			switch (*it) {
			case BOX:
				add_box_directions(dim, dirs);
				break;
			case OCTAGONAL:
				add_octagonal_directions(dim, dirs);
				break;
			case UNIFORM:
				nb = my_number;
				if (nb == 0) {
					nb = (unsigned int) std::pow(double(dim), 3.0) * 2;
				}
				add_uniform_directions(dim, nb, dirs);
				break;
			case POLYH:
				// add directions from the direction_store
				for (direction_store::const_iterator dit = my_store.begin(); dit
						!= my_store.end(); ++dit) {
					math::vdom_vector<T> vec = dit->first.template convert_to<T> ();
					vec.remap(dom);
					res_store.insert(vec.get_vector());
				}
				break;
			default:
				throw basic_exception("unknown direction type "
						+ to_string(*it));
				break;
			}
		}
		//LOGGER(DEBUG4,"choose_directions","obtained " +  to_string(res_store.size()) + " directions for requested domain.");
		// add dirs to res_store to eliminate doubles
		for (typename std::list<math::vector<T> >::const_iterator it =
				dirs.begin(); it != dirs.end(); ++it) {
			res_store.insert(*it);
		}
		// produce the final list
		dirs.clear();
		for (typename result_direction_store::const_iterator it = res_store.begin(); it
				!= res_store.end(); ++it) {
			dirs.push_back(it->first);
		}
		LOGGER(DEBUG4,"choose_directions","obtained " +  to_string(res_store.size()) + " non-redundant directions for requested domain.");
		return dirs;
	}
	;
private:
	static type_store my_types;
	static direction_store my_store;
	static unsigned int my_number;
	// forbid copying and assigning
	direction_chooser(const direction_chooser&);
	direction_chooser& operator =(const direction_chooser&);
};

/** Choose directions for the support function hyperplanes. */
template<class T> void choose_directions(const positional_vdomain& dom,
		std::list<math::vector<T> >& directions) {
	directions = direction_chooser::get_directions<T>(dom);
}
;

/** Choose directions for the support function hyperplanes. */
template<class T> void choose_directions(continuous_set::const_ptr cset,
		std::list<math::vector<T> >& directions) {
	unsigned int dim = cset->get_dim();
	directions = direction_chooser::get_directions<T>(dim);
}
;

/** Choose directions for the support function hyperplanes. */
template<class T> void choose_directions(unsigned int dim, std::list<
		math::vector<T> >& directions) {
	directions = direction_chooser::get_directions<T>(dim);
}
;

}
}

#endif /* CHOOSE_DIRECTIONS_H_ */
