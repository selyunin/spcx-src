#ifndef continuous_dynamics_COMPOSITION_H_
#define continuous_dynamics_COMPOSITION_H_

#include <stdexcept>
#include <sstream>
#include "utility/shared_ptr_output.h"
#include "utility/basic_warning.h"
#include "core/predicates/node_print_visitor.h"
#include "math/vdom/affine_map_utility.h"
#include "core/continuous/continuous_dynamics/continuous_dynamics.h"

namespace continuous {

// --------------------------------------------
/** \name Interface for parallel composition of continuous_dynamics
 *  \{ */
// --------------------------------------------

/** Computes the parallel composition of the two dynamics t1 and t2, and puts the result in t_ret.
 * Returns true iff the composition was successful. If unsuccessful, t_ret is the null pointer. */
bool parallel_compose(continuous_dynamics::ptr& t_ret,
		const continuous_dynamics::const_ptr& t1,
		const continuous_dynamics::const_ptr& t2);

/** Computes the parallel composition of the two dynamics t1 and t2, and puts the result in t_ret.
 * If unsuccessful, returns the null pointer.
 *
 * @note This is simply an alternative interface to the above. */
continuous_dynamics::ptr parallel_compose(
		const continuous_dynamics::const_ptr& t1,
		const continuous_dynamics::const_ptr& t2);

/* \} */
// --------------------------------------------
/** \name Default implementation
 *  \{ */
// --------------------------------------------

/** Default implementation of parallel composition. */
template<typename T1, typename T2> class continuous_dynamics_composition_operator {
public:
	static continuous_dynamics::ptr implement(const T1* t1, const T2* t2) {
		assert(t1);
		assert(t2);
		std::cerr << "offended by dynamics:" << std::endl;
		std::cerr << "t1:" << t1->get_predicate() << std::endl;
		std::cerr << "t2:" << t2->get_predicate() << std::endl;
		std::runtime_error(
				"Missing implementation; could not compose dynamics.");
		return continuous_dynamics::ptr(); // failed to find composition
	}
	;
};

/* \} */
// --------------------------------------------
/** \name Specializations for fully derived classes
 *  \{ */
// --------------------------------------------

/** Compose ode_affine_dynamics with ode_affine_dynamics.
 *
 * If the dynamics are conflicting (not the same coefficients
 * in at least one of the variables), the empty affine map is
 * returned. The intended meaning of this is that no
 * time elapse is to be performed. */
template<typename T> class continuous_dynamics_composition_operator<
		ode_affine_dynamics<T> , ode_affine_dynamics<T> > {
public:
	static continuous_dynamics::ptr implement(
			const ode_affine_dynamics<T> * t1,
			const ode_affine_dynamics<T> * t2) {
		continuous_dynamics::ptr p = continuous_dynamics::ptr();
		try {
			math::info_resolver<T> resolv;
			p = continuous_dynamics::ptr(new ode_affine_dynamics<T> (
					math::compose(*t1, *t2, resolv)));
			if (resolv.has_conflict()) {
				basic_warning("Composing dynamics", "Created empty dynamics",
						basic_warning::UNUSUAL_INPUT);
				positional_vdomain cdom = compose(t1->domain(), t2->domain());
				p = continuous_dynamics::ptr(new ode_affine_dynamics<T> (
						math::affine_map<T>::void_map(cdom)));
			}
		} catch (std::exception& e) {
			std::stringstream s;
			logger::copyfmt_to(s);
			s << *t1 << "\nand \n" << *t2 << std::endl;
			throw basic_exception("Could not compose the following dynamics:\n"
					+ s.str(), e);
		}
		return p;
	}
	;
};

/** Compose constant_bound_dynamics with constant_bound_dynamics. */
template<> class continuous_dynamics_composition_operator<
		constant_bound_dynamics, constant_bound_dynamics> {
public:
	static continuous_dynamics::ptr implement(
			const constant_bound_dynamics* t1,
			const constant_bound_dynamics* t2);
};

/** Compose relation_dynamics with relation_dynamics. */
template<> class continuous_dynamics_composition_operator<relation_dynamics,
		relation_dynamics> {
public:
	static continuous_dynamics::ptr implement(const relation_dynamics* t1,
			const relation_dynamics* t2);
};

/* \} */

}

#endif /*continuous_dynamics_COMPOSITION_H_*/
