/*
 * INTV_formatter.h
 *
 */

#ifndef INTV_FORMATTER_H_
#define INTV_FORMATTER_H_

#include <stdexcept>
#include <typeinfo>
#include "utility/basic_warning.h"
#include "io/common_output/dispatch_output_formatter.h"
#include "math/numeric/interval.h"

#include "core/continuous/support_function_provider_upcaster.h"

#include "core/continuous/support_function/sfm/sfm_cont_set.h"


/** Forward declarations */
namespace continuous {
template<typename T> class polyhedron;
template<typename T> class spacetime_flowpipe;
}

namespace io {

class INTV_formatter;

template<typename T>
class INTV_continuous_set_formatter {
public:
	static void output(output_formatter& of, const T& c) {
		std::string tname = typeid(&c).name();
		throw std::runtime_error(
				"INTV_continuous_set_formatter : can't format object pointed to by "
						+ tname + ".");
	}
	;
};

template<>
class INTV_continuous_set_formatter<continuous::support_function_provider> {
public:
	static void output(output_formatter& of,
			const continuous::support_function_provider& c);
};

template<typename scalar_type>
class INTV_continuous_set_formatter<continuous::support_function::sfm_cont_set<
		scalar_type> > {
public:
	static void output(output_formatter& of,
			const continuous::support_function::sfm_cont_set<scalar_type> & c) {
		//		continuous::support_function::sfm_cont_set<scalar_type>::set_output_format(
		//				continuous::support_function::sfm_cont_set<scalar_type>::DOUBLE_GENERATORS);
		//		c.print(os);
		continuous::polyhedron_collection<scalar_type> coll=c.get_outer_polytope_collection();
		of.output(coll);
//		if (!c.is_empty()) {
//			for (size_t i = 0; i < c.get_size(); ++i) {
//				continuous::constr_polyhedron<scalar_type> p =
//						c.get_polytope(i);
//				if (!p.is_empty()) {
//					of.output(p);
//					of.get_os() << std::endl << std::endl;
//				}
//			}
//		} else {
//			// generate a warning that the set is empty
//			basic_warning("GEN format output", "cannot represent empty set in GEN format",
//					basic_warning::AMBIGUOUS_OUTPUT);
//		}
	}
	;
};

/** Interval format.
 *
 * Outputs the global bounds [min,max] on each system variable and the
 * bounds per location
 *
 * Tries to first upcast to a polyhedron, then to a sfm_cont_set,
 * then to a support_function_provider.
 */
template<typename T> class INTV_formatter_caster: public dispatching::caster_list<
		continuous::support_function_provider_upcaster<T>,
		dispatching::convertible_caster<T,
					continuous::support_function::sfm_cont_set<
							global_types::float_type> > > {
};

class INTV_formatter: public dispatch_output_formatter<
		INTV_continuous_set_formatter, INTV_formatter_caster> {
	//		GEN_continuous_set_formatter, castercontinuous::polyhedron_upcaster> {
public:
	typedef dispatch_output_formatter<INTV_continuous_set_formatter,
			INTV_formatter_caster> base_class;

	typedef math::numeric::interval<double> bound_interval;

	/** Create an INTV_formatter that sends its output to the stream os. */
	explicit INTV_formatter(std::ostream& os) :
		base_class(os) {
	}
	;

	using base_class::output;
	/** Output a symbolic state.
	 *
	 * Same as for the base class, but issues a warning if discrete set is empty.
	 */
	void output(const hybrid_automata::symbolic_state& sstate);

	/** Output a symbolic state collection */
	void output(const hybrid_automata::symbolic_state_collection& sstates);

	/** Output a spacetime_flowpipe
	 *
	 * Computes the max in the axis direction outputs that. */
	void output(
			const continuous::spacetime_flowpipe<
					global_types::float_type>& flowpipe);

	/** Output a continuous set */
	void output(const continuous::continuous_set& c);
};

/*----------------------------------
 Implementations
 ----------------------------------
 */


}

#endif /* INTV_FORMATTER_H_ */
