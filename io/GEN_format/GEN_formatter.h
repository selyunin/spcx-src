/*
 * GEN_formatter.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef GEN_FORMATTER_H_
#define GEN_FORMATTER_H_

#include <stdexcept>
#include <typeinfo>
#include "utility/basic_warning.h"
#include "io/common_output/dispatch_output_formatter.h"

// Sets can't be incomplete types
#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/continuous/polyhedra/polyhedron_upcaster.h"
#include "core/continuous/support_function_provider_upcaster.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"
#include "core/continuous/support_function/spacetime_flowpipe.h"

#include "core/continuous/support_function/approximation/vertice_approx_2D_sv.h"
#include "core/continuous/support_function/approximation/vertice_approx_2D.h"

#include "math/vdom/trajectory.h"

#include "math/ode_solving/traj_simu/continuous_set_simulation.h"
#include "io/INTV_format/INTV_formatter.h"

/** Forward declarations */
namespace continuous {
template<typename T> class polyhedron;
}

namespace io {

class GEN_formatter;

template<typename T>
class GEN_continuous_set_formatter {
public:
	static void output(output_formatter& of, const T& c) {
		std::string tname = typeid(&c).name();
		throw std::runtime_error(
				"GEN_continuous_set_formatter : can't format object pointed to by "
						+ tname + ".");
	}
	;
};

template<>
class GEN_continuous_set_formatter<continuous::support_function_provider> {
public:
	static void output(output_formatter& of,
			const continuous::support_function_provider& c);
};

template<typename scalar_type>
class GEN_continuous_set_formatter<continuous::polyhedron<scalar_type> > {
public:
	static void output(output_formatter& of,
			const continuous::polyhedron<scalar_type> & c);
};

template<typename scalar_type>
class GEN_continuous_set_formatter<continuous::support_function::sfm_cont_set<
		scalar_type> > {
public:
	static void output(output_formatter& of,
			const continuous::support_function::sfm_cont_set<scalar_type> & c) {
		//		continuous::support_function::sfm_cont_set<scalar_type>::set_output_format(
		//				continuous::support_function::sfm_cont_set<scalar_type>::DOUBLE_GENERATORS);
		//		c.print(os);
		//continuous_set_collection coll=c.get_outer_polytop_collection(i);
		//of.output(coll);
		if (!math::definitely(c.is_empty())) {
			for (size_t i = 0; i < c.get_size(); ++i) {
				continuous::constr_polyhedron<scalar_type> p =
						c.get_polytope(i);
				if (!math::definitely(p.is_empty())) {
					of.output(p);
					of.get_os() << std::endl << std::endl;
				}
			}
		} else {
			// generate a warning that the set is empty
			basic_warning("GEN format output",
					"cannot represent empty set in GEN format",
					basic_warning::AMBIGUOUS_OUTPUT);
		}
	}
	;
};

/** Generator format.
 *
 * Outputs a list of vertices for each continuous_set, overapproximated
 * up to eps>=0 accuracy.
 *
 * Tries to first upcast to a polyhedron, then to a sfm_cont_set,
 * then to a support_function_provider.
 */
template<typename T> class GEN_formatter_caster: public dispatching::caster_list<
		continuous::support_function_provider_upcaster<T>,
		dispatching::caster_list<continuous::polyhedron_upcaster<T>,
				dispatching::convertible_caster<T,
						continuous::support_function::sfm_cont_set<
								global_types::float_type> > > > {
};

class GEN_formatter: public dispatch_output_formatter<
		GEN_continuous_set_formatter, GEN_formatter_caster> {
	//		GEN_continuous_set_formatter, castercontinuous::polyhedron_upcaster> {
public:
	typedef dispatch_output_formatter<GEN_continuous_set_formatter,
			GEN_formatter_caster> base_class;

	/** Construct generator formatter for a given stream, error bound */
	GEN_formatter(std::ostream& os, double eps) :
		base_class(os), my_eps(eps) {
	}
	;
	using base_class::output;

	/** Output a symbolic state collection.
	 *
	 * Same as for the base class, but issues a warning if discrete set is empty.
	 */
	void output(const hybrid_automata::symbolic_state_collection& sstates);

	/** Output a symbolic state.
	 *
	 * Same as for the base class, but issues a warning if discrete set is empty.
	 */
	void output(const hybrid_automata::symbolic_state& sstate);

	/** Output a trajectory in the state space.
	 *
	 * If the trajectory does not have one of the output variables, zeros are output instead.
	 * If none of the output variables are present, then nothing is output.
	 * */
	virtual void output(const math::trajectory<global_types::float_type>& traj);

	/** Output a continuous_set_simulation. (outputs the trajectory part, and not the root part) */
	virtual void output(
			const continuous::continuous_set_simulation<
					global_types::float_type>& cset_simu);

	/** Output a polyhedron_collection */
	void output(
			const continuous::polyhedron_collection<
					global_types::float_type>& polys);

	/** Output a spacetime_flowpipe
	 *
	 * Computes the outer polyhedral approximation and outputs that. */
	void output(
			const continuous::spacetime_flowpipe<
					global_types::float_type>& flowpipe);

	void output(const continuous::continuous_set& c) {
		// dynamic cast to continuous set classes that
		// are not in the typelist

		if (const continuous::continuous_set_simulation<global_types::float_type> * p =
				dynamic_cast<const continuous::continuous_set_simulation<global_types::float_type> *>(&c)) {
			output(*p);
		} else if (const continuous::polyhedron_collection<global_types::float_type> * p =
				dynamic_cast<const continuous::polyhedron_collection<global_types::float_type> *>(&c)) {
			output(*p);
		}  else if (const continuous::spacetime_flowpipe<global_types::float_type> * p =
				dynamic_cast<const continuous::spacetime_flowpipe<global_types::float_type> *>(&c)) {
			output(*p);
		}else {
			// otherwise call base class output to try double dispatch
			base_class::output(c);
		}
	}
	;

	double my_eps;
};

class TRAJ_formatter: public GEN_formatter {
public:
	typedef GEN_formatter base_class;

	TRAJ_formatter(std::ostream& os) :
		base_class(os,0.) {
	}
	;

	/** Output a trajectory in the state space using traj format*/
	virtual void output(const math::trajectory<global_types::float_type>& traj) {
		get_os() << traj;
	};

	/** Output a continuous_set_simulation. (outputs the trajectory part, and not the root part) */
	virtual void output(const continuous::continuous_set_simulation<global_types::float_type>& cset_simu){
		std::cout << "Printing trajectory format TRAJ" << std::endl;
		for (continuous::continuous_set_simulation<global_types::float_type>::const_iterator
				it = cset_simu.begin(); it != cset_simu.end(); ++it) {
			output(it->second);
		}
	};

};






/*----------------------------------
 Implementations
 ----------------------------------
 */

template<typename scalar_type>
void GEN_continuous_set_formatter<continuous::polyhedron<scalar_type> >::output(
		output_formatter& of, const continuous::polyhedron<scalar_type> & c) {
	/*
	 if (of.get_output_variables().empty()) {
	 c.print_double_generators(of.get_os());
	 } else {
	 typename continuous::polyhedron<scalar_type>::ptr cset(c.clone());
	 cset->project_to_variables(of.get_output_variables());
	 cset->print_double_generators(of.get_os());
	 } */

	if (!math::definitely(c.is_empty())) {
		variable_id x_id = 0;
		variable_id y_id = 0;
		bool found = false;
		if (of.get_output_variables().empty()) {
			variable_id_set vis = c.get_variable_ids();
			if (vis.size() == 2) {
				x_id = *vis.begin();
				y_id = *vis.rbegin();
				found = true;
			}
		} else if (of.get_output_variables().size() == 2) {
			x_id = *of.get_output_variables().begin();
			y_id = *of.get_output_variables().rbegin();
			found = true;
		}
		if (!found) {
			c.print_double_generators(of.get_os());
		} else {
			//std::cout << "approx:";
			double eps = static_cast<GEN_formatter&> (of).my_eps;

			try {
				typename continuous::support_function::vertice_approx_2D_sv<
						scalar_type>::sample_list pts;
				pts = continuous::support_function::vertice_approx_2D_sv<
						scalar_type>::approx(c, x_id, y_id, scalar_type(eps));
				//std::cout << pts.size() << std::endl;
				for (typename continuous::support_function::vertice_approx_2D_sv<
						scalar_type>::sample_list::const_iterator it =
						pts.begin(); it != pts.end(); ++it) {
					of.get_os() << convert_element<double> (it->a) << " "
							<< convert_element<double> (it->b) << std::endl;
				}
				// repeat the first point to close the curve
				of.get_os() << convert_element<double> (pts.begin()->a) << " "
						<< convert_element<double> (pts.begin()->b)
						<< std::endl;

			} catch (continuous::unbounded_set_exception& e) {
				// prepare a user-readable output of the variable bounds
				std::stringstream ss;
				INTV_formatter intvf(ss);
				intvf.set_output_variables(of.get_output_variables());
				intvf.set_context(of.get_context());
				intvf.output(c);
				// generate a warning that the set is unbounded
				basic_warning(
						"GEN format output",
						"cannot represent unbounded set in GEN format. Bounds:\n"
								+ ss.str(), basic_warning::AMBIGUOUS_OUTPUT);

			} catch (continuous::empty_set_exception& e) {
				// generate a warning that the set is empty
				basic_warning("GEN format output",
						"cannot represent empty set in GEN format",
						basic_warning::AMBIGUOUS_OUTPUT);
			}
		}
		// terminate poly with an empty line
		of.get_os() << std::endl;
	} else {
		// generate a warning that the set is empty
		basic_warning("GEN format output",
				"cannot represent empty set in GEN format",
				basic_warning::AMBIGUOUS_OUTPUT);
	}
}

}

#endif /* GEN_FORMATTER_H_ */
