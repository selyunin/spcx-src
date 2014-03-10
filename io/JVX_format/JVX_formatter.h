/*
 * JVX_formatter.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef JVX_FORMATTER_H_
#define JVX_FORMATTER_H_

#include <stdexcept>
#include <typeinfo>
#include "io/common_output/dispatch_output_formatter.h"

//#include "core/discrete/discrete_set.h"
// Sets can't be incomplete types
#include "core/continuous/support_function_provider_utility.h"
#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/continuous/polyhedra/polyhedron_upcaster.h"
#include "core/continuous/support_function_provider_upcaster.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"
#include "core/continuous/support_function/sf_base/sf_set.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "math/vdom/trajectory.h"
#include "core/continuous/support_function/spacetime_flowpipe.h"

#include "math/ode_solving/traj_simu/continuous_set_simulation.h"

/** Forward declarations */
namespace continuous {
template<typename T> class polyhedron;
}

namespace io {

template<typename T>
class JVX_continuous_set_formatter {
public:
	static void output(output_formatter& of, const T& c) {
		std::string tname = typeid(&c).name();
		throw std::runtime_error(
				"JVX_continuous_set_formatter : can't format object pointed to by "
						+ tname + ".");
	}
	;
};

template<>
class JVX_continuous_set_formatter<continuous::support_function_provider > {
public:
	static void output(output_formatter& of,
			const continuous::support_function_provider& c) {
		positional_vdomain dom(c.get_variable_ids());
		if (!of.get_output_variables().empty()) {
			dom = positional_vdomain(of.get_output_variables());
		}
		std::list<math::vector<global_types::float_type> > dir_list;

		continuous::support_function::choose_directions(dom, dir_list);
		typedef std::set<math::vdom_vector<global_types::float_type>,
				math::numeric::lex_comp_less<global_types::float_type,
						math::vdom_vector> > dir_set_type;

		dir_set_type dir_set;

		for (std::list<math::vector<global_types::float_type> >::const_iterator
				it = dir_list.begin(); it != dir_list.end(); ++it) {
			math::vdom_vector<global_types::float_type> d(dom, *it);
			dir_set.insert(d);
		}

		of.output(compute_outer_poly(c, dir_set));
	}
	;
};

template<typename scalar_type>
class JVX_continuous_set_formatter<continuous::polyhedron<scalar_type> > {
public:
	static void output(output_formatter& of, const continuous::polyhedron<
			scalar_type> & c) {
		if (of.get_output_variables().empty()) {
			c.print_JVX(of.get_os());
		} else {
			positional_vdomain dom(of.get_output_variables());
			// check if we need to perform projection
			if (set_contains(list_to_set(of.get_output_variables()),c.get_variable_ids())) {
				// no projection necessary
				c.print_JVX(of.get_os(),of.get_output_variables());
			} else {
				// output as support_function
				//continuous::support_function_provider::ptr sup=continuous::support_function_provider::ptr(c.clone());
				//continuous::support_function::sf_unary<scalar_type> s(sup);
				//of.output(*sup);
				JVX_continuous_set_formatter<continuous::support_function_provider>::output(of,c);
			}
		}
	}
	;
};

template<typename scalar_type>
class JVX_continuous_set_formatter<continuous::support_function::sfm_cont_set<scalar_type> > {
public:
	static void output(output_formatter& of, const continuous::support_function::sfm_cont_set<
			scalar_type> & c) {
		for (size_t i = 0; i < c.get_size(); ++i) {
			continuous::constr_polyhedron<scalar_type> p = c.get_polytope(i);
			of.output(p);
			//of.get_os() << std::endl << std::endl;
		}
	}
	;
};

/** Textual format, readable by sspaceex
 *
 * Tries to upcast to a polyhedron.
 */
template<typename T> class JVX_formatter_caster: public dispatching::caster_list<
		continuous::support_function_provider_upcaster<T>,
		dispatching::caster_list<continuous::polyhedron_upcaster<T>,
				dispatching::convertible_caster<T,
						continuous::support_function::sfm_cont_set<
								global_types::float_type> > > > {
};

class JVX_formatter: public dispatch_output_formatter<
		JVX_continuous_set_formatter, JVX_formatter_caster> {
public:
	typedef dispatch_output_formatter<JVX_continuous_set_formatter,
			JVX_formatter_caster> base_class;
	explicit JVX_formatter(std::ostream& os) :
		base_class(os) {
	}
	;
	virtual ~JVX_formatter() {
	}
	;
	virtual void prologue() {
		get_os()
				<< "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>"
				<< std::endl;
		get_os()
				<< "<!DOCTYPE jvx-model SYSTEM \"http://www.javaview.de/rsrc/jvx.dtd\"> "
				<< std::endl;
		get_os() << "<jvx-model>" << std::endl;
	}
	;
	virtual void epilogue() {
		get_os() << "</jvx-model>" << std::endl;
	}
	;
	using base_class::output;
	virtual void output(const discrete::discrete_set& d) {
	}
	;
	virtual void symbolic_state_element_separator() {
	}
	;
	virtual void output(const hybrid_automata::symbolic_state& sstate) {
		base_class::output(sstate);
	}
	;
	virtual void symbolic_state_collection_element_separator() {
		get_os() << "\n";
	}
	;
	virtual void output(const hybrid_automata::symbolic_state_collection& sstates) {
		get_os() << "<geometries>" << std::endl;
		base_class::output(sstates);
		get_os() << "</geometries>" << std::endl;
	}
	;
	void output(
			const continuous::polyhedron_collection<global_types::float_type>& polys) {
		for (continuous::polyhedron_collection<global_types::float_type>::const_iterator it =
				polys.begin(); it != polys.end(); ++it) {
			output(**it);
		}
	};
	/** Output a trajectory in the state space. */
	void output(const math::trajectory<global_types::float_type>& traj) {
		variable_id_list vars;
		if (get_output_variables().empty()) {
			const variable_id_set& vis=traj.get_variable_ids();
			for (variable_id_set::const_iterator it=vis.begin();it!=vis.end();++it) {
				vars.push_back(*it);
			}
		} else {
			vars=get_output_variables();
		}
		// Construct a vector of indices
		std::vector<unsigned int> indices(vars.size());
		unsigned int i=0;
		for (variable_id_list::const_iterator it=vars.begin();it!=vars.end();++it) {
			indices[i]=traj.domain().pos(variable(*it));
			++i;
		}

		get_os()  << "<geometry name=\"trajectory\">" << std::endl;
		get_os()  << "<pointSet dim=\"";
		get_os()  << vars.size();
		get_os()  << "\" point=\"hide\">" << std::endl;
		get_os()  << "   <points> " << std::endl;

		for (i=0;i<traj.size();++i) {
			get_os() << "      <p> ";
			for (unsigned int j=0;j<indices.size();++j) {
				get_os() << traj.get_states()(i,indices[j]) << " ";
			}
			get_os()  << " </p>";
			get_os() << std::endl;
		}
		get_os() << std::endl << "   </points>" << std::endl;
		get_os() << "</pointSet>" << std::endl;

		get_os() << "<lineSet line=\"show\">" << std::endl;
		get_os() << "   <lines>" << std::endl;
		for (i=0;i+1<traj.size();++i) {
				get_os() << "      <l>";
				get_os() << i << " " << i+1;
				get_os() << " </l>" << std::endl;
		}
		get_os() << std::endl << "   </lines>" << std::endl;
		get_os() << "</lineSet>" << std::endl;

		get_os() << "</geometry>" << std::endl;
	}
	;
	/** Output a continuous_set_simulation. (outputs the trajectory part, and not the root part) */
	void output(const continuous::continuous_set_simulation<global_types::float_type>& cset_simu) {
		for (continuous::continuous_set_simulation<global_types::float_type>::const_iterator
				it = cset_simu.begin(); it != cset_simu.end(); ++it) {
			output(it->second);
		}
		get_os() << std::endl;
	}
	;
	void output(
			const continuous::spacetime_flowpipe<global_types::float_type>& flowp) {
		typedef continuous::spacetime_flowpipe<global_types::float_type> flowpipe;
		flowpipe::cut_point_method m;
		double my_eps = 0.001;
		if (my_eps > 1e-6) {
			m.type = flowpipe::cut_point_method::MIN_CONCAVE_PIECES;
			m.approx_error = flowpipe::error_type(0, my_eps);
			m.lower_is_error_reference = true;
		} else {
			m.type = flowpipe::cut_point_method::ALL_PIECES;
		}
		continuous::polyhedron_collection<global_types::float_type> polys =
				flowp.compute_outer_polyhedra(m);

		double old_my_eps = my_eps;
		// set tolerance to zero before outputting, otherwise error will accumulate
		my_eps = 0;
		output(polys);
		my_eps = old_my_eps;
	}
	;
	void output(const continuous::continuous_set& c) {
		// dynamic cast to continuous set classes that
		// are not in the typelist
		if (const continuous::continuous_set_simulation<global_types::float_type> * p =
				dynamic_cast<const continuous::continuous_set_simulation<global_types::float_type> *>(&c)) {
			output(*p);
		} else if (const continuous::polyhedron_collection<global_types::float_type> * p =
				dynamic_cast<const continuous::polyhedron_collection<global_types::float_type> *>(&c)) {
			output(*p);
		} else if (const continuous::spacetime_flowpipe<global_types::float_type> * p =
				dynamic_cast<const continuous::spacetime_flowpipe<global_types::float_type> *>(&c)) {
			output(*p);
		} else{
			// otherwise call base class output to try double dispatch
			base_class::output(c);
		}
	}
	;
};


class JVX_OFFLINE_formatter: public JVX_formatter {
	typedef JVX_formatter base_class;
public:
	explicit JVX_OFFLINE_formatter(std::ostream& os) :
			base_class(os) {
		}
		;
		virtual ~JVX_OFFLINE_formatter() {
		}
		;

	virtual void prologue() {
			get_os()
					<< "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>"
					<< std::endl;
			get_os() << "<jvx-model>" << std::endl;
		}
		;
};

}

#endif /* JVX_FORMATTER_H_ */
