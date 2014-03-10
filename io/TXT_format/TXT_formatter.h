/*
 * TXT_formatter.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef TXT_FORMATTER_H_
#define TXT_FORMATTER_H_

#include <stdexcept>
#include <typeinfo>
#include "io/common_output/dispatch_output_formatter.h"

#include "core/discrete/discrete_set.h"
// Sets can't be incomplete types
#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC.h"
#include "core/continuous/polyhedra/polyhedron_upcaster.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"
#include "core/continuous/support_function/spacetime_flowpipe.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "math/vdom/context_variable_formatter.h"
#include "core/hybrid_automata/automaton_name_formatter.h"

/** Forward declarations */
namespace continuous {
template<typename T> class polyhedron;
}

namespace io {

class TXT_formatter;

template<typename T>
class TXT_continuous_set_formatter {
public:
	static void output(output_formatter& of, const T& c) {
		std::string tname = typeid(&c).name();
		throw std::runtime_error(
				"TXT_continuous_set_formatter : can't format object pointed to by "
						+ tname + ".");
	}
	;
};

template< >
class TXT_continuous_set_formatter<continuous::predicate_continuous_set > {
public:
	static void output(output_formatter& of, const continuous::predicate_continuous_set & c) {
		c.print(of.get_os());
	}
	;
};

template<typename scalar_type>
class TXT_continuous_set_formatter<continuous::polyhedron<scalar_type> > {
public:
	static void output(output_formatter& of, const continuous::polyhedron<
			scalar_type> & c) {
		c.print_textual(of.get_os());
	}
	;
};

/** Textual format, readable by sspaceex
 *
 * Tries to upcast to a polyhedron.
 */
class TXT_formatter: public dispatch_output_formatter<
		TXT_continuous_set_formatter, continuous::polyhedron_upcaster> {
public:
	typedef dispatch_output_formatter<TXT_continuous_set_formatter,
			continuous::polyhedron_upcaster> base_class;
	explicit TXT_formatter(std::ostream& os) :
		base_class(os), my_parentheses_around_symbolic_state(false) {
	}
	;
	virtual ~TXT_formatter() {
	}
	;
	virtual void prologue() {
	}
	;
	virtual void epilogue() {
	}
	;
	virtual void output(const discrete::discrete_set& d) {
		hybrid_automata::context_automaton_name_formatter form("",get_context());
		get_os() << form;
		get_os() << d;
	}
	;
	virtual void symbolic_state_element_separator() {
		get_os() << " & ";
	}
	;
	virtual void output(const hybrid_automata::symbolic_state& sstate) {
		if (my_parentheses_around_symbolic_state)
			get_os() << "(";
		base_class::output(sstate);
		if (my_parentheses_around_symbolic_state)
			get_os() << ")";
	}
	;
	virtual void symbolic_state_collection_element_separator() {
		get_os() << "|";
	}
	;
	virtual void output(const continuous::continuous_set& c) {
		context_variable_formatter form(get_context());
		get_os() << form;
		base_class::output(c);
	}
	;
	virtual void output(
			const hybrid_automata::symbolic_state_collection& sstates) {
		// set my_parentheses_around_symbolic_state according to
		// whether there is more than one element in sstates
		my_parentheses_around_symbolic_state = false;
		hybrid_automata::symbolic_state_collection::const_iterator it =
				sstates.begin();
		if (it != sstates.end()) {
			// size of sstates >=1
			++it;
			if (it != sstates.end()) {
				// size of sstates > 1
				my_parentheses_around_symbolic_state = true;
			}
		}
		get_os() << "{";
		base_class::output(sstates);
		get_os() << "}";
		my_parentheses_around_symbolic_state = false;
	}
	;
	static void set_sfm_format(const std::string& opt) {
		my_sfm_format = opt;
	};
	static const std::string& get_sfm_format() {
		return my_sfm_format;
	};
protected:
	bool my_parentheses_around_symbolic_state;
	static std::string my_sfm_format;
};

template<typename scalar_type>
class TXT_continuous_set_formatter<continuous::support_function::sfm_cont_set<scalar_type> > {
public:
	static void output(output_formatter& of, const continuous::support_function::sfm_cont_set<
			scalar_type> & c) {
		TXT_formatter& tf=static_cast<TXT_formatter&>(of);

		if (tf.get_sfm_format().empty() || tf.get_sfm_format() == "code") {
			continuous::support_function::sfm_cont_set<scalar_type>::set_output_format(
					continuous::support_function::sfm_cont_set<scalar_type>::TEXTUAL);
			c.print(of.get_os());
		} else if (tf.get_sfm_format() == "outer") {
			continuous::polyhedron_collection<scalar_type> coll =
					c.get_outer_polytope_collection();
			// use default output format (textual with parenthesis)
			//		continuous::polyhedron_collection<scalar_type>::set_output_format(typename continuous::polyhedron_collection<scalar_type>::output_format());
			//		of.output(coll);

			of.get_os() << "(";
			for (typename continuous::polyhedron_collection<scalar_type>::const_iterator
					it = coll.begin(); it != coll.end(); ++it) {
				if (it != coll.begin())
					of.get_os() << ")|(";
				of.get_os() << *it;
			}
			of.get_os() << ")";

		} else if (tf.get_sfm_format() == "hull") {
			continuous::constr_polyhedron<scalar_type> hull=c.compute_template_hull();
			of.output(hull);
		} else {
			throw basic_exception(
					"Unknown entry for option output-format-sfm:"
							+ tf.get_sfm_format() + ".");
		}
}
	;
};

template<typename scalar_type>
class TXT_continuous_set_formatter<continuous::spacetime_flowpipe<scalar_type> > {
public:
	static void output(output_formatter& of, const continuous::spacetime_flowpipe<
			scalar_type> & c) {
		using namespace continuous;

		TXT_formatter& tf=static_cast<TXT_formatter&>(of);

		if (tf.get_sfm_format().empty() || tf.get_sfm_format() == "code") {
			c.print(of.get_os());
		} else if (tf.get_sfm_format() == "outer") {
			typename spacetime_flowpipe<scalar_type>::cut_point_method m; // use default
			polyhedron_collection<scalar_type> coll =
					c.compute_outer_polyhedra(m);

			of.get_os() << "(";
			for (typename polyhedron_collection<scalar_type>::const_iterator
					it = coll.begin(); it != coll.end(); ++it) {
				if (it != coll.begin())
					of.get_os() << ")|(";
				of.get_os() << *it;
			}
			of.get_os() << ")";

		} else if (tf.get_sfm_format() == "hull") {
			typename spacetime_flowpipe<scalar_type>::cut_point_method m(
					spacetime_flowpipe<scalar_type>::cut_point_method::FIXED_COUNT_SAME_SIZE_PIECES,
					1); // use template hull
			polyhedron_collection<scalar_type> coll =
					c.compute_outer_polyhedra(m);

			if (coll.begin()!=coll.end()) {
				of.output(**coll.begin());
			} else {
				constr_polyhedron<scalar_type> poly = constr_polyhedron<scalar_type>::empty_poly();
				of.output(poly);
			}
		} else {
			throw basic_exception(
					"Unknown entry for option output-format-sfm:"
							+ tf.get_sfm_format() + ".");
		}
}
	;
};




}

#endif /* TXT_FORMATTER_H_ */
