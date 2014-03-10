/*
 * SX_automaton_formatter.cpp
 *
 *  Created on: Sep 16, 2010
 *      Author: frehse
 */

#include "SX_automaton_formatter.h"

#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/transition.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"
#include "utility/shared_ptr_output.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/hybrid_automata/automaton_name_formatter.h"
#include "core/hybrid_automata/hybrid_automaton_utility.h"
#include "math/vdom/context_variable_formatter.h"

namespace io {

SX_automaton_formatter::SX_automaton_formatter(std::ostream& new_os) :
	os(new_os) {
}

SX_automaton_formatter::~SX_automaton_formatter() {
}

void SX_automaton_formatter::file_prologue() {
	os << "<?xml version=\"1.0\" encoding=\"iso-8859-1\" ?>" << std::endl;
	os << "<sspaceex xmlns=\"http://www-verimag.imag.fr/xml-namespaces/sspaceex\" version=\"0.2\" math=\"SpaceEx\">" << std::endl;
}

void SX_automaton_formatter::file_epilogue() {
	os << "</sspaceex>" << std::endl;
}

std::string SX_automaton_formatter::format_ident(const std::string& ident,const hybrid_automata::hybrid_automaton& h) {
	std::string new_str=ident;
	replace(new_str,h.get_name()+".","");
	return new_str;
}

void SX_automaton_formatter::prologue(hybrid_automata::hybrid_automaton& h) {
	os << "<component id=\"" << h.get_name() << "\">" << std::endl;
	// set formatting for variables and locations
	context_variable_formatter vform(h.get_name());
	os << vform;
	hybrid_automata::context_automaton_name_formatter aform("",h.get_name());
	os << aform;

	variable_id_set vars = h.get_variable_ids();
	variable_id_set contrvars = h.get_variable_ids();
	variable_id_set inpvars = h.get_input_variables();
	variable_id_set constvars = h.get_const_variables();
	set_difference_assign(contrvars, inpvars);
	for (variable_id_set::const_iterator it = vars.begin(); it != vars.end(); ++it) {
		os << "<param name=\"" << variable(*it) << "\" type=\"real\" d1=\"1\" d2=\"1\" local=\"false\" dynamics=\"";
		if (constvars.find(*it)!=constvars.end())
			os << "const";
		else
			os << "any";
		os << "\" controlled=\"";
		if (contrvars.find(*it)!=contrvars.end())
			os << "true";
		else
			os << "false";
		os << "\"/>" << std::endl;
	}

	// output labels
	hybrid_automata::label_id_set labels = h.get_labels();
	for (hybrid_automata::label_id_set::const_iterator it = labels.begin(); it != labels.end(); ++it) {
		os << "<param name=\"" << format_ident(hybrid_automata::named_label::get_name(*it),h) << "\" type=\"label\" local=\"false\"/>" << std::endl;
	}
}

void SX_automaton_formatter::epilogue(hybrid_automata::hybrid_automaton& h) {
	os << "</component>" << std::endl;
}

void SX_automaton_formatter::visit(hybrid_automata::hybrid_automaton& h, hybrid_automata::location& l, hybrid_automata::location_id l_id) {
	using namespace math;
	using namespace hybrid_automata;
	typedef global_types::float_type float_type;

	// set formatting for variables and locations
	std::string hname = h.get_name();

	context_variable_formatter vform(hname);
	os << vform;
	// need to make a new formatter because the one from the prologue
	// gets destroyed at the end of the prologue
	hybrid_automata::context_automaton_name_formatter lform("",
			hname);
	os << lform;

	os << "<location id=\"" << l_id << "\"";
	os << " name=\"" << l.get_name() << "\">" << std::endl;

	os << "   <note>";
	// @todo: Show canonlicalized location constraints
	hybrid_automaton::const_ptr p = hybrid_automaton_cache::get_automaton(h.get_name());
	if (p) {
		os << "Location constraints: ";
		location_constraint_set lcs = canonicalize_location_constraint(p, l_id);
		os << lcs << std::endl << std::endl;
	}

	const hybrid_automata::time_constraints& tcons = l.get_time_constraints();
	if (tcons.get_dynamics()) {
		// if dynamics in matrix form, show eigenvalues and vectors.
		typedef continuous::ode_affine_dynamics<float_type> ode_aff_dyn;
		ode_aff_dyn::const_ptr dp = boost::dynamic_pointer_cast<
				const ode_aff_dyn>(tcons.get_dynamics());

		if (dp && dp->get_A().size1()>0 && dp->get_A().size2()>0) {
			os << "State variables: " << dp->get_A().codomain() << std::endl;

			using namespace math;
			vector<float_type> v_real_d, v_imag_d;
			matrix<float_type> V_d;
			compute_eigenvalues(dp->get_A().get_matrix(), v_real_d, v_imag_d, V_d);
			os << "real parts of Eigenvalues: "
					<< v_real_d << ", imag parts of Eigenvalues: " << v_imag_d << std::endl;
			os << "Eigenvectors (columns): " << V_d << std::endl;
		}
	}
	os << " </note>" << std::endl;
	if (tcons.get_invariant()) {
		os << "   <invariant>";
		std::stringstream invstream;
		invstream << vform << l.get_time_constraints().get_invariant();
		std::string invstring = invstream.str();
		os << string_to_xml(invstring) << " </invariant>" << std::endl;
	}
	if (tcons.get_dynamics()) {
		os << "   <flow>";
		std::stringstream dynstream;
		dynstream << vform << l.get_time_constraints().get_dynamics();
		std::string dynstring = dynstream.str();
		os << string_to_xml(dynstring) << " </flow>" << std::endl;
	}
	os << "</location>" << std::endl;
}

void SX_automaton_formatter::visit(hybrid_automata::hybrid_automaton& h, hybrid_automata::transition& t, hybrid_automata::transition_id t_id) {
	// set formatting for variables and locations
	context_variable_formatter vform(h.get_name());
	os << vform;

	os << "<transition ";
	os << "source=\"" << t.get_source() << "\" ";
	os << "target=\"" << t.get_target() << "\">" << std::endl;
	if (t.get_label()) {
		os << "   <label>" << format_ident(hybrid_automata::named_label::get_name(t.get_label()),h) << "</label>"
			<< std::endl;
	}
	if (t.get_jump_constraints().get_guard()) {
		os << "   <guard> ";
		std::stringstream gstream;
		gstream << vform << t.get_jump_constraints().get_guard();
		std::string gstring = gstream.str();
		os << string_to_xml(gstring) << " </guard>" << std::endl;
	}
	if (t.get_jump_constraints().get_transform()) {
		os << "   <assignment> ";
		std::stringstream tstream;
		tstream << vform << t.get_jump_constraints().get_transform();
		std::string tstring=tstream.str();
		os << string_to_xml(tstring) << " </assignment>" << std::endl;
	}
	os << "</transition>" << std::endl;
}

SX_automaton_formatter::network_visits SX_automaton_formatter::get_network_visits() {
	return only_composition;
}

}


