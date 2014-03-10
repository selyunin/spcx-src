/*
 * CIF_automaton_formatter.cpp
 *
 *  Created on: May 11, 2011
 *      Author: goyal
 */

#include "CIF_automaton_formatter.h"
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/transition.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"
#include "utility/shared_ptr_output.h"
#include "core/hybrid_automata/automaton_name_formatter.h"
#include "math/vdom/context_variable_formatter.h"
#include "core/hybrid_automata/automaton_cache.h"

namespace io {

CIF_automaton_formatter::CIF_automaton_formatter(std::ostream& new_os) :
	os(new_os) {
}

CIF_automaton_formatter::~CIF_automaton_formatter() {
}

void CIF_automaton_formatter::file_prologue() {
	os << "";
}

void CIF_automaton_formatter::file_epilogue() {
	os << "" << std::endl;
}

hybrid_automata::hybrid_automaton_visitor::visiting_order CIF_automaton_formatter::get_visiting_order()
{
	return hybrid_automata::hybrid_automaton_visitor::outgoing_transitions;
}

std::string CIF_automaton_formatter::format_ident(const std::string& ident,const hybrid_automata::hybrid_automaton& h) {
	std::string new_str=ident;
	replace(new_str,h.get_name()+".","");
	return new_str;
}

std::string CIF_automaton_formatter::spaceex_to_cif(const std::string& str, CIF_automaton_formatter::predicate_type pred)
{
	std::string res;
	if(pred == CIF_automaton_formatter::inv_predicate)
	{
		res = str;
		replace(res, "==", "=");
		replace(res, "&", ",\n");
	}

	else if (pred == CIF_automaton_formatter::guard_predicate)
	{
		res = str;
		replace(res, "&", "and");
	}

	else if (pred == CIF_automaton_formatter::reset_map_predicate)
	{
		std::string param_ss, val_ss;

		boost::char_separator<char> emp_sep("&");
		boost::tokenizer<boost::char_separator<char> > emp_tokens(str, emp_sep);
		BOOST_FOREACH(std::string emp_t, emp_tokens)
		{
			if(emp_t != *emp_tokens.begin()) {

			param_ss.append(",");
			val_ss.append(",");
			}

			boost::char_separator<char> eql_sep("==");
			boost::tokenizer<boost::char_separator<char> > eql_tokens(emp_t, eql_sep);
			boost::tokenizer<boost::char_separator<char> >::const_iterator it = eql_tokens.begin();
			param_ss.append(*it);
			it++;
			val_ss.append(*it);
		}
		res.append(param_ss + ":=" + val_ss);
	}
	return res;
}

std::string CIF_automaton_formatter::format_with_context(const std::string& ident, const std::string& context) {
	std::string new_str = ident;
	replace(new_str, context+".", "");
	return new_str;
}

void CIF_automaton_formatter::output(hybrid_automata::hybrid_automaton& h)
{
	//output network automaton
	if (hybrid_automata::hybrid_automaton_network* n = dynamic_cast<hybrid_automata::hybrid_automaton_network*>(&h)) {

	std::string context = hybrid_automata::hybrid_automaton_cache::get_automaton(h.get_id())->get_name();

	os << "model " << context << " |[ = " << std::endl;

	variable_id_set net_vars = h.get_variable_ids();
	variable_id_set temp_vars, vars;

	hybrid_automata::label_id_set net_labels = h.get_labels();
	hybrid_automata::label_id_set temp_labels, labels;

	for (variable_id_set::const_iterator it = net_vars.begin(); it != net_vars.end(); ++it) {
			if(it != net_vars.begin())
				os << "; " ;
			os << "var " << variable(*it).get_name() << std::endl;
	}

	for (hybrid_automata::label_id_set::const_iterator it = net_labels.begin(); it != net_labels.end(); ++it) {
		if(!net_vars.empty()) {
			if(hybrid_automata::named_label::is_not_silent(*it))
				os << "; " << "act " << hybrid_automata::named_label::get_name(*it) << std::endl;
		}
		else {
			if(hybrid_automata::named_label::is_not_silent(*it)) {
				if(it != net_labels.begin())
					os << "; ";
			 os << "act " << hybrid_automata::named_label::get_name(*it) << std::endl;
			}
		}
	}
	os << ":: ";
	hybrid_automata::automaton_id_set aut_ids = n->get_automata();

		for (hybrid_automata::automaton_id_set::const_iterator it =
					aut_ids.begin(); it != aut_ids.end(); ++it) {
				temp_vars = vars = net_vars;
				temp_labels = labels = net_labels;
				hybrid_automata::hybrid_automaton::ptr a =
						hybrid_automata::hybrid_automaton_cache::get_automaton(*it);

				set_difference_assign(temp_vars, a->get_variable_ids());
				set_difference_assign(vars, temp_vars);

				set_difference_assign(temp_labels, a->get_labels());
				set_difference_assign(labels, temp_labels);

				if(it != aut_ids.begin())
					os << std::endl << "|| ";
				os << a->get_name() << "(";
				for (variable_id_set::const_iterator it = vars.begin(); it != vars.end(); ++it) {
					if(it != vars.begin())
						os << ", ";
					os << variable(*it).get_name();
				}


				for (hybrid_automata::label_id_set::const_iterator it = labels.begin(); it != labels.end(); ++it) {
					if(!vars.empty()) {
						if(hybrid_automata::named_label::is_not_silent(*it))
								os << ", " << hybrid_automata::named_label::get_name(*it);
					}
					else {
						if(hybrid_automata::named_label::is_not_silent(*it)) {
							if(it != labels.begin())
								os << ", ";
							os << hybrid_automata::named_label::get_name(*it);
							}
					}
				}
				os << ")";
		}

		os << std::endl << "]|" << std::endl;

			//print each base automaton
			for (hybrid_automata::automaton_id_set::const_iterator it =
								aut_ids.begin(); it != aut_ids.end(); ++it) {
				hybrid_automata::hybrid_automaton::ptr a =
										hybrid_automata::hybrid_automaton_cache::get_automaton(*it);
				a->accept(*this);
			}
	} else {
		// output a base automaton
		h.accept(*this);
	}
}

void CIF_automaton_formatter::prologue(hybrid_automata::hybrid_automaton& h)
{
		std::string h_name = h.get_name();
		variable_id_set vars = h.get_variable_ids();
		os << std::endl << "automaton " << h_name + "(";
		for (variable_id_set::const_iterator it = vars.begin(); it != vars.end(); ++it) {
			if(it != vars.begin())
				os << "; ";
			os << "var ";

			os << variable(*it).get_name();
		}

		hybrid_automata::label_id_set labels = h.get_labels();

		for (hybrid_automata::label_id_set::const_iterator it = labels.begin(); it != labels.end(); ++it) {

			if(!vars.empty()) {
				if(hybrid_automata::named_label::is_not_silent(*it))
						os << "; " << "inout act sync " << hybrid_automata::named_label::get_name(*it);
			}
			else {
				if(hybrid_automata::named_label::is_not_silent(*it)) {
					if(it != labels.begin())
						os << "; ";
					os << "inout act sync " << hybrid_automata::named_label::get_name(*it);
				}
			}
		}

		os <<  ") = |(" << std::endl;
}

void CIF_automaton_formatter::epilogue(hybrid_automata::hybrid_automaton& h)
{
	os << std::endl << ")|" << std::endl;
}

void CIF_automaton_formatter::visit(hybrid_automata::hybrid_automaton& h, hybrid_automata::location& l, hybrid_automata::location_id l_id)
{
	context_variable_formatter vform("");
		os << vform;

	if(h.get_location_id(l.get_name()) != 1)
				os << std::endl << ", ";

	os << "mode " << l.get_name() << " = ";

	if (l.get_time_constraints().get_invariant()) {
		std::stringstream invstream;
		invstream << vform << l.get_time_constraints().get_invariant();
		std::string invstring = invstream.str();
		if(invstring.compare("true") != 0) {
			os << "inv " << spaceex_to_cif(invstring,CIF_automaton_formatter::inv_predicate);

			/* In original CIF, inv and flow predicates are not separated with ','. We add comma for simplicity.
			 * The point to be discussed in forthcoming group meeting in Denmark.
			 */
			if(l.get_time_constraints().get_dynamics()) {
			std::stringstream tempstream;
			tempstream << l.get_time_constraints().get_dynamics();
			std::string tempstring = tempstream.str();
			if(tempstring.compare("true") != 0)
				os << " , ";
			}
		}
	}

	if (l.get_time_constraints().get_dynamics()) {
		std::stringstream dynstream;
		dynstream << vform << l.get_time_constraints().get_dynamics();
		std::string dynstring = dynstream.str();
		if(dynstring.compare("true") != 0)
			os << "flow " << spaceex_to_cif(dynstring,CIF_automaton_formatter::inv_predicate);
	}
}

void CIF_automaton_formatter::visit(hybrid_automata::hybrid_automaton& h, hybrid_automata::transition& t, hybrid_automata::transition_id t_id)
	{

	context_variable_formatter vform("");
	os << vform;

	os << "(";
	std::stringstream lstream;
	if (t.get_jump_constraints().get_guard()) {
		std::stringstream gstream;
		gstream << vform << t.get_jump_constraints().get_guard();
		std::string gstring = gstream.str() ;
		if(gstring.compare("true") != 0)
			os << "when " << vform << spaceex_to_cif(gstring,CIF_automaton_formatter::guard_predicate);
	}

	if(t.get_label())
		lstream << vform << hybrid_automata::named_label::get_name(t.get_label());

	if (t.get_jump_constraints().get_transform()) {
		std::stringstream tstream;
		tstream << vform << t.get_jump_constraints().get_transform();
		std::string tstring=tstream.str();
		if(tstring.compare("true") != 0)
			os << " now act " << lstream.str() << " do " << spaceex_to_cif(tstring, CIF_automaton_formatter::reset_map_predicate);
		else
			os << "act " << lstream.str();
	}
	else
		os << " act "<< lstream.str();

	os << ") goto " << h.get_location(t.get_target())->get_name();

	}

}


