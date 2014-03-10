/*
 * CIF_parser.cpp
 *
 *  Created on: Jul 7, 2009
 *      Author: gvincent
 */
#include <sstream> // needed to convert string into an input stream
#include <vector> // used to store the words of the string
#include <iterator>
#include <iostream>
#include <fstream>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>

#include "io/CIF_format/CIF_parser.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/predicates/valuation_function_tree_utility.h"
#include "io/common_input/scalar_node_creation.h"
#include "core/continuous/predicate_continuous_set_constructors.h"
#include "core/discrete/singleton_set.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/symbolic_states/symbolic_state_collection_stl_list.h"
#include "core/discrete/discrete_set_operators.h"
#include "io/common_input/symbol_table.h"
#include "io/common_input/symbol_table_operators.h"

using namespace parser;
using namespace hybrid_automata;

namespace CIF_predicate_parser {

dyn_ast& dyn_ast::operator+=(std::string const& inv_or_flow) {
	if (inv_or_flow.compare("inv") == 0)
		my_type = "inv";
	else if (inv_or_flow.compare("flow") == 0)
		my_type = "flow";
	return *this;
}

dyn_ast& dyn_ast::operator-=(predicate_ast const& pred) {
	my_pred = pred;
	return *this;
}
/* += operator to assign guard predicate to the edge */
edge_ast& edge_ast::operator+=(predicate_ast const& pred) {
	my_guard = pred;
	return *this;
}

/* -= operator to assign update predicate to the edge */
edge_ast& edge_ast::operator-=(predicate_ast const& pred) {
	my_update = pred;
	return *this;
}

/* *= operator to assign target location to the edge */
edge_ast& edge_ast::operator*=(std::string const& target) {
	my_target = target;
	return *this;
}

/* >= operator to assign label to the edge */
edge_ast& edge_ast::operator>=(std::string const& label) {
	my_label = label;
	return *this;
}

/* += operator to assign modes vector to the automaton*/
automaton_ast& automaton_ast::operator+=(std::vector<mode_ast> const& modes) {
	my_modes = modes;
	return *this;
}

/* -= operator to assign name to the automaton */
automaton_ast& automaton_ast::operator-=(std::string const& name) {
	my_name = name;
	return *this;
}

/* >= operator to assign local declarations to the automaton */
automaton_ast& automaton_ast::operator>=(
		std::vector<declaration_ast>& declaration) {
	for (std::vector<declaration_ast>::size_type i = 0; i < declaration.size(); ++i)
		declaration[i].is_local = true;

	if (my_declaration.empty())
		my_declaration = declaration;
	else
		my_declaration.insert(my_declaration.end(), declaration.begin(),
				declaration.end());

	return *this;
}

/* *= operator to assign non-local declarations to the automaton */
automaton_ast& automaton_ast::operator*=(
		std::vector<declaration_ast>& declaration) {
	for (std::vector<declaration_ast>::size_type i = 0; i < declaration.size(); ++i)
		declaration[i].is_local = false;

	if (my_declaration.empty())
		my_declaration = declaration;
	else
		my_declaration.insert(my_declaration.end(), declaration.begin(),
				declaration.end());
	return *this;
}

/* += operator to assign instantiation name to the automaton_dec */
automaton_dec_ast& automaton_dec_ast::operator+=(
		std::string const& instance_name) {
	my_def_name = my_instance_name = instance_name;
	return *this;
}

/* -= operator to assign automaton definition name to the automaton_dec */
automaton_dec_ast automaton_dec_ast::operator-=(std::string const& def_name) {
	my_def_name = def_name;
	return *this;
}

/* *= operator to assign parameters vector to the automaton_dec */
automaton_dec_ast automaton_dec_ast::operator*=(
		std::vector<expr_ast> const& parameters) {
	my_parameters = parameters;
	return *this;
}

/* += operator to assign automata vector to the network */
network_ast& network_ast::operator+=(std::vector<automaton_ast> const& automata) {
	my_automata = automata;
	return *this;
}

/* -= operator to assign name to the network */
network_ast& network_ast::operator-=(std::string const& name) {
	my_name = name;
	return *this;
}

/* *= operator to assign declarations vector to the network */
network_ast& network_ast::operator*=(
		std::vector<declaration_ast> const& declaration) {
	my_declaration = declaration;
	return *this;
}

/* *= operator to assign automata declarations to the network */
network_ast& network_ast::operator>=(
		std::vector<automaton_dec_ast> const& automata_dec) {
	my_automata_dec = automata_dec;
	return *this;
}

/* += operator to assign parameter type to the declaration*/
declaration_ast& declaration_ast::operator+=(std::string const& type) {
	my_dec.append(type);
	my_dec.append(" ");
	return *this;
}

/* -= operator to assign comma separated expressions to the declaration*/
declaration_ast& declaration_ast::operator-=(
		std::vector<expr_ast> const& expression) {
	std::vector<expr_ast>::size_type i = 0;
	for (; i < expression.size() - 1; i++) {
		my_dec.append(expression[i].my_exp);
		my_dec.append(" , ");
	}
	my_dec.append(expression[i].my_exp);
	return *this;
}

/* += operator to assign name of the parameter to the expression*/
expr_ast& expr_ast::operator+=(std::string const& param_name) {
	my_exp.append(param_name);
	return *this;
}

class generic_visitor: public boost::static_visitor<void> {
	std::ostringstream& stm_;
public:
	generic_visitor(std::ostringstream& stm) :
		stm_(stm) {
	}

	template<typename T> void operator()(T& t) const {
		stm_ << t;
	}
};

/* -= operator to assign value of the parameter to the expression*/
expr_ast& expr_ast::operator-=(boost::variant<int, double>& param_val) {
	std::ostringstream stm;
	my_exp.append(" = ");
	boost::apply_visitor(generic_visitor(stm), param_val);
	my_exp.append(stm.str());
	return *this;
}

/* += operator either to assign name or to set the initial attribute true*/
mode_ast& mode_ast::operator+=(std::string const& name) {
	if (name.compare("initial") == 0)
		initial = true;
	else
		my_name = name;
	return *this;
}

/* -= operator to assign dynamics to the mode*/
mode_ast& mode_ast::operator-=(std::vector<dyn_ast> const& dyns) {
	dyn_ast dyn;
	for (std::vector<dyn_ast>::size_type i = 0; i < dyns.size(); ++i) {

		if (dyns[i].my_type.compare("inv") == 0)
			dyn.my_type = "inv";
		else if (dyns[i].my_type.compare("flow") == 0)
			dyn.my_type = "flow";

		dyn.my_pred = dyns[i].my_pred;

		my_dyns.push_back(dyn);
	}
	return *this;
}

/* *= operator to assign transition to the mode*/
mode_ast& mode_ast::operator*=(std::vector<edge_ast> const& edges) {
	my_edges = edges;
	return *this;
}

CIF_definition::CIF_definition() {
	using qi::_val;
	using qi::double_;
	using qi::int_;
	using qi::_1;
	using qi::lit;
	using qi::_a;
	using qi::alnum;
	using ascii::string;
	using boost::phoenix::push_back;

	/**
	 * - Parsing declarations in the network -
	 * Only Real/Integer values are considered
	 */
	value = double_ | int_;

	inv_flow_token %= lit("inv") | lit("flow");

	static_type %= lit("real") | lit("nat") | lit("int");

	param_type %= lit("cont") | lit("var") | lit("disc") | lit("clock") | lit(
			"const") | lit("act");

	controlled %= lit("control");

	initial %= lit("initial");

	//expression in a declaration is of type: A or A = 2
	dec_expression = token[_val += _1] >> -(lit('=') >> value[_val -= _1]);

	// expressions are separated by ','

	dec_expressions = dec_expression[push_back(_val, _1)] % ',';

	declaration = (param_type[_val += _1] >> -controlled[_val += _1]
			>> static_type[_val += _1] >> dec_expressions[_val -= _1]) // local parameter declaration
			| (-lit("inout") >> param_type[_val += _1] >> -static_type[_val
					+= _1] >> -lit("sync") >> dec_expressions[_val -= _1]) // formal parameter declaration
	;

	// declarations are separated by ';'

	declarations = -declaration[push_back(_val, _1)] % ';';

	update %= lit("do") >> or_equa;

	guard %= lit("when") >> or_equa;

	edge = -lit("(") >> -guard[_val += _1] >> -lit("now") >> -(lit("act")
			>> token[_val >= _1]) >> -update[_val -= _1] >> -lit(")") >> lit(
			"goto") >> token[_val *= _1];

	//dyn %= (string("inv") >> or_equa) | (string("flow") >> or_equa) ;

	/**
	 * An empty string - string("") - is only to keep dyn structure intact.
	 * dynamics type is being taken care of later depending on the expression
	 * variable whether it is primed (flow) or not.
	 */
	dyn = -(inv_flow_token[_val += _1]) >> or_equa[_val -= _1];

	/**
	 * For the declarations of type: inv/flow V' = Qi - Qo,
	 * 								 , Qi = n * 5.0
	 */
	dyns = *(dyn[push_back(_val, _1)] % ',');

	/**
	 * ****** If required later ******
	 * For the declarations where flow and invariants are declared separately
	 * i.e., the declarations of type: inv V' = Qi - Qo,
	 * 								 , Qi = n * 5.0
	 * 								 , flow A' = V
	 * 								 , Qo = 2 * n
	 */
	//dyns = *(((string("inv") | string("flow")) >> (dyn[push_back(_val, _1)] % ',')) % ',');


	edges = *(edge[push_back(_val, _1)]);

	mode = token[_val += _1] >> lit("=") >> -initial[_val += _1] >> dyns[_val
			-= _1] >> edges[_val *= _1];

	modes = mode[push_back(_val, _1)] % ',';

	automaton = (token[_val -= _1] >> lit(":") >> lit("|(")
			>> -declarations[_val >= _1] >> lit("mode") >> modes[_val += _1]
			>> lit(")|")) | (lit("automaton") >> token[_val -= _1] >> lit("(")
			>> -declarations[_val *= _1] >> lit(")") >> lit("=") >> lit("|(")
			>> -declarations[_val >= _1] >> lit("mode") >> modes[_val += _1]
			>> lit(")|"));

	automaton_dec = token[_val += _1] >> -(lit(":") >> token[_val -= _1])
			>> lit("(") >> -dec_expressions[_val *= _1] >> lit(")");

	// automata are separated by "||"
	automata = automaton[push_back(_val, _1)] % lit("||");

	automata_defs = *(automaton[push_back(_val, _1)]);

	// automata declarations are seperated by "||"
	automata_decs = automaton_dec[push_back(_val, _1)] % lit("||");

	network = (lit("model") >> token[_val -= _1] >> lit("(") >> lit(")")
			>> lit("=") >> lit("|[") >> declarations[_val *= _1] >> lit("::")
			>> automata[_val += _1] >> lit("]|")) | (lit("model") >> token[_val
			-= _1] >> lit("(") >> lit(")") >> lit("=") >> lit("|[")
			>> -declarations[_val *= _1] >> lit("::") >> automata_decs[_val
			>= _1] >> lit("]|") >> automata_defs[_val += _1]) | automata[_val
			+= _1] //if there is only one automaton.
	;

}

/**
 * Return formal parameters as keys in a function definition.
 * The mapping between the key and its value is positional.
 */
std::vector<std::string> get_keys(std::vector<declaration_ast> declarations) {
	std::vector<std::string> keys;
	for (std::vector<declaration_ast>::size_type i = 0; i < declarations.size(); i++) {
		if (!declarations[i].is_local) {
			std::string buf; // Have a buffer string
			std::vector<std::string> tokens; // Create vector to hold our tokens

			std::stringstream ss(declarations[i].my_dec); // Insert the string into a stream


			while (ss >> buf)
				tokens.push_back(buf);

			for (std::vector<std::string>::size_type j = 0; j < tokens.size(); j++) {
				if ((tokens[j].compare("real") == 0) || (tokens[j].compare(
						"int") == 0) || (tokens[j].compare("cont") == 0)
						|| (tokens[j].compare("disc") == 0)
						|| (tokens[j].compare("var") == 0)
						|| (tokens[j].compare("clock") == 0)
						|| (tokens[j].compare("act") == 0)
						|| (tokens[j].compare(",") == 0)) {// do nothing
				} else
					keys.push_back(tokens[j]);
			}
		}
	}
	return keys;
}

/** (1)  Create a vector of symbols from the declaration. A declaration is of type: cont control real A = 2.0, B;
 * where the symbols' declarations are separated by ','. The symbols are later stuffed into symbol table.
 * (2) 	 Set the initial predicate - the continuous part of an initial symbolic state collection for the network.
 * (3) 	 Push the clock symbols into the 'clocks' vector.
 **/
std::vector<parser::symbol> get_new_symbols(declaration_ast const& declaration,
		tree::node::ptr& init_pred, std::vector<std::string>& clocks,
		parse_policy const& ppol) {

	tree::node::ptr pred_left, pred_right, pred_op;

	std::string s_name;
	boost::any s_value;
	bool is_clock = false;
	bool local = declaration.is_local;
	symbol::symbol_type s_type;
	symbol::data_type d_type = symbol::REAL;
	symbol::dynamics_type dyn_type = symbol::ANY;
	bool controlled = false;
	std::vector<parser::symbol> s_vector;
	std::vector<std::string>::size_type i;

	std::string buf; // Have a buffer string
	std::vector<std::string> tokens; // Create vector to hold our tokens

	std::stringstream ss(declaration.my_dec); // Insert the string into a stream

	while (ss >> buf)
		tokens.push_back(buf);

	for (i = 0; i < tokens.size(); i++) {
		std::string token = tokens[i];

		if (token.compare("const") == 0) {
			dyn_type = symbol::CONSTANT;
			s_type = symbol::CONST_VALUE;
		} else if ((token.compare("var") == 0) || (token.compare("disc") == 0)
				|| (token.compare("cont") == 0))
			s_type = symbol::VARIABLE;

		else if (token.compare("clock") == 0) {
			s_type = symbol::VARIABLE;
			is_clock = true;
		}

		else if (token.compare("control") == 0)
			controlled = true;

		else if ((token.compare("int") == 0) || (token.compare("nat") == 0))
			d_type = symbol::INT;

		else if ((token.compare("real") == 0))
			d_type = symbol::REAL;

		else if ((token.compare("act") == 0))
			s_type = symbol::LABEL;

		else
			break;
	}

	for (; i < tokens.size(); i++) {
		if (tokens[i].compare(",") != 0) {
			s_name = tokens[i];
			if (i + 1 < tokens.size() && (tokens[i + 1].compare("=") == 0)) {
				s_value = tokens[i + 2];

				pred_left = tree::node::ptr(
						valuation_functions::variable_node_creator::create(
								tokens[i], "", 0, ppol));
				pred_right = predicate_parser::create_scalar_const_node(
						tokens[i + 2]);
				pred_op = tree::node_ptr(
						new valuation_functions::comparison_node(EQ, pred_left,
								pred_right));
// the following is now done in add_to_init_pred
//				if (init_pred)
//					init_pred = tree::node_ptr(
//							new valuation_functions::boolean_node(AND,
//									init_pred, pred_op));
//				else
//					init_pred = pred_op;

				i = i + 2;
			}

			s_vector.push_back(
					symbol(s_name, s_type, d_type, dyn_type, s_value, local,
							controlled));

			// if clock, push it into the vector
			if (is_clock)
				clocks.push_back(s_name);
		}
	}
	return s_vector;
}

/** This is a copy of get_new_symbols, which generates the predicate using the symbol mapping from the
 *  symbol table
 **/
void add_to_init_pred(declaration_ast const& declaration, tree::node::ptr& init_pred,
		const parser::symbol_table& s_table, parse_policy const& ppol) {
	tree::node::ptr pred_left, pred_right, pred_op;

	std::string s_name;
	boost::any s_value;
	bool is_clock = false;
	bool local = declaration.is_local;
	symbol::symbol_type s_type;
	symbol::data_type d_type = symbol::REAL;
	symbol::dynamics_type dyn_type = symbol::ANY;
	bool controlled = false;
	std::vector<parser::symbol> s_vector;
	std::vector<std::string>::size_type i;

	std::string buf; // Have a buffer string
	std::vector<std::string> tokens; // Create vector to hold our tokens

	std::stringstream ss(declaration.my_dec); // Insert the string into a stream

	while (ss >> buf)
		tokens.push_back(buf);

	for (i = 0; i < tokens.size(); i++) {
		std::string token = tokens[i];

		if (token.compare("const") == 0) {
			dyn_type = symbol::CONSTANT;
			s_type = symbol::CONST_VALUE;
		} else if ((token.compare("var") == 0) || (token.compare("disc") == 0)
				|| (token.compare("cont") == 0))
			s_type = symbol::VARIABLE;

		else if (token.compare("clock") == 0) {
			s_type = symbol::VARIABLE;
			is_clock = true;
		}

		else if (token.compare("control") == 0)
			controlled = true;

		else if ((token.compare("int") == 0) || (token.compare("nat") == 0))
			d_type = symbol::INT;

		else if ((token.compare("real") == 0))
			d_type = symbol::REAL;

		else if ((token.compare("act") == 0))
			s_type = symbol::LABEL;

		else
			break;
	}

	for (; i < tokens.size(); i++) {
		if (tokens[i].compare(",") != 0) { // if not all characters in tokens[i] are commas
			s_name = tokens[i];
			if (i + 1 < tokens.size() && (tokens[i + 1].compare("=") == 0)) {
				s_value = tokens[i + 2];

				// don't add const_value symbols to initial predicate
				if (s_table.get_symbol(tokens[i]).my_symbol_type
						!= symbol::CONST_VALUE) {
					//				pred_left = tree::node::ptr(
					//						valuation_functions::variable_node_creator::create(
					//								tokens[i], "", 0, ppol));
					variable_id v_id = variable::get_variable_id(
							s_table.get_symbol(tokens[i]).my_name);
					tree::node::ptr pred_left = tree::node::ptr(
							new valuation_functions::variable_node(v_id));
					pred_right = predicate_parser::create_scalar_const_node(
							tokens[i + 2]);
					pred_op = tree::node_ptr(
							new valuation_functions::comparison_node(EQ,
									pred_left, pred_right));
					init_pred = valuation_functions::boolean_and(init_pred, pred_op);
				}

				i = i + 2;
			}


		}
	}
}

/** Creates a universe set of discrete states
 *
 * The universe set consists of all locations of all automata. */
discrete::discrete_set::ptr create_universe_discrete_states() {
	discrete::discrete_set::ptr temp_dset;
	// universe location constraint
	hybrid_automata::location_constraint_set temp_lcset;
	temp_dset = discrete::discrete_set::ptr(
		new discrete::singleton_set(temp_lcset));
	return temp_dset;
}

/** Creates a set of discrete states from the modes of an automaton
 *
 * If there are no modes, the universe location constraint is returned.
 * */
discrete::discrete_set::ptr create_discrete_states(
		const hybrid_automata::hybrid_automaton::ptr& new_automaton,
		std::vector<std::string>& initial_modes) {
	/** Create a location constraint only for initial mode and
	 * add it to the constraint set for this automaton */
	discrete::discrete_set::ptr temp_dset;

	for (std::vector<std::string>::size_type j = 0; j < initial_modes.size(); j++) {

		hybrid_automata::location_constraint_set temp_lcset;
		temp_lcset = location_constraint_set(new_automaton->get_id(),
				new_automaton->get_location_id(initial_modes[j]));

		// Discrete part
		if (temp_dset)
			temp_dset = discrete::compute_or_assign_union(
					temp_dset,
					discrete::discrete_set::ptr(
							new discrete::singleton_set(temp_lcset)));
		else
			temp_dset = discrete::discrete_set::ptr(
					new discrete::singleton_set(temp_lcset));
	}

	// if nothing was defined, create the universe
	if (!temp_dset) {
		temp_dset = create_universe_discrete_states();
	}
	return temp_dset;
}

/** Create a symbolic state collection from a set of discrete states and a predicate */
symbolic_state_collection::ptr create_symbolic_state_collection(
		discrete::discrete_set::ptr temp_dset, tree::node::ptr initial_pred) {
	// Continuous part
	continuous::continuous_set_ptr init_cset =
			continuous::construct_predicate_continuous_set(initial_pred);

	// Create the symbolic state
	symbolic_state::ptr init_sstate = symbolic_state::ptr(
			new symbolic_state(temp_dset, init_cset));

	// Create new symbolic_state_collection and add symbolic state
	symbolic_state_collection::ptr init_scol = symbolic_state_collection::ptr(
			new symbolic_state_collection_stl_list());

	init_scol->add(init_sstate);
	return init_scol;
}

/** Create a symbolic state collection from the modes of an automaton and a predicate */
symbolic_state_collection::ptr create_symbolic_state_collection(
		const hybrid_automata::hybrid_automaton::ptr& new_automaton,
		std::vector<std::string>& initial_modes, tree::node::ptr initial_pred) {
	discrete::discrete_set::ptr temp_dset = create_discrete_states(
			new_automaton, initial_modes);
	symbolic_state_collection::ptr scol = create_symbolic_state_collection(
			temp_dset, initial_pred);
	return scol;
}

/** Create an automaton from its ast structure and adds initial locations to the vector initial_modes */
hybrid_automata::hybrid_automaton::ptr create_automaton_from_ast(
		automaton_ast& aut_ast, std::string& aut_name, bool bind,
		symbol_table& s_table, parse_policy const& ppol) {

	std::vector<std::string> initial_modes;
	tree::node::ptr initial_pred;

	hybrid_automata::hybrid_automaton::ptr new_automaton =
			hybrid_automata::hybrid_automaton::ptr(
					new hybrid_automata::explicit_automaton(aut_name));

	std::vector<std::string> clocks;

	// Instantiate symbol table of the automaton
	for (std::vector<declaration_ast>::size_type i = 0; i
			< aut_ast.my_declaration.size(); i++) {
		std::vector<parser::symbol> s_vector = get_new_symbols(
				aut_ast.my_declaration[i], initial_pred, clocks, ppol);
		for (std::vector<parser::symbol>::size_type j = 0; j < s_vector.size(); j++)
			symbol res_symbol = instantiate_symbol(s_vector[j], s_table, bind);
	}
	// lock the symbol table
	s_table.set_locked();

	// Add symbols to the automaton
	std::vector<std::string> symbol_list = s_table.get_symbol_list(false);
	for (std::vector<std::string>::size_type k = 0; k < symbol_list.size(); k++) {
		parser::automaton_parser::add_symbol(new_automaton,
				s_table.get_symbol(symbol_list[k]));
		// std::cout << "adding symbol " << s_table.get_symbol(symbol_list[k]) << " to automaton " << new_automaton->get_name() << std::endl;
	}

	// get the initial predicate
	for (std::vector<declaration_ast>::size_type i = 0; i
			< aut_ast.my_declaration.size(); i++) {
		add_to_init_pred(aut_ast.my_declaration[i], initial_pred, s_table, ppol);
	}

	tree::node::ptr clock_pred = tree::node::null_node();

	// Create predicate t' == 1 for each clock variable "t" in the mode and add it to clock_pred
	if (!clocks.empty()) {
		for (std::vector<std::string>::size_type i = 0; i < clocks.size(); ++i) {
			// checks if clock variable is present in the symbol table
			if (s_table.is_symbol(clocks[i])) {
				variable_id v_id = variable::get_variable_id(
						s_table.get_symbol(clocks[i]).my_name);
				variable_id primed_id = variable::get_id_primedness_increased(
						v_id);
				tree::node::ptr pred_left = tree::node::ptr(
						new valuation_functions::variable_node(primed_id));
				tree::node::ptr pred_right =
						predicate_parser::create_scalar_const_node("1");
				tree::node::ptr pred_op = tree::node_ptr(
						new valuation_functions::comparison_node(EQ, pred_left,
								pred_right));

				clock_pred = valuation_functions::boolean_and(clock_pred, pred_op);
			} else
				throw std::runtime_error("Unknown clock variable " + clocks[i]);
		}
	}

	// Create and add location to the automaton for each mode
	for (std::vector<mode_ast>::size_type i = 0; i < aut_ast.my_modes.size(); ++i) {

		mode_ast mode = aut_ast.my_modes[i];

		// if initial, add it to the vector
		if (mode.initial)
			initial_modes.push_back(mode.my_name);

		tree::node::ptr inv = tree::node::null_node();
		tree::node::ptr flow = clock_pred; // initialize flow with clock_pred

		for (std::vector<dyn_ast>::size_type j = 0; j < mode.my_dyns.size(); ++j) {

			dyn_ast dyn = mode.my_dyns[j];
			tree::node::ptr pred = make_predicate_tree(dyn.my_pred, s_table,
					ppol);

			//if invariant
			if (dyn.my_type.compare("inv") == 0) {
				if (inv == tree::node::null_node())
					inv = pred;
				else
					inv = tree::node::ptr(
							new valuation_functions::boolean_node(AND, inv,
									pred));
			}

			//else flow
			else {
				if (flow == tree::node::null_node())
					flow = pred;
				else
					flow = tree::node::ptr(
							new valuation_functions::boolean_node(AND, flow,
									pred));
			}
		}

		new_automaton->add_location(
				parser::automaton_parser::create_location(mode.my_name, inv,
						flow,new_automaton->get_const_variables()));
	}

	// Create and add transitions to each mode
	for (std::vector<mode_ast>::size_type i = 0; i < aut_ast.my_modes.size(); ++i) {

		mode_ast mode = aut_ast.my_modes[i];

		for (std::vector<edge_ast>::size_type j = 0; j < mode.my_edges.size(); ++j) {

			hybrid_automata::location_id sloc = new_automaton->get_location_id(
					mode.my_name);
			hybrid_automata::location_id tloc = new_automaton->get_location_id(
					mode.my_edges[j].my_target);
			tree::node::ptr guard = make_predicate_tree(
					mode.my_edges[j].my_guard, s_table, ppol);
			tree::node::ptr update = make_predicate_tree(
					mode.my_edges[j].my_update, s_table, ppol);

			std::string label = mode.my_edges[j].my_label;

			if (s_table.is_symbol(label))
				label = s_table.get_symbol(label).my_name;

			else if (label.compare("UNLABELED") == 0)
				label = named_label::silent_name();
			//					label = aut_name+".UNLABELED";

			/**
			 * error handle : Unknown label
			 */
			else
				throw std::runtime_error("Unknown label " + label);

			// prefix label with automaton name so the label is globally unique.
			// this avoids unwanted synchronization with other automata
			// @note Label needs to be added before transition
			new_automaton->add_label(label);

			new_automaton->add_transition(
					parser::automaton_parser::create_transition(sloc, label,
							guard, update, tloc,
							new_automaton->get_controlled_variables(),
							new_automaton->get_input_variables(),
							new_automaton->get_const_variables()), false);
		}
	}

	symbolic_state_collection::ptr ini_scol = create_symbolic_state_collection(new_automaton,initial_modes,initial_pred);
	new_automaton->set_initial_states(ini_scol);

	return new_automaton;
}

hybrid_automata::hybrid_automaton_network::ptr parse_CIF(
		std::vector<hybrid_automata::hybrid_automaton_ptr>& list_automatons,
		std::istreambuf_iterator<char> beg_iter,
		std::istreambuf_iterator<char> end_iter) {
	// convert input iterator to buffered iterator
	typedef std::istreambuf_iterator<char> base_iterator_type;
	//    boost::spirit::multi_pass<base_iterator_type> first =
	//    		boost::spirit::make_default_multi_pass(base_iterator_type(beg_iter));
	//    boost::spirit::multi_pass<base_iterator_type> last =
	//    		boost::spirit::make_default_multi_pass(base_iterator_type(end_iter));

	hybrid_automata::hybrid_automaton_network::ptr aut_net;

	if (beg_iter != end_iter) {

		CIF_grammar g;

		network_ast my_network_ast;
		parser::symbol_table s_table;

		symbolic_state_collection::ptr init_scol;
		std::vector<std::string> clocks;
		std::ostringstream sout;

		// Initialize the parser policy to SX_policy
		parse_policy ppol = parse_policy::SX_policy();

		std::copy(beg_iter, end_iter, std::ostreambuf_iterator<char>(sout));
		std::string str = sout.str();

		grammar_type first = str.begin();
		grammar_type last = str.end();

		bool r = boost::spirit::qi::phrase_parse(first, last, g,
				boost::spirit::ascii::space, my_network_ast);

		aut_net = parser::automaton_parser::create_automaton_network(
				my_network_ast.my_name); // create a network of automata

		if (r && first == last) {

			tree::node::ptr init_pred;
			discrete::discrete_set::ptr init_dset;

			s_table.set_context(my_network_ast.my_name);

			// Instantiate symbol table of the network
			for (std::vector<declaration_ast>::size_type i = 0; i
					< my_network_ast.my_declaration.size(); i++) {

				std::vector<parser::symbol> s_vector = get_new_symbols(
						my_network_ast.my_declaration[i], init_pred, clocks,
						ppol);
				for (std::vector<parser::symbol>::size_type j = 0; j
						< s_vector.size(); j++)
					symbol res_symbol = instantiate_symbol(s_vector[j],
							s_table, false);
			}

			if (my_network_ast.my_automata_dec.empty()) {

				for (std::vector<automaton_ast>::size_type i = 0; i
						< my_network_ast.my_automata.size(); ++i) {

					hybrid_automata::hybrid_automaton::ptr new_aut;
					std::vector<std::string> init_modes; // initial mode of an automaton
					std::string bind_name =
							my_network_ast.my_automata[i].my_name;
					std::string aut_name = my_network_ast.my_name + "."
							+ bind_name;

					/** New symbol table for the instantiation. */
					// symbol_table aut_s_table(aut_name);
					// Note: The symbol table for the automaton needs to have the same variables
					// as those in the network, because in the CIF format variables can
					// be declared in the network, and are then available in all subcomponents
					symbol_table aut_s_table = s_table;
					aut_s_table.set_context(aut_name);

//					// Add symbols from s_table to aut_s_table
//					std::vector<std::string> s_list = s_table.get_symbol_list(
//							false);
//
//					for (std::vector<std::string>::size_type j = 0; j
//							< s_list.size(); ++j)
//						aut_s_table.add_symbol(s_list[j],
//								s_table.get_symbol(s_list[j]));

					new_aut = create_automaton_from_ast(
							my_network_ast.my_automata[i], aut_name, true,
							aut_s_table, ppol);
					list_automatons.push_back(new_aut);

					aut_net = aut_net->compute_or_assign_composition(new_aut);
					if (init_scol) {
						init_scol->intersection_assign(new_aut->get_initial_states());
					} else
						init_scol = symbolic_state_collection::ptr(new_aut->get_initial_states()->clone());

				} // end for
			} // end if(my_network_ast.my_automata_dec.empty())

			else {

				for (std::vector<automaton_dec_ast>::size_type i = 0; i
						< my_network_ast.my_automata_dec.size(); ++i) {

					hybrid_automata::hybrid_automaton::ptr new_aut;
					std::vector<std::string> init_modes; // initial mode of an automaton
					std::string bind_name =
							my_network_ast.my_automata_dec[i].my_instance_name;
					std::string aut_name = my_network_ast.my_name + "."
							+ bind_name;

					/** New symbol table for the instantiation. */
					symbol_table aut_s_table(aut_name);

					automaton_ast automaton = my_network_ast.get_automaton_ast(
							my_network_ast.my_automata_dec[i].my_def_name);

					// Get the keys
					std::vector<std::string> my_keys = get_keys(
							automaton.my_declaration);

					// Mapping between key and its value is positional
					for (std::vector<std::string>::size_type j = 0; j
							< my_keys.size(); ++j) {
						std::string
								value =
										my_network_ast.my_automata_dec[i].my_parameters[j].my_exp;
						aut_s_table.add_symbol(my_keys[j],
								s_table.get_symbol(value));
					}

					new_aut = create_automaton_from_ast(automaton, aut_name,
							true, aut_s_table, ppol);
					list_automatons.push_back(new_aut);

					aut_net = aut_net->compute_or_assign_composition(new_aut);
					if (init_scol) {
						init_scol->intersection_assign(new_aut->get_initial_states());
					} else
						init_scol = symbolic_state_collection::ptr(new_aut->get_initial_states()->clone());

				} // end for

			} // end else

			// compute the initial predicate
			// Note: This needs to be done after all automata have been instantiated, so that
			// the variables are instantiated! Otherwise, we wouldn't be allowed to use them.
			for (std::vector<declaration_ast>::size_type i = 0; i
					< my_network_ast.my_declaration.size(); i++) {
				add_to_init_pred(my_network_ast.my_declaration[i], init_pred, s_table, ppol);
			}
			// if there are constraints on the initial states of the network,
			// add them to the initial states
			if (init_pred) {
				discrete::discrete_set::ptr temp_dset = create_universe_discrete_states();
				symbolic_state_collection::ptr temp_scol;
				// Create a symbolic state collection from the initial predicate:
				// the null pointer temp_dset will be replaced by the universe constraint
				temp_scol = create_symbolic_state_collection(temp_dset,
						init_pred);
				if (init_scol) {
					init_scol->intersection_assign(temp_scol);
				} else
					init_scol = temp_scol;
			}
			aut_net->set_initial_states(init_scol);

		} // end if (r && first==last)
		else {
			/**
			 * error handle : precise as to when the parser fails
			 * (the location fails)
			 */
			throw std::runtime_error(
					"Could not parse CIF file from '"
							+ std::string(first, last) + "'.");
			//throw std::runtime_error("Could not parse CIF file from '" + std::string(beg_iter, end_iter) + "'.");
		}
	} //end if(beg_iter!=end_iter)

	if (aut_net) {
		list_automatons.push_back(aut_net);
	}

	return aut_net;
}

void parse_CIF(std::string source,
		std::vector<hybrid_automata::hybrid_automaton_ptr>& list_automatons,
		std::string component_to_create) {

	std::ifstream f(source.c_str(), std::ifstream::in);

	hybrid_automata::hybrid_automaton_network::ptr a =
			CIF_predicate_parser::parse_CIF(list_automatons,
					std::istreambuf_iterator<char>(f),
					std::istreambuf_iterator<char>());

	//std::cout << " double checking ini states: " << a->get_initial_states() << std::endl;
}

}
