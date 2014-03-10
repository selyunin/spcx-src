/*
 * CIF_parser.h
 *
 *  Created on: Jul 3, 2009
 *      Author: gvincent
 */

#ifndef CIF_PARSER_H_
#define CIF_PARSER_H_

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/shared_ptr.hpp>

#include "io/common_input/automaton_creation.h"
#include "core/hybrid_automata/explicit_automaton.h"
#include "io/CIF_format/CIF_predicate_grammar.h"
#include "io/CIF_format/CIF_predicate_parser.h"
#include "utility/calc_string.h"
#include "core/predicates/node_print_visitor.h"
#include "core/hybrid_automata/hybrid_automaton_network.h"
#include "io/common_input/network_creation.h"

namespace CIF_predicate_parser {

struct dyn_ast {
	std::string my_type; // type (inv or flow)
	predicate_ast my_pred;
	dyn_ast& operator+=(std::string const& inv_or_flow);
	dyn_ast& operator-=(predicate_ast const& pred);
};

struct edge_ast {

	edge_ast() : my_guard(), my_update(), my_target(), my_label("UNLABELED"){}

	edge_ast& operator+=(predicate_ast const& pred); //set my_guard
	edge_ast& operator-=(predicate_ast const& pred); //set my_update
	edge_ast& operator*=(std::string const& target); //set my_target
	edge_ast& operator>=(std::string const& label); //set my_label

	predicate_ast my_guard;
	predicate_ast my_update;
	std::string my_target;
	std::string my_label;
};

struct mode_ast{

	std::string my_name;
	std::vector<dyn_ast> my_dyns;
	std::vector<edge_ast> my_edges;
	bool initial;

	mode_ast& operator+=(std::string const& n);
	mode_ast& operator-=(std::vector<dyn_ast> const& dyns);
	mode_ast& operator*=(std::vector<edge_ast> const& edges);

};

struct expr_ast{

	std::string my_exp;

	expr_ast& operator+=(std::string const& param_name);
	expr_ast& operator-=(boost::variant<int, double>& param_val);
};

struct declaration_ast{

	std::string my_dec;
	bool is_local;

	declaration_ast& operator+=(std::string const& type);
	declaration_ast& operator-=(std::vector<expr_ast> const& expression);
};

struct automaton_ast{

	automaton_ast& operator+=(std::vector<mode_ast> const& modes);
	automaton_ast& operator-=(std::string const& n);
	automaton_ast& operator>=(std::vector<declaration_ast>& declaration);
	automaton_ast& operator*=(std::vector<declaration_ast>& declaration);

	std::string my_name;
	std::vector<mode_ast> my_modes;
	std::vector<declaration_ast> my_declaration;
};

struct automaton_dec_ast{

	std::string my_instance_name;
	std::string my_def_name;
	std::vector<expr_ast> my_parameters;

	automaton_dec_ast& operator+=(std::string const& instance_name);
	automaton_dec_ast operator-=(std::string const& def_name);
	automaton_dec_ast operator*=(std::vector<expr_ast> const& parameters);
};

struct network_ast{

	std::vector<automaton_ast> my_automata;
	std::string my_name;
	std::vector<declaration_ast> my_declaration;
	std::vector<automaton_dec_ast> my_automata_dec;

	network_ast& operator+=(std::vector<automaton_ast> const& automata);
	network_ast& operator-=(std::string const& name);
	network_ast& operator*=(std::vector<declaration_ast> const& declaration);
	network_ast& operator>=(std::vector<automaton_dec_ast> const& automata_dec);

	automaton_ast get_automaton_ast(std::string def_name)
	{
		std::vector<automaton_ast>::size_type i;
		for( i = 0; i < my_automata.size(); i++){
			if(my_automata[i].my_name.compare(def_name) == 0)
				break;
		}

		if(i>my_automata.size())
			throw std::runtime_error("Could not find automaton " + def_name);

		return my_automata[i];
	}
};

}

BOOST_FUSION_ADAPT_STRUCT(
		CIF_predicate_parser::expr_ast,
		(std::string, my_exp)
)

BOOST_FUSION_ADAPT_STRUCT(
	CIF_predicate_parser::edge_ast,
	(CIF_predicate_parser::predicate_ast, my_guard)
	(CIF_predicate_parser::predicate_ast, my_update)
	(std::string, my_target)
	(std::string, my_label)
)

BOOST_FUSION_ADAPT_STRUCT(
	CIF_predicate_parser::dyn_ast,
    (std::string, my_type)
    (CIF_predicate_parser::predicate_ast, my_pred)
)

BOOST_FUSION_ADAPT_STRUCT(
		CIF_predicate_parser::declaration_ast,
		(std::string, my_dec)
		(bool, is_local)
)

BOOST_FUSION_ADAPT_STRUCT(
	CIF_predicate_parser::mode_ast,
    (std::string, my_name)
    (std::vector<CIF_predicate_parser::dyn_ast>, my_dyns)
    (std::vector<CIF_predicate_parser::edge_ast>, my_edges)
    (bool, initial)
)

BOOST_FUSION_ADAPT_STRUCT(
	CIF_predicate_parser::automaton_dec_ast,
	(std::string, my_instance_name)
	(std::string, my_def_name)
	(std::vector<CIF_predicate_parser::expr_ast>, my_parameters)
)

BOOST_FUSION_ADAPT_STRUCT(
	CIF_predicate_parser::automaton_ast,
	(std::string, my_name)
	(std::vector<CIF_predicate_parser::mode_ast>, my_modes)
	(std::vector<CIF_predicate_parser::declaration_ast>, my_declaration)
)

BOOST_FUSION_ADAPT_STRUCT(
	CIF_predicate_parser::network_ast,
	(std::vector<CIF_predicate_parser::automaton_ast>, my_automata)
	(std::string, my_name)
	(std::vector<CIF_predicate_parser::declaration_ast>, my_declaration)
	(std::vector<CIF_predicate_parser::automaton_dec_ast>, my_automata_dec)
)

namespace CIF_predicate_parser {

/** Create a network structure from the input string
 *
 * All generated automata are stored in list_automatons. */
hybrid_automata::hybrid_automaton_network::ptr parse_CIF(std::vector<
		hybrid_automata::hybrid_automaton_ptr>& list_automatons, std::istreambuf_iterator<char> beg_iter,std::istreambuf_iterator<char> end_iter);

/** Read the CIF file with name source, and return the vector of all created automata in list_automatons.
 *
 * @todo: component_to_create is ignored for now. All automata are created */
void parse_CIF(std::string source, std::vector<
		hybrid_automata::hybrid_automaton_ptr>& list_automatons, std::string component_to_create);

//  CIF grammar
struct CIF_definition : predicate_definition
{
	CIF_definition();

	qi::rule<grammar_type, boost::variant<int, double>(), ascii::space_type> value;
	qi::rule<grammar_type, std::string(), ascii::space_type> inv_flow_token;
	qi::rule<grammar_type, std::string(), ascii::space_type> static_type;
	qi::rule<grammar_type, std::string(), ascii::space_type> param_type;
	qi::rule<grammar_type, std::string(), ascii::space_type> controlled;
	qi::rule<grammar_type, std::string(), ascii::space_type> initial;
	qi::rule<grammar_type, expr_ast(), ascii::space_type> dec_expression;
	qi::rule<grammar_type, std::vector<expr_ast>(),ascii::space_type> dec_expressions;
	qi::rule<grammar_type, declaration_ast(), ascii::space_type> declaration;
	qi::rule<grammar_type, std::vector<declaration_ast>(),ascii::space_type> declarations;

	qi::rule<grammar_type, predicate_ast(), ascii::space_type > update, guard;
	qi::rule<grammar_type, dyn_ast(), ascii::space_type > dyn;
	qi::rule<grammar_type, std::vector<dyn_ast>(), ascii::space_type > dyns;
	qi::rule<grammar_type, edge_ast(), ascii::space_type > edge;
	qi::rule<grammar_type, std::vector<edge_ast>(), ascii::space_type > edges;
	qi::rule<grammar_type, mode_ast(), ascii::space_type > mode;
	qi::rule<grammar_type, std::vector<mode_ast>(), ascii::space_type > modes;
	qi::rule<grammar_type, automaton_dec_ast(), ascii::space_type > automaton_dec;
	qi::rule<grammar_type, std::vector<automaton_dec_ast>(), ascii::space_type > automata_decs;
	qi::rule<grammar_type, automaton_ast(), ascii::space_type > automaton;
	qi::rule<grammar_type, std::vector<automaton_ast>(), ascii::space_type > automata;
	qi::rule<grammar_type, std::vector<automaton_ast>(), ascii::space_type > automata_defs;
	qi::rule<grammar_type, network_ast(), ascii::space_type > network;
};

struct CIF_grammar : qi::grammar<grammar_type, network_ast(), ascii::space_type>, CIF_definition
    {

		CIF_grammar() : CIF_grammar::base_type(network){};
    };

}
#endif /* CIF_PARSER_H_ */
