/*
 * symbol_table_operators.h
 *
 *  Created on: Mar 1, 2011
 *      Author: frehse
 */

#ifndef SYMBOL_TABLE_OPERATORS_H_
#define SYMBOL_TABLE_OPERATORS_H_

#include <boost/variant.hpp>
#include "symbol_table.h"

namespace parser {

/**
 * Take dimension in a string and return a int.
 * If the dimension is a token, replace it, when we know its value (use symbol_table)
 * 1 is default value.
 */
unsigned int instantiate_dimension(std::string const& dim,
		symbol_table s_table);

/**
 * Instantiation of a symbol according to a given symbol table
 *
 * If the symbol is not yet in the table, it is added.
 * If the token name already exists in the symbol table,
 * the symbol is checked for compatibility. If the symbol table
 * maps the token to a constant value, the value is checked
 * for compatibility and instantiated as a matrix of corresponding
 * dimensions.
 *  *
 * If bind is true, then the symbol must be mapped by list_arguments.
 */
symbol instantiate_symbol(const symbol& unv_symbol, symbol_table& s_table, bool bind);

}

#endif /* SYMBOL_TABLE_OPERATORS_H_ */
