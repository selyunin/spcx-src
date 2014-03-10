/*
 * symbol_table_operators.cpp
 *
 *  Created on: Mar 1, 2011
 *      Author: frehse
 */

#include "symbol_table_operators.h"
#include "parser_basics.h"

namespace parser {

unsigned int instantiate_dimension(std::string const& dim,
		symbol_table s_table) {
	unsigned int dim_int;
	if (!is_number(dim)) {
		if (s_table.is_symbol(dim)) {
//std::cout << "searching for dimension " << dim << " in symbol table "<< s_table << std::endl;
			symbol& sym = s_table.get_symbol(dim);
			std::string s_d1;
			/** If the symbol is attributed a constant value,
			 * take this value. Otherwise, take dimension 1.
			 */
			if (sym.my_symbol_type==symbol::CONST_VALUE) {
				s_d1=boost::any_cast<std::string>(sym.my_value(0,0));
			}
			if (is_number(s_d1)) {
				dim_int = from_string<unsigned int> (s_d1);
			} else
				dim_int = 1;
		} else
			dim_int = 1;
	} else {
		dim_int = from_string<unsigned int> (dim);
	}
	return dim_int;
}

symbol instantiate_symbol(
		const symbol& unv_symbol,
		symbol_table& s_table, bool bind ) {

	symbol res_symbol;

	std::string token_name = unv_symbol.my_name;
	bool local = unv_symbol.is_local;
	bool controlled = unv_symbol.is_controlled;
	symbol::symbol_type s_type = unv_symbol.my_symbol_type;
	symbol::data_type d_type = unv_symbol.my_data_type;
	symbol::dynamics_type dyn_type = unv_symbol.my_dynamics_type;
	unsigned int dim1 = unv_symbol.dim1;
	unsigned int dim2 = unv_symbol.dim2;

	std::string local_name(s_table.get_context() + "." + token_name);

	bool is_mapped = s_table.is_symbol(token_name);
	bool is_symbol = false;

	/**
	 * Verifications
	 */
	if (bind) {
		if (local) {
			if (is_mapped) {
				std::string s=s_table.get_symbol(token_name).my_name;
				throw std::runtime_error("Cannot map local parameter " + token_name
						+ " to " + s + " in component " + s_table.get_context()
						+ ". Local parameters can not be mapped.");
			}

			/** Add local (new) symbol to symbol table */
			res_symbol = unv_symbol;
			res_symbol.my_name = local_name;
			s_table.add_symbol(token_name,res_symbol);
		} else {
			if (!is_mapped)
				throw std::runtime_error(
						"Parameter " + token_name
								+ " must be mapped to either a value or another parameter.");

			/** Check compatibility of mapped symbol */
			symbol& mapped_symbol = s_table.get_symbol(token_name);
			if (mapped_symbol.my_symbol_type == symbol::CONST_VALUE) {
				/** The symbol is mapped to numeric values, so attribute
				 * them to the matrix my_value.
				 */
				std::string symbol_value = boost::any_cast<std::string>(mapped_symbol.my_value(0, 0));

				try {
					math::matrix<boost::any> values(dim1,dim2);
					unsigned int i = 0;
					std::string num;
				    std::stringstream valueStream(symbol_value);
					while (getline(valueStream, num, ' ')) {
						// If it's a vector, treat as column vector instead
						// of row vector
						if (dim2 > 1)
							values(i/dim2,i%dim1) = boost::any(num);
						else
							values(i,0) = boost::any(num);
						++i;
					}
					if (i != dim1*dim2)
						throw basic_exception("Expected "+to_string(dim1*dim2)+" values, found "+to_string(i)+".");

					res_symbol = unv_symbol;
					res_symbol.my_symbol_type = symbol::CONST_VALUE;
					res_symbol.my_value = values;
					mapped_symbol = res_symbol;
//std::cout << "affected symbol " << mapped_symbol << std::endl;
				} catch ( std::exception& e ) {
					if (dim1 > 1 || dim1 > 1)
						throw basic_exception("Can't convert \"" + symbol_value
								+ "\" to (" + to_string(dim1) + "," + to_string(dim2) + ") matrix.",
								e);
					else
						throw basic_exception("Can't convert \"" + symbol_value
								+ "\" to scalar value.",
								e);
				}
			} else {
				if (mapped_symbol.dim1 != dim1 || mapped_symbol.dim2 != dim2)
					throw basic_exception("Cannot map parameter " + token_name
							+ " to " + mapped_symbol.my_name
							+ " because they have different dimensions.");

				if (mapped_symbol.my_data_type != d_type)
					throw basic_exception("Cannot map parameter " + token_name
							+ " to " + mapped_symbol.my_name
							+ " because they have different data types.");
				if (!mapped_symbol.is_controlled && controlled)
					throw basic_exception("Cannot map controlled parameter "
							+ token_name + " to uncontrolled parameter "
							+ mapped_symbol.my_name + ".");

				/** The token mapped to a symbol */
				res_symbol = mapped_symbol;

				/** Adopt the controlled status of the component
				 * Even if the variable is controlled in the network, might
				 * still be uncontrolled in one of the subcomponents,
				 * in which it must be instantiated as uncontrolled. */
				res_symbol.is_controlled = controlled;
			}
		}
	} else {
		if (is_mapped)
			throw std::runtime_error("Cannot map parameter " + token_name
					+ ". Check if parameter is listed twice. Note that only bind parameters can be mapped to other parameters.");

		res_symbol = unv_symbol;
		res_symbol.my_name = local_name;
		s_table.add_symbol(token_name,res_symbol);
	}
	return res_symbol;
}

}
