/*
 * symbol_table.h
 *
 *  Created on: Apr 22, 2010
 *      Author: gvincent
 */

#ifndef SYMBOL_TABLE_H_
#define SYMBOL_TABLE_H_
#include <boost/any.hpp>
#include <string>
#include <map>
#include <vector>
#include "math/matrix.h"

namespace parser {

class symbol {
public:
	typedef enum {
		CONSTANT, ANY, EXPLICIT
	} dynamics_type;
	typedef enum {
		REAL, INT, UNINTERPR_STRING
	} data_type;
	typedef enum {
		VARIABLE, LABEL, CONST_VALUE
	} symbol_type;

	/** Create a symbol from a scalar value.
	 *
	 * The symbol is of dimension (1,1).
	 * The matrix my_value is of
	 * dimension (1,1) so that it can store a
	 * constant value if necessary.
	 */
	symbol(std::string const& name, symbol_type type, data_type number_type,
			dynamics_type dyn_type, boost::any const& value,
			bool local = false, bool controlled = true);

	/** Create a symbol from a matrix value.
	 *
	 * The value of the symbol is the (d1,d2) matrix value.
	 */
	symbol(std::string const& name, symbol_type type, data_type number_type,
			dynamics_type dyn_type, unsigned int const& d1,
			unsigned int const& d2, math::matrix<boost::any> const& value,
			bool local = false, bool controlled = true);

	/** Copy all properties of s except the name. */
	symbol(std::string const& name, symbol const& s);

	symbol();

	/** Equality test */
	bool operator==(const symbol& s) const;

	/** Inequality test */
	bool operator!=(const symbol& s) const;

	/** Output as ASCII stream */
	void print(std::ostream& os) const;

	std::string my_name;
	data_type my_data_type;
	dynamics_type my_dynamics_type;
	symbol_type my_symbol_type;
	unsigned int dim1, dim2;
	math::matrix<boost::any> my_value;
	bool is_local, is_controlled;
};

/** Output symbol to stream
 */
std::ostream& operator<<(std::ostream& os, const symbol& s);

/** The symbol_table attributes to a token a symbol. */
class symbol_table {
public:

	symbol_table(const std::string& its_context = "");

	/** Associate a symbol with a token of the same name. */
	void add_symbol(symbol const& symb);

	/** Associate a symbol with a token. */
	void add_symbol(std::string const& token_name, symbol const& symb);

	const boost::any& get_value(std::string const& name,
			std::string const& context_error = "");

	const symbol& get_symbol(std::string const& name) const;
	symbol& get_symbol(std::string const& name);

	const boost::any& get_value(std::string const& name,
			unsigned int const& d1, unsigned int const& d2,
			std::string const& context_error = "");

	symbol::symbol_type get_type(std::string const& name);
	symbol::data_type get_data_type(std::string const& name);

	bool get_local(std::string const& name);
	bool get_controlled(std::string const& name);

	unsigned int get_dim1(std::string const& name);
	unsigned int get_dim2(std::string const& name);

	void set_type(std::string const& name, symbol::symbol_type type);

	void set_value(std::string const& name, boost::any const& value);

	void set_value(std::string const& name, unsigned int const& d1,
			unsigned int const& d2, boost::any const& value);

	void delete_symbol(std::string const& name);

	bool is_symbol(std::string const& name) const;

	/** Returns the symbols in the symbol table.
	 *
	 * If typed is false, all symbols are returned.
	 * Otherwise, only symbols of the given type are returned
	 * (default type is variable). */
	std::vector<std::string> get_symbol_list(bool typed, symbol::symbol_type typ = symbol::VARIABLE);

	std::string get_context() const;

	void set_context(std::string const& s);

	/** Throws if there is no symbol with name s in the symbol table. */
	void throw_if_not_symbol(std::string const& s,
			std::string const& context_msg = "") const;

	bool is_locked() const;
	void set_locked(bool b = true);

protected:
	std::map<std::string, symbol> my_symbol_table;
	std::string context;
	bool my_locked;

	friend std::ostream& operator<<(std::ostream& os, const symbol_table& s);
};

/** Output symbol to stream
 */
std::ostream& operator<<(std::ostream& os, const symbol_table& s);



}

#endif /* SYMBOL_TABLE_H_ */
