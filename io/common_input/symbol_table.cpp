/*
 * symbol_table.cpp
 *
 *  Created on: May 6, 2010
 *      Author: gvincent
 */

#include "symbol_table.h"
#include <stdexcept>

#include "utility/stl_helper_functions.h"

namespace parser {

symbol::symbol(std::string const& name, symbol_type type,
		data_type number_type, dynamics_type dyn_type, unsigned int const& d1, unsigned int const& d2,
		math::matrix<boost::any> const& value, bool local, bool controlled) {
	my_name = name;
	my_symbol_type = type;
	my_data_type = number_type;
	my_dynamics_type = dyn_type;
	dim1 = d1;
	dim2 = d2;
	my_value = value;
	is_local = local;
	is_controlled = controlled;
}

symbol::symbol(std::string const& name, symbol_type type,
		data_type number_type, dynamics_type dyn_type,
		boost::any const& value, bool local, bool controlled) {
	my_name = name;
	my_symbol_type = type;
	my_data_type = number_type;
	my_dynamics_type = dyn_type;
	dim1 = 1;
	dim2 = 1;
	my_value = math::matrix<boost::any>(1, 1, value);
	is_local = local;
	is_controlled = controlled;
}

symbol::symbol(std::string const& name, symbol const& s) {
	my_name = name;
	my_symbol_type = s.my_symbol_type;
	my_data_type = s.my_data_type;
	dim1 = s.dim1;
	dim2 = s.dim2;
	my_value = s.my_value;
	is_local = s.is_local;
	is_controlled = s.is_controlled;
}

symbol::symbol() {
	my_name = "";
	my_symbol_type = CONST_VALUE;
	my_data_type = REAL;
	my_dynamics_type = CONSTANT;
	dim1 = 1;
	dim2 = 1;
	my_value = math::matrix<boost::any>();
	is_local = false;
	is_controlled = true;
}

void symbol::print(std::ostream& os) const {
	os << my_name << ": ";
	if (my_symbol_type == VARIABLE)
		os << "var ";
	else if (my_symbol_type == LABEL)
		os << "label ";
	else if (my_symbol_type == CONST_VALUE)
		os << "const value ";
	else
		os << "unknown symbol type ";
	if (my_dynamics_type == CONSTANT)
		os << "const ";
	else if (my_dynamics_type == ANY)
		os << "any ";
	else if (my_dynamics_type == EXPLICIT)
		os << "explicit ";
	else
		os << "unknown dynamics ";
	if (my_data_type == INT)
		os << "int ";
	else if (my_data_type == REAL)
		os << "real ";
	else if (my_data_type == UNINTERPR_STRING)
		os << "uninterpreted ";
	else
		os << "unknown data type ";
	os << "(" << dim1 << "," << dim2 << ") ";
	if (is_local)
		os << "local ";
	if (is_controlled)
		os << "controlled ";
	else
		os << "uncontrolled ";

	if (my_symbol_type == CONST_VALUE && my_data_type == UNINTERPR_STRING) {
		std::string symbol_value = boost::any_cast<std::string>(my_value(0, 0));
		os << "value: " << std::flush;
		os << symbol_value;
	}
}

bool symbol::operator==(const symbol& s) const {
	// if there is a constant value attributed, we
	// don't know how to compare (because we don't know the
	// data type. So we report the symbols as being
	// different.
	if ((my_value.size1()>0) || (s.my_value.size1()>0))
		return false;
	return my_name == s.my_name && my_data_type == s.my_data_type
			&& my_dynamics_type == s.my_dynamics_type && my_symbol_type
			== s.my_symbol_type && dim1 == s.dim1 && dim2 == s.dim2 && is_local
			== s.is_local && is_controlled == s.is_controlled;
}

bool symbol::operator!=(const symbol& s) const {
	return !operator==(s);
}

std::ostream& operator<<(std::ostream& os, const symbol& s) {
	s.print(os);
	return os;
}

symbol_table::symbol_table(const std::string& its_context) : context(its_context),my_locked(false) {
}

void symbol_table::add_symbol(symbol const& symb) {
	if (is_symbol(symb.my_name)) {
		//  && my_symbol_table[symb.my_name]!=symb)
		std::stringstream ss;
		ss << "existing symbol: " << my_symbol_table[symb.my_name] << std::endl;
		ss << "symbol to add  : " << symb << std::endl;
		throw std::runtime_error("Symbol "+symb.my_name+ " exists already, with different definition.\n" + ss.str());
	}
	else
		my_symbol_table[symb.my_name] = symb;
}

void symbol_table::add_symbol(std::string const& token_name, symbol const& symb) {
	if (is_symbol(token_name)) {
		//  && my_symbol_table[token_name]!=symb
		std::stringstream ss;
		ss << "existing symbol: " << my_symbol_table[token_name] << std::endl;
		ss << "symbol to add  : " << symb << std::endl;
		throw std::runtime_error("Symbol "+token_name+" exists already, with different definition.\n" + ss.str());
	}
	else
		my_symbol_table[token_name] = symb;
}

const boost::any& symbol_table::get_value(std::string const& name,
		std::string const& context_error) {
	throw_if_not_symbol(name, context_error);
	if (my_symbol_table[name].dim1 > 1 || my_symbol_table[name].dim2 > 1)
		throw std::runtime_error(name + " is a matrix. Specify coordinates"
				+ context_error);
	else {
		return my_symbol_table[name].my_value(0, 0);
	}
}

const symbol& symbol_table::get_symbol(std::string const& name) const {
	throw_if_not_symbol(name);
	return my_symbol_table.find(name)->second;
}

symbol& symbol_table::get_symbol(std::string const& name) {
	throw_if_not_symbol(name);
	return my_symbol_table[name];
}

const boost::any& symbol_table::get_value(std::string const& name,
		unsigned int const& d1, unsigned int const& d2,
		std::string const& context_error) {
	throw_if_not_symbol(name, context_error);
	if (my_symbol_table[name].dim1 < d1 || my_symbol_table[name].dim2 < d2)
		throw std::runtime_error("exceeded limits of the matrix : " + name
				+ context_error);
	else
		return my_symbol_table[name].my_value((d1 - 1), (d2 - 1));
}

symbol::symbol_type symbol_table::get_type(std::string const& name) {
	throw_if_not_symbol(name);
	return my_symbol_table[name].my_symbol_type;
}

void symbol_table::set_type(std::string const& name, symbol::symbol_type type) {
	throw_if_not_symbol(name);
	my_symbol_table[name].my_symbol_type = type;
}

symbol::data_type symbol_table::get_data_type(std::string const& name) {
	throw_if_not_symbol(name);
	return my_symbol_table[name].my_data_type;
}

bool symbol_table::get_local(std::string const& name) {
	throw_if_not_symbol(name);
	return my_symbol_table[name].is_local;
}

bool symbol_table::get_controlled(std::string const& name) {
	throw_if_not_symbol(name);
	return my_symbol_table[name].is_controlled;
}

unsigned int symbol_table::get_dim1(std::string const& name) {
	throw_if_not_symbol(name);
	return my_symbol_table[name].dim1;
}

unsigned int symbol_table::get_dim2(std::string const& name) {
	throw_if_not_symbol(name);
	return my_symbol_table[name].dim2;
}

void symbol_table::set_value(std::string const& name, boost::any const& value) {
	throw_if_not_symbol(name);
	if (my_symbol_table[name].dim1 > 1 || my_symbol_table[name].dim2 > 1)
		throw std::runtime_error(name + "is a matrix. Precise coordinates");
	else {
		my_symbol_table[name].my_value(0, 0) = value;
	}
}

void symbol_table::set_value(std::string const& name, unsigned int const& d1,
		unsigned int const& d2, boost::any const& value) {
	throw_if_not_symbol(name);
	if (my_symbol_table[name].dim1 < d1 || my_symbol_table[name].dim2 < d2)
		throw std::runtime_error("exceeded limits of the matrix : " + name);
	else {
		my_symbol_table[name].my_value((d1 - 1), (d2 - 1)) = value;
	}
}

void symbol_table::delete_symbol(std::string const& name) {
	my_symbol_table.erase(name);
}

bool symbol_table::is_symbol(std::string const& name) const {
	if (my_symbol_table.find(name) != my_symbol_table.end())
		return true;
	else
		return false;
}

std::vector<std::string> symbol_table::get_symbol_list(bool typed,
		symbol::symbol_type typ) {
	std::vector<std::string> list_var;
	for (std::map<std::string, symbol>::const_iterator it =
			my_symbol_table.begin(); it != my_symbol_table.end(); ++it) {
		if (typed) {
			if (it->second.my_symbol_type == typ)
				list_var.push_back(it->first);
		} else
			list_var.push_back(it->first);
	}
	return list_var;
}

std::string symbol_table::get_context() const {
	return context;
}

void symbol_table::set_context(std::string const& s) {
	context = s;
}

void symbol_table::throw_if_not_symbol(std::string const& s, std::string const& context_msg) const {
	if (!is_symbol(s)) {
		throw basic_exception(
				"Unknown symbol " + s + context_msg
						+ ".\n Note that all symbols need to be declared as formal parameters of a component.");
	}
}

bool symbol_table::is_locked() const {
	return my_locked;
}

void symbol_table::set_locked(bool b) {
	my_locked = b;
}


std::ostream& operator<<(std::ostream& os, const symbol_table& s) {
	os << "context: " << s.get_context();
	os << std::endl << "symbols: " << std::endl;
	for (std::map<std::string, symbol>::const_iterator it =
			s.my_symbol_table.begin(); it != s.my_symbol_table.end(); ++it) {
		os << it->first << " -> " << it->second << std::endl;
	}
	return os;
}

}
