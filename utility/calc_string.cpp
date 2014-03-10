#include "utility/calc_string.h"
//#include "calc_string_operators.h"
#include "utility/stl_helper_functions.h"

#include <iostream>

calc_string::calc_string() {
	my_string = "";
}

calc_string::calc_string(double b) {
	my_string = to_string(b);
}

calc_string::calc_string(int b) {
	my_string = to_string(b);
}

calc_string::calc_string(bool b) {
	my_string = to_string(b);
}

calc_string::calc_string(const std::string& s) {
	my_string = s;
}

calc_string::calc_string(const calc_string& cs) {
	my_string = cs.get_my_string();
}

calc_string::~calc_string() {
}
;

const std::string& calc_string::get_my_string() const {
	return my_string;
}

bool calc_string::parentheses = false;

bool calc_string::get_parentheses() {
	return calc_string::parentheses;
}

void calc_string::set_parentheses(bool b) {
	calc_string::parentheses = b;
}



