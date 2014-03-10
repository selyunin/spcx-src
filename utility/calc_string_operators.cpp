/*
 * calc_string_operators.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: frehse
 */

#include "utility/calc_string_operators.h"

// Operators
calc_string operator<(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " < "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " < " + s2.get_my_string());
}

calc_string operator<=(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " <= "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " <= " + s2.get_my_string());
}

calc_string operator>(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " > "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " > " + s2.get_my_string());
}

calc_string operator>=(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " >= "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " >= " + s2.get_my_string());
}

calc_string operator==(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " == "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " == " + s2.get_my_string());
}

calc_string operator!=(const calc_string& s1, const calc_string& s2) {
	return !(s1 == s2);
}

calc_string operator&&(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " && "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " && " + s2.get_my_string());
}

calc_string operator||(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " || "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " || " + s2.get_my_string());
}

calc_string operator!(const calc_string& s) {
	if (calc_string::get_parentheses() == true)
		return calc_string("( !" + s.get_my_string() + ")");
	else
		return calc_string(" !" + s.get_my_string());
}

calc_string operator-(const calc_string& s) {
	if (calc_string::get_parentheses() == true)
		return calc_string("( -" + s.get_my_string() + ")");
	else
		return calc_string(" -" + s.get_my_string());
}
calc_string operator/(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " / "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " / " + s2.get_my_string());
}

calc_string operator*(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " * "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " * " + s2.get_my_string());
}

calc_string operator-(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " - "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " - " + s2.get_my_string());
}

calc_string operator+(const calc_string& s1, const calc_string& s2) {
	if (calc_string::get_parentheses() == true)
		return calc_string("(" + s1.get_my_string() + " + "
				+ s2.get_my_string() + ")");
	else
		return calc_string(s1.get_my_string() + " + " + s2.get_my_string());
}
