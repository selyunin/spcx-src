/*
 * calc_string_operators.h
 *
 *  Created on: Nov 11, 2009
 *      Author: frehse
 */

#ifndef CALC_STRING_OPERATORS_H_
#define CALC_STRING_OPERATORS_H_

#include "utility/calc_string.h"

// Operators
calc_string operator<(const calc_string& s1, const calc_string& s2);
calc_string operator<=(const calc_string& s1, const calc_string& s2);
calc_string operator>(const calc_string& s1, const calc_string& s2);
calc_string operator>=(const calc_string& s1, const calc_string& s2);
calc_string operator==(const calc_string& s1, const calc_string& s2);
calc_string operator!=(const calc_string& s1, const calc_string& s2);

calc_string operator&&(const calc_string& s1, const calc_string& s2);
calc_string operator||(const calc_string& s1, const calc_string& s2);
calc_string operator!(const calc_string& s);
calc_string operator-(const calc_string& s);
calc_string operator/(const calc_string& s1, const calc_string& s2);
calc_string operator*(const calc_string& s1, const calc_string& s2);
calc_string operator-(const calc_string& s1, const calc_string& s2);
calc_string operator+(const calc_string& s1, const calc_string& s2);


template<class T> T from_string(const std::string& s);

/** Convert a string to a calc_string.
 *
 */
template<> inline
calc_string from_string<calc_string>(const std::string& s) {
	return calc_string(s);
}

#include "utility/stl_helper_functions.h"
#include "math/type_conversion.h"


/** Specialize conversion so as to create exceptions when necessary */
template<> inline double convert_element<double,calc_string>(const calc_string& x) {
	throw std::runtime_error("calc_string can not be converted to this type");
	return 0;
}
;

template<> inline std::string convert_element<std::string,calc_string>(const calc_string& x) {
	return x.get_my_string();
}
;

#endif /* CALC_STRING_OPERATORS_H_ */
