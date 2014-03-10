/*
 * global_types.h
 *
 *  Created on: Sep 23, 2009
 *      Author: frehse
 */

#ifndef GLOBAL_TYPES_H_
#define GLOBAL_TYPES_H_

#include <iostream>
#include <math.h>

/** Stream output of __float128 */
template<class charT, class traits>
std::basic_ostream<charT, traits>& operator<<(
		std::basic_ostream<charT, traits>& os, const __float128& f) {
			long double d(f);
			os << d;
			return os;
		};

/** abs for __float128 */
inline __float128 abs(const __float128& f) {
	if (f<__float128(0))
	return -f;
	else
	return f;
};

/** sqrt for __float128 */
inline __float128 sqrt(const __float128& f) {
	double d(f);
	double d_sqrt = sqrt(d);
	return __float128(d_sqrt);
}

//namespace math {
///** Stream output of __float128 */
//template<class charT, class traits>
//std::basic_ostream<charT, traits>& operator<<(
//		std::basic_ostream<charT, traits>& os, const __float128& f) {
//			long double d(f);
//			os << d;
//			return os;
//		};
//
//		///** abs for __float128 */
//		//inline __float128 abs(const __float128& f) {
//		//	if (f<__float128(0))
//		//		return -f;
//		//	else
//		//		return f;
//		//};
//	}

//#include "../math/rational.h"
//#include "../utility/calc_string.h"
#include "utility/stl_helper_functions.h"

template<> inline __float128 from_string<__float128 > (const std::string& s) {
	// @todo do proper implementation
	return __float128(from_string<long double> (s));
}

class Rational;
class calc_string;

namespace global_types {

/** These are the types considered in the project. */
enum coefficient_type {
	UNDEFINED_TYPE,
	STD_BOOL,
	STD_INT,
	STD_DOUBLE,
	STD_LONG_DOUBLE,
	STD_FLOAT128,
	GMP_RATIONAL,
	CALC_STR
};

template<coefficient_type t>
struct type_selector {
};

template<>
struct type_selector<STD_BOOL> {
	typedef bool type;
};

template<>
struct type_selector<STD_INT> {
	typedef int type;
};

template<>
struct type_selector<STD_DOUBLE> {
	typedef double type;
};

template<>
struct type_selector<STD_LONG_DOUBLE> {
	typedef long double type;
};

template<>
struct type_selector<STD_FLOAT128> {
	typedef __float128 type;
};

template<>
struct type_selector<GMP_RATIONAL> {
	typedef Rational type;
};

template<>
struct type_selector<CALC_STR> {
	typedef calc_string type;
};

template<typename T>
struct type_identifier {
	static const coefficient_type coeff;
	static const std::string name;
};

template<typename T> const coefficient_type type_identifier<T>::coeff =
		UNDEFINED_TYPE;

template<typename T> const std::string type_identifier<T>::name = std::string(
		"UNDEFINED_TYPE");

template<> const coefficient_type
		type_identifier<type_selector<STD_BOOL>::type>::coeff;
template<> const coefficient_type
		type_identifier<type_selector<STD_INT>::type>::coeff;
template<> const coefficient_type type_identifier<
		type_selector<STD_DOUBLE>::type>::coeff;
template<> const coefficient_type type_identifier<
		type_selector<STD_LONG_DOUBLE>::type>::coeff;
template<> const coefficient_type type_identifier<
		type_selector<STD_FLOAT128>::type>::coeff;
template<> const coefficient_type type_identifier<
		type_selector<GMP_RATIONAL>::type>::coeff;
template<> const coefficient_type
		type_identifier<type_selector<CALC_STR>::type>::coeff;

template<> const std::string
		type_identifier<type_selector<STD_BOOL>::type>::name;
template<> const std::string
		type_identifier<type_selector<STD_INT>::type>::name;
template<> const std::string
		type_identifier<type_selector<STD_DOUBLE>::type>::name;
template<> const std::string type_identifier<
		type_selector<STD_LONG_DOUBLE>::type>::name;
template<> const std::string
		type_identifier<type_selector<STD_FLOAT128>::type>::name;
template<> const std::string
		type_identifier<type_selector<GMP_RATIONAL>::type>::name;
template<> const std::string
		type_identifier<type_selector<CALC_STR>::type>::name;

template<typename base_type, template<typename > class T>
class template_type_factory {
public:
	static base_type* create(coefficient_type c) {
		switch (c) {
		case STD_BOOL:
			return new T<type_selector<STD_BOOL>::type> ();
		case STD_INT:
			return new T<type_selector<STD_INT>::type> ();
		case STD_DOUBLE:
			return new T<type_selector<STD_DOUBLE>::type> ();
		case STD_LONG_DOUBLE:
			return new T<type_selector<STD_LONG_DOUBLE>::type> ();
//		case STD_FLOAT128:
//			return new T<type_selector<STD_FLOAT128>::type> ();
		case GMP_RATIONAL:
			return new T<type_selector<GMP_RATIONAL>::type> ();
		case CALC_STR:
			return new T<type_selector<CALC_STR>::type> ();
		default:
			throw std::runtime_error(
					"template_type_factory : unsupported type " + int2string(c));
		}
	}
	;
	static base_type* create_number(coefficient_type c) {
		switch (c) {
		case STD_DOUBLE:
			return new T<type_selector<STD_DOUBLE>::type> ();
		case STD_LONG_DOUBLE:
			return new T<type_selector<STD_LONG_DOUBLE>::type> ();
//		case STD_FLOAT128:
//			return new T<type_selector<STD_FLOAT128>::type> ();
		case GMP_RATIONAL:
			return new T<type_selector<GMP_RATIONAL>::type> ();
		default:
			throw std::runtime_error(
					"template_type_factory : unsupported type " + int2string(c));
		}
	}
	;
	template<typename arg_type>
	static base_type* create_number(coefficient_type c,const arg_type& arg) {
		switch (c) {
		case STD_DOUBLE:
			return new T<type_selector<STD_DOUBLE>::type> (arg);
		case STD_LONG_DOUBLE:
			return new T<type_selector<STD_LONG_DOUBLE>::type> (arg);
//		case STD_FLOAT128:
//			return new T<type_selector<STD_FLOAT128>::type> ();
		case GMP_RATIONAL:
			return new T<type_selector<GMP_RATIONAL>::type> (arg);
		default:
			throw std::runtime_error(
					"template_type_factory : unsupported type " + int2string(c));
		}
	}
	;
};

/** Returns R=T<C>::implement(S), where C is a type determined by coeff, and R
 * and S are arbitrary types.
 * T must implement R implement(S).
 */
template<typename R, typename S, template<typename > class T>
class coefficient_type_caller {
public:
	static R call(S s, coefficient_type c) {
		switch (c) {
		case STD_DOUBLE:
			return T<type_selector<STD_DOUBLE>::type>::implement(s);
		case STD_LONG_DOUBLE:
			return T<type_selector<STD_LONG_DOUBLE>::type>::implement(s);
//		case STD_FLOAT128:
//			return T<type_selector<STD_FLOAT128>::type>::implement(s);
		case GMP_RATIONAL:
			return T<type_selector<GMP_RATIONAL>::type>::implement(s);
		default:
			throw std::runtime_error(
					"coefficient_type_caller : unsupported type " + int2string(
							c));
		}
	}
	;
};

const coefficient_type default_rational_type_id = GMP_RATIONAL;
const coefficient_type default_float_type_id = STD_DOUBLE;
const coefficient_type precise_float_type_id=STD_LONG_DOUBLE;
//const coefficient_type precise_float_type_id = STD_FLOAT128;

typedef type_selector<default_rational_type_id>::type rational_type;
typedef type_selector<default_float_type_id>::type float_type;
typedef type_selector<precise_float_type_id>::type precise_float_type;

/** Measure the decimal precision of a scalar type
 *
 * Returns -1 if the precision is 500 bits or more.
 */
template<typename scalar_type>
int measure_precision() {
	int i = 0;
	int max_i = 500;
	scalar_type one(1), eps, last_res;
    eps = scalar_type(1);
    last_res = one;
	while (one + eps / scalar_type(2) != last_res && i < max_i) {
		++i;
		last_res = one + eps;
		eps = eps / scalar_type(2);
	}
	// check for substraction in case rounding
	// messes things up
	int j = 0;
    eps = scalar_type(1);
    last_res = one;
	while (one - eps / scalar_type(2) != last_res && j < max_i) {
		++j;
		last_res = one - eps;
		eps = eps / scalar_type(2);
	}
	// take the min for addition and substraction
	i = std::min(i,j);

	if (i >= max_i)
		i = -1;
	return i;
}

/** Returns a description of the floating point number types used. */
std::string coefficient_type_summary();

}

#endif /* GLOBAL_TYPES_H_ */
