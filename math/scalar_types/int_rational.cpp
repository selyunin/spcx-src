#include "math/scalar_types/int_rational.h"
#include <climits> /* for INT_MAX -wsc */

const double int_rational::safety_margin = 0.95;
const double int_rational::max_double=safety_margin*double(INT_MAX);


bool is_bigger(const double& d){
	return (fabs(d) > int_rational::max_double);
}

//boolean functions to simplify comparison operators
bool is_bigger_than_zero(const int_rational& r){
	if(r.get_use_rational())
			return (*(r.get_rational()) > Rational(0));
	else
			return (r.get_num_if_exist() > 0);
}

bool is_zero(const int_rational& r){
	if(r.get_use_rational())
			return (*(r.get_rational()) == Rational(0));
	else
			return (r.get_num_if_exist() == 0);
}

bool is_bigger_or_equal_than_zero(const int_rational& r){
	return (is_bigger_than_zero(r) || is_zero(r));
}

int_rational operator+(const int_rational& r1, const int_rational& r2) {
	if (r1.get_use_rational() || r2.get_use_rational())
			return int_rational(*(r1.get_rational()) + *(r2.get_rational()));
	else {
		double n1 = double(r1.get_num_if_exist())*double(r2.get_den_if_exist());
		double n2 = double(r2.get_num_if_exist())*double(r1.get_den_if_exist());
		double d = double(r1.get_den_if_exist())*double(r2.get_den_if_exist());
		if (is_bigger(n1 * 2)|| is_bigger(n2 * 2)|| is_bigger(d)) {
			int_rational::slow_rational_type num = r1.get_rational_num()*r2.get_rational_den() + r2.get_rational_num()*r1.get_rational_den();
			int_rational::slow_rational_type den = r1.get_rational_den()*r2.get_rational_den();
			return int_rational(num/den);
		} else {
			int_rational::fast_integer_type num = r1.get_num_if_exist()*r2.get_den_if_exist() + r2.get_num_if_exist()*r1.get_den_if_exist();
			int_rational::fast_integer_type den = r1.get_den_if_exist()*r2.get_den_if_exist();
			
			int_rational::fast_integer_type div = boost::math::gcd(num, den);
			den /= div;
			num /= div;
			return int_rational(num, den);
		}
	}
}

int_rational operator-(const int_rational& r1, const int_rational& r2) {
	if (r1.get_use_rational() || r2.get_use_rational())
			return int_rational(*(r1.get_rational()) - *(r2.get_rational()));
	else {
		double n1 = double(r1.get_num_if_exist())*double(r2.get_den_if_exist());
		double n2 = double(r2.get_num_if_exist())*double(r1.get_den_if_exist());
		double d = double(r1.get_den_if_exist())*double(r2.get_den_if_exist());
		if (is_bigger(n1 * 2)|| is_bigger(n2 * 2)|| is_bigger(d)) {
			int_rational::slow_rational_type num = r1.get_rational_num()*r2.get_rational_den() - r2.get_rational_num()*r1.get_rational_den();
						int_rational::slow_rational_type den = r1.get_rational_den()*r2.get_rational_den();
			return int_rational(num/den);
		} else {
			int_rational::fast_integer_type num = r1.get_num_if_exist()*r2.get_den_if_exist() - r2.get_num_if_exist()*r1.get_den_if_exist();
			int_rational::fast_integer_type den = r1.get_den_if_exist()*r2.get_den_if_exist();
			boost::math::gcd_evaluator<int_rational::fast_integer_type> gcd;
			int_rational::fast_integer_type div = gcd(num, den);
			den = den/div;
			num = num/div;
			return int_rational(num, den);
		}
	}
}

int_rational operator*(const int_rational& r1, const int_rational& r2) {
	if (r1.get_use_rational() || r2.get_use_rational())
			return int_rational((*(r1.get_rational())) * (*(r2.get_rational())));
	else {
		double c = double(r1.get_num_if_exist())*double(r2.get_num_if_exist());
		double d = double(r1.get_den_if_exist())*double(r2.get_den_if_exist());
		if (is_bigger(c) || is_bigger(d)) {
			int_rational::slow_rational_type num = r1.get_rational_num()*r2.get_rational_num();
			int_rational::slow_rational_type den = r1.get_rational_den()*r2.get_rational_den();
			return int_rational(num/den);
		} else {
			int_rational::fast_integer_type num = r1.get_num_if_exist()*r2.get_num_if_exist();
			int_rational::fast_integer_type den = r1.get_den_if_exist()*r2.get_den_if_exist();
			boost::math::gcd_evaluator<int_rational::fast_integer_type> gcd;
			int_rational::fast_integer_type div = gcd(num, den);
			den = den/div;
			num = num/div;
			return int_rational(num, den);
		}
	}
}

int_rational operator/(const int_rational& r1, const int_rational& r2) {
	int_rational r;
	if (r1.get_use_rational() || r2.get_use_rational())
			r = int_rational((*(r1.get_rational())) / (*(r2.get_rational())));
	else {
		double c = double(r1.get_num_if_exist())*double(r2.get_den_if_exist());
		double d = double(r1.get_den_if_exist())*double(r2.get_num_if_exist());
		if (is_bigger(c) || is_bigger(d)) {
			int_rational::slow_rational_type num = r1.get_rational_num()*r2.get_rational_den();
			int_rational::slow_rational_type den = r1.get_rational_den()*r2.get_rational_num();
			r = int_rational(num/den);
		} else {
			int_rational::fast_integer_type num = r1.get_num_if_exist()*r2.get_den_if_exist();
			int_rational::fast_integer_type den = r1.get_den_if_exist()*r2.get_num_if_exist();
			boost::math::gcd_evaluator<int_rational::fast_integer_type> gcd;
			int_rational::fast_integer_type div = gcd(num, den);
			den = den/div;
			num = num/div;
			r = int_rational(num, den);
		};
	}
	return r;
}

int_rational operator-(const int_rational& r) {
	if (r.get_use_rational())
		return int_rational( - *(r.get_rational()));
	else
		return int_rational(-r.get_num_if_exist(), r.get_den_if_exist());
}

// Operators

bool operator<(const int_rational& r1, const int_rational& r2) {
	return(is_bigger_than_zero(r2-r1));
}

bool operator<=(const int_rational& r1, const int_rational& r2) {
	return(is_bigger_or_equal_than_zero(r2-r1));
}

bool operator>=(const int_rational& r1, const int_rational& r2) {
	return(is_bigger_or_equal_than_zero(r1-r2));
}

bool operator>(const int_rational& r1, const int_rational& r2) {
	return(is_bigger_than_zero(r1-r2));
}

bool operator==(const int_rational& r1, const int_rational& r2) {
	return(is_zero(r2-r1));
}

bool operator!=(const int_rational& r1, const int_rational& r2) {
	return(!is_zero(r2-r1));
}
