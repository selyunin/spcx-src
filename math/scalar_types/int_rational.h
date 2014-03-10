#ifndef INT_RATIONAL_H_
#define INT_RATIONAL_H_

#include <boost/shared_ptr.hpp>
#include <boost/math/common_factor_rt.hpp>
#include "math/scalar_types/rational.h"

typedef boost::shared_ptr<Rational> rational_ptr;
typedef boost::shared_ptr<const Rational> rational_const_ptr;

/**
 * Representation of unbounded rationals. Use two int (numerator and denominator)
 * if they are lower than INT_MAX, and a Rational else.
 */
class int_rational {
public:
	typedef int fast_integer_type;
	typedef Rational slow_rational_type;
	int_rational() :
		num(0), den(1), my_rational(rational_ptr()), use_rational(false) { assert(!my_rational);
		assert(!use_rational);
	}
	;

	explicit int_rational(const rational_ptr& p) {
		my_rational = p;
		use_rational = true;
	}
	;
	explicit int_rational(const slow_rational_type& r) {
		my_rational = rational_ptr(new slow_rational_type(r));
		use_rational = true;
	}
	;
	explicit int_rational(const fast_integer_type& n, const fast_integer_type& d = 1) : my_rational(rational_ptr()) {
			if (d < 0) {
				num = -n;
				den = -d;
			} else if (d > 0) {
				num = n;
				den = d;
			} else
				throw std::runtime_error("division by zero");
			use_rational = false;
			assert(!my_rational);
			assert(!use_rational);
		}
		;
	int_rational(const int_rational& int_r) {
			if (int_r.get_use_rational()) {
				// @todo copy the rational to get a new pointer
				my_rational = rational_ptr(new slow_rational_type(*int_r.get_rational()));
				use_rational = true;
			} else {
				num = int_r.get_num_if_exist();
				den = int_r.get_den_if_exist();
				use_rational = false;
				my_rational=rational_ptr();
			}
		}
		;
	
	virtual ~int_rational(){}
	
	const slow_rational_type get_rational_num() const {
					if(get_use_rational())
						return slow_rational_type(my_rational->get_num());
					else
						return slow_rational_type(num);
				}
				;
	
	const slow_rational_type get_rational_den() const {
					if(get_use_rational())
						return slow_rational_type(my_rational->get_den());
					else
						return slow_rational_type(den);
				}
				;
	
	const rational_const_ptr get_rational() const {
		if(use_rational)
			return my_rational;
		else
			return rational_ptr(new slow_rational_type(get_num_if_exist(), get_den_if_exist()));
	}
	;

	const bool get_use_rational() const {
		return use_rational;
	}
	;
	
	friend bool is_bigger(const double& d);
	friend int_rational operator-(const int_rational& r1, const int_rational& r2);
	friend int_rational operator-(const int_rational& r);
	friend int_rational operator+(const int_rational& r1, const int_rational& r2);
	friend int_rational operator*(const int_rational& r1, const int_rational& r2);
	friend int_rational operator/(const int_rational& r1, const int_rational& r2);
	
	friend bool is_bigger_than_zero(const int_rational& r);
	friend bool is_zero(const int_rational& r);
	friend inline std::ostream& operator<<(std::ostream& os, const int_rational& r);
private:
	
		const fast_integer_type& get_num_if_exist() const {
			if(!get_use_rational())
				return num;
			else
				throw std::runtime_error("Invalid numerator");
		}
		;
		const fast_integer_type& get_den_if_exist() const {
			if(!get_use_rational())
				return den;
			else
				throw std::runtime_error("Invalid denominator");
		}
		;
		
	fast_integer_type num;
	fast_integer_type den; // >0
	rational_ptr my_rational;
	bool use_rational;
	static const double safety_margin; 
	static const double max_double;
};

inline std::ostream& operator<<(std::ostream& os, const int_rational& r) {
	if(r.get_use_rational())
		os << *(r.get_rational());
	else
	{	
		os << r.get_num_if_exist();
		if (r.get_den_if_exist()!= 1)
			os << "/" << r.get_den_if_exist();
	}
	return os;
}

// Operators

bool operator<(const int_rational& r1, const int_rational& r2);
bool operator<=(const int_rational& r1, const int_rational& r2);
bool operator>(const int_rational& r1, const int_rational& r2);
bool operator>=(const int_rational& r1, const int_rational& r2);
bool operator==(const int_rational& r1, const int_rational& r2);
bool operator!=(const int_rational& r1, const int_rational& r2);

int_rational operator-(const int_rational& r);
int_rational operator/(const int_rational& r1, const int_rational& r2);
int_rational operator*(const int_rational& r1, const int_rational& r2);

int_rational operator-(const int_rational& r1, const int_rational& r2);
int_rational operator+(const int_rational& r1, const int_rational& r2);

void simplify(double& a, double& b);

#endif /*INT_RATIONAL_H_*/
