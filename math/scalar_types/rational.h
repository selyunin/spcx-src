#ifndef RATIONAL_H_
#define RATIONAL_H_

#include <iostream>
#include <stdexcept>
//#include <cmath>
#include <math.h>
#include "gmp.h"
#include "gmpxx.h"
#include "utility/stl_helper_functions.h"
#include "math/type_conversion.h"

typedef mpz_class Integer;

class Rational;

/** Representation of unbounded rationals. Canonical form in which denominator is always >0. */
class Rational {
public:
	Rational() :
		my_q(Integer(0), Integer(1)) {
		my_q.canonicalize();
	}
	;
	explicit Rational(const mpq_class& q) :
		my_q(q) {
	}
	;
	Rational(const Integer& n, const Integer& d) :
		my_q(n, d) {
		// if (d<0) my_q=mpq_class(-n,-d); if (d==0) throw ...
		// From the man pages: "In all the following constructors, if a fraction is given then it should be in canonical form, or if not then `mpq_class::canonicalize' called."
		my_q.canonicalize();
	}
	;
	explicit Rational(const Integer& n) :
		my_q(n, Integer(1)) {
		my_q.canonicalize();
	}
	;
	explicit Rational(const int& n) :
		my_q(Integer(n), Integer(1)) {
		my_q.canonicalize();
	}
	;
	explicit Rational(const double& d) :
		my_q(d) {
		my_q.canonicalize();
	}
	;
	explicit Rational(const long double& d) :
		// @todo do a proper conversion
		my_q((double)d) {
		my_q.canonicalize();
	}
	;
	explicit Rational(const __float128& d) :
		// @todo do a proper conversion
		my_q((double)d) {
		my_q.canonicalize();
	}
	;
	//	explicit Rational(const calc_string& c) {
	//		init_with_string(c.get_my_string());
	//	}
	//	;
	explicit Rational(const std::string& str) {
		init_with_string(str);
	}
	;

	virtual ~Rational();

	double get_double() const {
		return my_q.get_d();
	}
	;
	/** double typecast overloaded
	 * @todo this seems highly suspicious ans may lead to implicit casts!
	 */
	/*
	 operator double() const {
	 return get_double();
	 }
	 */
	Integer get_Integer() const {
		return (my_q.get_num() / my_q.get_den());
	}
	;
	int get_int() const {
		//Integer a(get_num()/get_den());
		//mpz_t z(get_num()/get_den());
		//return mpz_get_si(z);
		return (int) get_double();
	}
	;

	const Integer& get_num() const {
		return my_q.get_num();
	}
	;
	const Integer& get_den() const {
		return my_q.get_den();
	}
	;
	const mpq_class& get_mpq_class() const {
		return my_q;
	}
	/**
	 * \brief Computes absolute value of the rational
	 */
	const Rational abs() const {
		if (get_num() < 0 && get_den() >= 0)
			return Rational(-1 * get_num(), get_den());
		else if (get_num() >= 0 && get_den() < 0)
			return Rational(get_num(), -1 * get_den());
		else
			return Rational(get_num(), get_den());
	}
	;
	Rational& operator+=(const Rational& r2);
	Rational& operator-=(const Rational& r2);
	Rational& operator*=(const Rational& r2);
	Rational& operator/=(const Rational& r2);
	void swap(Rational& r2);

private:
	mpq_class my_q;
	void init_with_string(const std::string& str);
};

// Operators

inline bool operator<(const Rational& r1, const Rational& r2) {
	return r1.get_mpq_class() < r2.get_mpq_class();
}
;
inline bool operator<=(const Rational& r1, const Rational& r2) {
	return r1.get_mpq_class() <= r2.get_mpq_class();
}
;
inline bool operator>(const Rational& r1, const Rational& r2) {
	return r1.get_mpq_class() > r2.get_mpq_class();
}
;
inline bool operator>=(const Rational& r1, const Rational& r2) {
	return r1.get_mpq_class() >= r2.get_mpq_class();
}
;
inline bool operator==(const Rational& r1, const Rational& r2) {
	return r1.get_mpq_class() == r2.get_mpq_class();
}
;
inline bool operator!=(const Rational& r1, const Rational& r2) {
	return r1.get_mpq_class() != r2.get_mpq_class();
}
;

Rational operator-(const Rational& r);
Rational operator/(const Rational& r1, const Rational& r2);
Rational operator*(const Rational& r1, const Rational& r2);

Rational operator-(const Rational& r1, const Rational& r2);
Rational operator+(const Rational& r1, const Rational& r2);

template<> inline Rational from_string<Rational> (const std::string& s) {
	return Rational(s);
}

inline std::ostream& operator<<(std::ostream& os, const Rational& r) {
	os << r.get_num();
	if (r.get_den() != Integer(1)) {
		os << "/" << r.get_den();
	}
	return os;
}

/** Specialize conversion */

template<> inline double convert_element<double, Rational> (const Rational& x) {
	return x.get_double();
}
;

template<> inline long double convert_element<long double, Rational> (const Rational& x) {
	long double num = x.get_num().get_d();
	long double den = x.get_den().get_d();
	return num / den;
}
;

template<> inline Rational convert_element<Rational, long double> (const long double& x) {
	return Rational(x);
}
;

template<> inline int convert_element<int, Rational> (const Rational& x) {
	return x.get_int();
}
;

template<>
class converter<Rational, Integer> {
public:
	static Rational convert(const Integer& x) {
		return Rational(x);
	}
	;
};

template<>
class converter<Rational, int> {
public:
	static Rational convert(const int& x) {
		return Rational(x);
	}
	;
};

template<>
class converter<Rational, double> {
public:
	static Rational convert(const double& x) {
		return Rational(x);
	}
	;
};

/** Convert __float128 to Rational */
template<>
class converter<Rational, __float128> {
public:
static Rational convert(const __float128& x) {
	Rational y(x);
	return y;
};
};

/** Convert Rational to __float128 */
template<>
class converter<__float128, Rational> {
public:
static double convert(const Rational& x) {
	__float128 y(x.get_double());
	return y;
};
};
template<>
class converter<Rational, std::string> {
public:
	static Rational convert(const std::string& x) {
		return Rational(x);
	}
	;
};

inline Rational abs(const Rational& r) {
	return r.abs();
}
;

using ::sqrt;
inline Rational sqrt(const Rational& r) {
	return Rational(sqrt(r.get_double()));
}
;


inline Rational log(const Rational& r) {
	return Rational(log(r.get_double()));
}
;

inline Rational exp(const Rational& r) {
	return Rational(exp(r.get_double()));
}
;

inline Rational sin(const Rational& r) {
	return Rational(sin(r.get_double()));
}
;

inline Rational tan(const Rational& r) {
	return Rational(tan(r.get_double()));
}
;

inline Rational cos(const Rational& r) {
	return Rational(cos(r.get_double()));
}
;

inline Rational pow(const Rational& r, const Rational& e) {
	return Rational(pow(r.get_double(),e.get_double()));
}
;

#endif /*RATIONAL_H_*/

