#include "math/scalar_types/rational.h"

using namespace std;

Rational::~Rational()
{
}

void Rational::init_with_string(const std::string& str) {
	string flt(""), expon(""), numstr(""), denstr("");

	// if there's a fraction,  put the numerator in flt and the denominator in d
	if (str.find("/")!=string::npos){
		string num=string_before(str, "/");
		string den=string_after(str, "/");
		*this=Rational(num)/Rational(den);
	}

	// if there's an exponent, find it and put the rest in flt
	int e=0;
	if (str.find("e")!=string::npos) {
		flt=string_before(str, "e");
		expon=string_after(str, "e");
		e=atoi(expon.c_str());
	} else if (str.find("E")!=string::npos) {
		flt=string_before(str, "E");
		expon=string_after(str, "E");
		e=atoi(expon.c_str());
	}
	else
		flt=str;

	numstr=string_before(flt, ".");
	denstr=string_after(flt, ".");

	// Conversion is simple:
	// Count the number of decimal places
	// and subtract them from the exponent
	e-=denstr.size();

	//cout << "str:" << str << "; n:" <<numstr <<";" << " d:"<< denstr <<  " Net Exponent:" << e << endl << flush;

	mpq_class q(numstr+denstr, 10); // 10 is the base
	my_q.canonicalize();

	//cout << endl << numstr+denstr << "->" << q << endl << flush;

	mpz_class den(1);
	// for a positive exponent, multiply the numerator
	while (e>0) {
		--e;
		q.get_num()*=mpz_class(10);
	};
	// for a negative exponent, multiply the denominator
	while (e<0) {
		++e;
		den*=mpz_class(10);
	};


	q.get_den()*=den;

	//cout << endl << str << "-e>" << q << endl << flush;

	my_q=q;
	my_q.canonicalize();

	//cout << endl << str << "==>" << *this << endl << flush;
}

// bool operator<(const Rational& r1, const Rational& r2)
// {
// 	// Note: only true if denominator is always greater than zero
// 	if (((r2.get_den() > 0) && (r1.get_den() > 0)) || ((r2.get_den() < 0) && (r1.get_den() < 0)))
// 		return (r2.get_den())*(r1.get_num())<(r1.get_den())*(r2.get_num());
// 	else if (((r2.get_den() < 0) && (r1.get_den() > 0)) || ((r2.get_den() > 0) && (r1.get_den() < 0)))
// 		return (r2.get_den())*(r1.get_num())>(r1.get_den())*(r2.get_num());
// 	else
// 	{
// 		throw_error("Rational operator<(const Rational& r1, const Rational& r2): denominator == 0!");
// 		return false;
// 	};
// }
//
// bool operator<=(const Rational& r1, const Rational& r2)
// {
// 	// Note: only true if denominator is always greater than zero
// 	if (((r2.get_den() > 0) && (r1.get_den() > 0)) || ((r2.get_den() < 0) && (r1.get_den() < 0)))
// 		return (r2.get_den())*(r1.get_num())<=(r1.get_den())*(r2.get_num());
// 	else if (((r2.get_den() < 0) && (r1.get_den() > 0)) || ((r2.get_den() > 0) && (r1.get_den() < 0)))
// 		return (r2.get_den())*(r1.get_num())>=(r1.get_den())*(r2.get_num());
// 	else
// 	{
// 		throw_error("Rational operator<(const Rational& r1, const Rational& r2): denominator == 0!");
// 		return false;
// 	};
// }
//
// bool operator>(const Rational& r1, const Rational& r2)
// {
// 	// Note: only true if denominator is always greater than zero
// 	if (((r2.get_den() > 0) && (r1.get_den() > 0)) || ((r2.get_den() < 0) && (r1.get_den() < 0)))
// 		return (r2.get_den())*(r1.get_num())>(r1.get_den())*(r2.get_num());
// 	else if (((r2.get_den() < 0) && (r1.get_den() > 0)) || ((r2.get_den() > 0) && (r1.get_den() < 0)))
// 		return (r2.get_den())*(r1.get_num())<(r1.get_den())*(r2.get_num());
// 	else
// 	{
// 		throw_error("Rational operator<(const Rational& r1, const Rational& r2): denominator == 0!");
// 		return false;
// 	};
// }
//
// bool operator>=(const Rational& r1, const Rational& r2)
// {
// 	// Note: only true if denominator is always greater than zero
// 	if (((r2.get_den() > 0) && (r1.get_den() > 0)) || ((r2.get_den() < 0) && (r1.get_den() < 0)))
// 		return (r2.get_den())*(r1.get_num())>=(r1.get_den())*(r2.get_num());
// 	else if (((r2.get_den() < 0) && (r1.get_den() > 0)) || ((r2.get_den() > 0) && (r1.get_den() < 0)))
// 		return (r2.get_den())*(r1.get_num())<=(r1.get_den())*(r2.get_num());
// 	else
// 	{
// 		throw_error("Rational operator<(const Rational& r1, const Rational& r2): denominator == 0!");
// 		return false;
// 	};
// }

Rational operator-(const Rational& r) {
	return Rational(-r.get_num(), r.get_den());
}

Rational operator/(const Rational& r1, const Rational& r2) {
	if (r2.get_num() == Integer(0) )
		throw std::runtime_error("Division of rational by zero!");
	return Rational(r1.get_num()*r2.get_den(), r2.get_num()*r1.get_den());
}

Rational operator*(const Rational& r1, const Rational& r2) {
	return Rational(r1.get_num()*r2.get_num(), r2.get_den()*r1.get_den());
}

Rational operator-(const Rational& r1, const Rational& r2) {
	return Rational(r1.get_num()*r2.get_den()-r1.get_den()*r2.get_num(),
			r2.get_den()*r1.get_den());
}

Rational operator+(const Rational& r1, const Rational& r2) {
	return Rational(r1.get_num()*r2.get_den()+r1.get_den()*r2.get_num(),
			r2.get_den()*r1.get_den());
}

Rational& Rational::operator+=(const Rational& r2) {
	/*
	my_q*=r2.get_den();
	my_q+=get_den()*r2.get_num();
	my_q/=r2.get_den();
	*/
	my_q+=r2.my_q;
	return *this;
}

Rational& Rational::operator-=(const Rational& r2) {
	my_q-=r2.my_q;
	return *this;
}

Rational& Rational::operator*=(const Rational& r2) {
	my_q*=r2.my_q;
	return *this;
}

Rational& Rational::operator/=(const Rational& r2) {
	my_q/=r2.my_q;
	return *this;
}

void Rational::swap(Rational &r1){
	mpq_class temp;
	temp = my_q;
	my_q = r1.my_q;
	r1.my_q = temp;
}
