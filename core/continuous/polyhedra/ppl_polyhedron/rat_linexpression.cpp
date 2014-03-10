/***************************************************************************
 *   Copyright (C) 2004 by Goran Frehse                                    *
 *   gfrehse@localhost                                                     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "core/continuous/polyhedra/ppl_polyhedron/rat_linexpression.h"

#include <vector>
#include <list>
#include <set>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <string>
//#include <stdlib.h>

//#include <ppl.hh>

#include "core/continuous/polyhedra/ppl_polyhedron/general.h"

namespace ppl_polyhedron {

using namespace std;
using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;

//Rational_Linear_Expression::Rational_Linear_Expression(const string& str) :
//	my_LE(0), get_denominator()(1) {
//	//	mpf_class f(str);
//	//	Rational r(f);
//	Rational r(str);
//	//	cout << str << "r: " << r << endl;
//	//	cout << r.get_num() << flush << endl;
//	//	cout << r.get_d() << endl;
//	//this->Linear_Expression::operator=(Linear_Expression(r.get_num()) );
//	my_LE=r.get_num();
//	get_denominator()=r.get_den();
//}

Rational_Linear_Expression operator+(const Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2) {
	Rational_Linear_Expression r;
	if (l2.get_denominator() != l1.get_denominator()) {
		r=Rational_Linear_Expression(l2.get_denominator()*l1.get_LE() + l1.get_denominator() *l2.get_LE(), l1.get_denominator()*l2.get_denominator());
	} else {
		r=Rational_Linear_Expression(l1.get_LE() + l2.get_LE(),l1.get_denominator());
	}

	return r;
}

Rational_Linear_Expression operator-(const Rational_Linear_Expression& l1) {
	Rational_Linear_Expression r;
	r=Rational_Linear_Expression(-l1.get_LE(), l1.get_denominator());

	return r;
}

Rational_Linear_Expression operator-(const Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2) {
	Rational_Linear_Expression r;
	if (l2.get_denominator() != l1.get_denominator()) {
		r=Rational_Linear_Expression(l2.get_denominator()*l1.get_LE()
				- l1.get_denominator() *l2.get_LE(), l1.get_denominator()*l2.get_denominator());
	} else {
		r=Rational_Linear_Expression(l1.get_LE() - l2.get_LE(), l1.get_denominator());
	}

	return r;
}

Rational_Linear_Expression& operator+=(Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2) {
	if (l2.get_denominator() != l1.get_denominator()) {
		l1.set_LE(l2.get_denominator()*l1.get_LE() + l1.get_denominator() *l2.get_LE());
		l1.access_denominator()*=l2.get_denominator();
	} else {
		l1.set_LE(l1.get_LE() + l2.get_LE());
	}

	return l1;
}

Rational_Linear_Expression& operator-=(Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2) {
	if (l2.get_denominator() != l1.get_denominator()) {
		l1.set_LE(l2.get_denominator()*l1.get_LE() - l1.get_denominator() *l2.get_LE());
		l1.access_denominator()*=l2.get_denominator();
	} else {
		l1.set_LE(l1.get_LE() - l2.get_LE());
	}

	return l1;
}

Rational_Linear_Expression operator*(const Rational& c,
		const Rational_Linear_Expression& l) {
	Rational_Linear_Expression r(c.get_num()*l.get_LE(),l.get_denominator()*c.get_den());
	return r;
}

bool attempt_multiply(Rational_Linear_Expression& l1, const
Rational_Linear_Expression& l2) { // multiply if l1 or l2 is a scalar and put the result in l1. Return true if successful.

	if (get_true_dimension(l1.get_LE())==0) // l1 is a scalar
	{
		l1.set_LE(l1.get_LE().inhomogeneous_term()*l2.get_LE());
		l1.access_denominator()*=l2.get_denominator();
		return true;
	} else if (get_true_dimension(l2.get_LE())==0) {
		l1.set_LE(l1.get_LE()*l2.get_LE().inhomogeneous_term());
		l1.access_denominator()*=l2.get_denominator();
		return true;
	} else
		return false;
}

bool attempt_division(Rational_Linear_Expression& l1, const
Rational_Linear_Expression& l2) { // divide if l2 is a scalar and put the result in l1. Return true if successful.

	if (l2.get_LE().inhomogeneous_term() == Integer(0))
		throw std::runtime_error("Division of rational by zero!");

	if (get_true_dimension(l2.get_LE())==0) {
		l1.set_LE(l2.get_denominator()*l1.get_LE());
		l1.access_denominator()*=l2.get_LE().inhomogeneous_term();
		return true;
	} else
		return false;
}

Rational_Linear_Expression operator*(const Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2) {
	Rational_Linear_Expression r(l1);
	bool success=attempt_multiply(r, l2);
	if (!success)
		throw std::runtime_error("Attempt to multiply nontrivial linear expressions.");
	return r;
}

Rational_Linear_Expression operator/(const Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2) {
	Rational_Linear_Expression r(l1);
	bool success=attempt_division(r, l2);
	if (!success)
		throw std::runtime_error("Attempt to divide nontrivial linear expressions.");
	return r;
}

Rational_Linear_Expression pow(const Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2) {

	if (get_true_dimension(l1.get_LE())==0 && get_true_dimension(l2.get_LE())==0) {
		Rational p = pow( l1.rat_inhomogeneous_term() , l2.rat_inhomogeneous_term() );
		return Rational_Linear_Expression(p);
	} else {
		throw std::runtime_error("Attempt to exponentiate nontrivial linear expressions.");
	}
	return Rational_Linear_Expression();
}

/*Rational_Linear_Expression sqrt(const Rational_Linear_Expression& l1)
{
	return l1;
}*/


std::ostream& operator<<(std::ostream& os, const Rational_Linear_Expression &r) {
	using Parma_Polyhedra_Library::IO_Operators::operator<<;
	os << "(" << r.get_LE() << ")/" << r.get_denominator();
	return os;
}

}
