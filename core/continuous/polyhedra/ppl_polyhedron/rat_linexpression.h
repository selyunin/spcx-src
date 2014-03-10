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
#ifndef GUARD_rat_linexpression_h
#define GUARD_rat_linexpression_h

//#include <stdio.h>
#include <iostream>
#include <utility>

#include "core/continuous/polyhedra/ppl_polyhedron/extended_ppl.h"
#include "math/scalar_types/rational.h"
//#include "../../utility/calc_string.h"
#include "math/type_conversion.h"

namespace ppl_polyhedron {

class Rational_Linear_Expression :
	public Parma_Polyhedra_Library::Linear_Expression {
public:
	typedef Parma_Polyhedra_Library::Linear_Expression LE_type;
	// Constructors
	Rational_Linear_Expression() :
		denominator(1), my_LE(Integer(0)) {
	}
	;
	explicit Rational_Linear_Expression(const LE_type& lexp) :
		denominator(1), my_LE(lexp) {
	}
	;
	Rational_Linear_Expression(const LE_type& lexp, const Integer& d) :
		denominator(d), my_LE(lexp) {
	}
	;
	explicit Rational_Linear_Expression(const Rational& r) :
		denominator(r.get_den()), my_LE(r.get_num()) {
	}
	;
	explicit Rational_Linear_Expression(const Integer& n) :
		denominator(1), my_LE(n) {
	}
	;
	explicit Rational_Linear_Expression(const int& n) :
		denominator(1), my_LE(Integer(n)) {
	}
	;
	explicit Rational_Linear_Expression(const double& n) :
		denominator(1), my_LE(Integer(0)) {
		Rational r(n);
		denominator=r.get_den();
		my_LE=LE_type(r.get_num());
	}
	;
//	explicit Rational_Linear_Expression(const calc_string& n) :
//		denominator(1), my_LE(Integer(0)) {
//		Rational r(n);
//		denominator=r.get_den();
//		my_LE=LE_type(r.get_num());
//	}
//	;
	//Rational_Linear_Expression(const string& str);

	// Copy Constructor
	Rational_Linear_Expression(const Rational_Linear_Expression& r) :
		denominator(r.denominator), my_LE(r.get_LE()) {
		//this->Parma_Polyhedra_Library::Linear_Expression::operator=( (Parma_Polyhedra_Library::Linear_Expression)r );
	}
	;

	Rational rat_coefficient(const Parma_Polyhedra_Library::Variable j) const {
		return Rational(my_LE.coefficient(j), denominator);
	}
	;

	Rational rat_inhomogeneous_term() const {
		return Rational(my_LE.inhomogeneous_term(), denominator);
	}
	;

	void set_LE(const Parma_Polyhedra_Library::Linear_Expression& lexp) {
		my_LE=lexp;
	}
	;
	const LE_type& get_LE() const {
		return my_LE;
	}
	;
	const Integer& get_denominator() const {
		return denominator;
	}
	;
	Integer& access_denominator() {
		return denominator;
	}
	;

private:
	// Properties
	Integer denominator;

private:
	LE_type my_LE;
};

bool attempt_multiply(Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2);
bool attempt_division(Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2);

// Operators
Rational_Linear_Expression operator+(const Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2);
Rational_Linear_Expression operator-(const Rational_Linear_Expression& l1);
Rational_Linear_Expression operator-(const Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2);
Rational_Linear_Expression& operator+=(Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2);
Rational_Linear_Expression& operator-=(Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2);
Rational_Linear_Expression operator*(const Rational& c,
		const Rational_Linear_Expression& l);
Rational_Linear_Expression operator*(const Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2);
Rational_Linear_Expression operator/(const Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2);
Rational_Linear_Expression pow(const Rational_Linear_Expression& l1,
		const Rational_Linear_Expression& l2);


//Rational_Linear_Expression sqrt(const Rational_Linear_Expression& l1);

std::ostream& operator<<(std::ostream& os, const Rational_Linear_Expression &r);

typedef std::pair <Parma_Polyhedra_Library::Constraint, std::pair <Rational,Rational> >
		ConstraintRatPair;

}

/** Specialize conversion */
template<> inline ppl_polyhedron::Rational_Linear_Expression convert_element<ppl_polyhedron::Rational_Linear_Expression,Rational>(const Rational& x) {
	return ppl_polyhedron::Rational_Linear_Expression(x);
}
;


#endif
