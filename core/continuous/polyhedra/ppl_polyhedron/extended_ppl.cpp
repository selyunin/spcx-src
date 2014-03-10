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

#include "core/continuous/polyhedra/ppl_polyhedron/extended_ppl.h"


//#include <stdio.h>
#include <iostream>
#include <vector>
//#include <list>
#include <set>
#include <map>
#include <algorithm>
#include <stdexcept>
#include "math/scalar_types/rational.h"

namespace ppl_polyhedron {

using namespace std;
using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;
using Parma_Polyhedra_Library::IO_Operators::operator<<;

Linear_Expression Linear_Expression_subset(const Linear_Expression& o,
		const clock_ref_set& crs) {
	Linear_Expression e;
	for (clock_ref_set::const_iterator it=crs.begin(); it!=crs.end(); ++it)
		e+=o.coefficient(Variable(*it))*Variable(*it);
	return e;
}

void print_Linear_Expression_hom(ostream& s, const Linear_Expression& c,
		const variable_id_map& vnvec, bool& first) {
	// does not print the inhomogeneous term
	// returns true if nothing was output
	const int num_variables = c.space_dimension();
	for (int v = 0; v < num_variables; ++v) {
		Integer cv = c.coefficient(Variable(v));
		if (cv != 0) { // put non-primed on the other side
			if (!first) {
				if (cv > 0) {
					s << " + ";
				} else {
					s << " - ";
					cv=-cv;
				}
			} else {
				first = false;
				if (cv == -1) {
					s << "- ";
					cv=-cv;
				}
			}
			if (cv != 1) {
				s << cv << "*";
			}
			if (v<(int)vnvec.size()) {
				s << vnvec.get_name(v);
			} else if (v<(int)(2*vnvec.size())) {
				s << vnvec.get_name(v-vnvec.size()) << "'";
			} else
				throw_error("print_constraint:: variable_id_map not long enough");
		}
	}
}

void print_relation_symbol(ostream& s,
		Parma_Polyhedra_Library::Constraint::Type t, bool reverse_sign) {
	const char* relation_symbol = 0;
	switch (t) {
	case Constraint::EQUALITY:
		relation_symbol = "==";
		break;
	case Constraint::NONSTRICT_INEQUALITY:
		if (reverse_sign)
			relation_symbol = "<=";
		else
			relation_symbol = ">=";
		break;
	case Constraint::STRICT_INEQUALITY:
		if (reverse_sign)
			relation_symbol = "<";
		else
			relation_symbol = ">";
		break;
	}
	s << relation_symbol;
}

void print_inhomogeneous_term(ostream& s, const Constraint& c,
		bool reverse_sign, bool& first) {
	// reverse_sign reverses sign,
	// first determines if a + sign is inserted before a positive term
	// if first, there is always an output, even when the coeff is zero
	Integer v=c.inhomogeneous_term();
	if (reverse_sign)
		v=-v;

	if (!first && v>0) {
		s << " + ";
	}
	if (!first && v<0) {
		s << " - ";
		v=-v;
	}
	if (first || v!=0) {
		s << v;
	}
	first=false;
}

void print_constraint(ostream& s, const Constraint& c,
		const variable_id_map& vnvec) {
	// Special Cases:
	if (c.is_tautological()) {
		s << "true";
	} else if (c.is_inconsistent()) {
		s << "false";
	} else {
		bool first=true;
		const int num_variables = c.space_dimension();
		clock_ref_set variables(vnvec.size());
		clock_ref_set primed_variables(vnvec.size(), num_variables);
		primed_variables.intersection_assign(get_nonzero_vars(c));
		clock_ref_set primed_neg=primed_variables;
		primed_neg.difference_assign(get_positive_vars(c));
		clock_ref_set primed_pos=primed_variables;
		primed_pos.intersection_assign(get_positive_vars(c));

		Linear_Expression e;

		bool reverse_sign=false;
		if (get_positive_vars(c).size()==0 || (primed_pos.size()==0 && primed_neg.size()>0)) {
			reverse_sign=true;
		}

		if (reverse_sign) {
			// negate linear_expression
			e=-ConstraintHom2Linear_Expression(&c);
		} else {
			e=ConstraintHom2Linear_Expression(&c);
		}

		// recompute pos and neg sets because of negation
		primed_neg=primed_variables;
		primed_neg.difference_assign(get_positive_vars(e));
		primed_pos=primed_variables;
		primed_pos.intersection_assign(get_positive_vars(e));

		clock_ref_set neg=variables;
		neg.difference_assign(get_positive_vars(e));
		clock_ref_set pos=variables;
		pos.intersection_assign(get_positive_vars(e));

		using namespace ppl_polyhedron;
		if (primed_variables.size()==0) // no primed variables
		{
			ppl_polyhedron::print_Linear_Expression_hom(s, ppl_polyhedron::Linear_Expression_subset(e, pos),
					vnvec, first);
			s << " ";
			ppl_polyhedron::print_relation_symbol(s, c.type(), reverse_sign);
			s << " ";
			first=true;
			ppl_polyhedron::print_Linear_Expression_hom(s, -ppl_polyhedron::Linear_Expression_subset(e, neg),
					vnvec, first);
		} else // if (primed_variables.size()!=0) // primed variables
		{
			ppl_polyhedron::print_Linear_Expression_hom(s, ppl_polyhedron::Linear_Expression_subset(e,
					primed_pos), vnvec, first);
			s << " ";
			ppl_polyhedron::print_relation_symbol(s, c.type(), reverse_sign);
			s << " ";
			first=true;
			ppl_polyhedron::print_Linear_Expression_hom(s, -ppl_polyhedron::Linear_Expression_subset(e,
					primed_neg), vnvec, first);
			ppl_polyhedron::print_Linear_Expression_hom(s, -ppl_polyhedron::Linear_Expression_subset(e,
					variables), vnvec, first);
		}

		ppl_polyhedron::print_inhomogeneous_term(s, c, !reverse_sign, first);
	}
}

void Constraint_List::print_phaver(ostream& s, const variable_id_map& vnvec) const {
	if (begin()==end())
		s << "true";

	// filter opposing inequalities and output them in one line
	Constraint_List cl=*this;
	Constraint_List::iterator it=cl.begin(), it2;
	clock_ref_set crs, crs2;
	bool found;
	Linear_Expression e, e2;
	bool first(true);
	bool first_constraint(true);
	while (it!=cl.end()) {
		found=false;
		crs=get_nonzero_vars(*it);
		if (crs.size()>0) {
			e=ConstraintHom2Linear_Expression(&(*it));
			// find matching one with negative coeff
			it2 = it;
			++it2;
			while (it2 != cl.end() && !found) {
				crs2=get_nonzero_vars(*it2);
				if (it2 != it && crs==crs2) {
					e2=ConstraintHom2Linear_Expression(&(*it2));
					if (get_nonzero_vars(e+e2).size()==0) {
						// found a match
						found=true;
					}
				}
				if (!found)
					++it2;
			}
			if (found) {
				// output the combined equations
				if (!first_constraint)
					s << " & ";
				else
					first_constraint=false;
				if (it->coefficient(Variable(*crs.rbegin()))>0) // count sign of the highest variable (->primed)
				{
					// ax+b >= 0 & -ax + c >= 0 --> -b <= ax <= c
					// the inhom term of the pos coeff is the lower bound
					first=true;
					print_inhomogeneous_term(s, *it, true, first); // reverse the sign
					s << " ";
					print_relation_symbol(s, it->type(), true); // reverse the sign
					s << " ";
					first=true;
					ppl_polyhedron::print_Linear_Expression_hom(s, e, vnvec, first);
					first=true;
					s << " ";
					print_relation_symbol(s, it2->type(), true); // reverse the sign
					s << " ";
					first=true;
					print_inhomogeneous_term(s, *it2, false, first); // reverse the sign
				} else {
					first=true;
					print_inhomogeneous_term(s, *it2, true, first); // reverse the sign
					s << " ";
					print_relation_symbol(s, it2->type(), true); // reverse the sign
					s << " ";
					first=true;
					ppl_polyhedron::print_Linear_Expression_hom(s, -e, vnvec, first);
					first=true;
					s << " ";
					print_relation_symbol(s, it->type(), true); // reverse the sign
					s << " ";
					first=true;
					print_inhomogeneous_term(s, *it, false, first); // reverse the sign
				}

				cl.erase(it2);
				it=cl.erase(it);
			} else
				++it;
		} else
			++it;
	}

	for (Constraint_List::const_iterator i = cl.begin(); i != cl.end(); i++) {
		if (!first_constraint)
			s << " & ";
		else
			first_constraint=false;
		ppl_polyhedron::print_constraint(s, *i, vnvec);
	}
}

int Constraint_List::get_memory() const {
	int m=0;
	Constraint_List::const_iterator i;
	for (i = begin(); i != end(); i++) {
		m+=i->total_memory_in_bytes();
	}
	return m;
}

clock_ref_set get_nonzero_vars(const Constraint& c) {
	clock_ref_set crs;
	for (dimension_type i=0; i<c.space_dimension(); ++i) {
		if (c.coefficient(Variable(i))!=0)
			crs.insert(i);
	}
	return crs;
}

clock_ref_set get_positive_vars(const Constraint& c) {
	clock_ref_set crs;
	for (dimension_type i=0; i<c.space_dimension(); ++i) {
		if (c.coefficient(Variable(i))>0)
			crs.insert(i);
	}
	return crs;
}

clock_ref_set get_nonzero_vars(const Linear_Expression& c) {
	clock_ref_set crs;
	for (dimension_type i=0; i<c.space_dimension(); ++i) {
		if (c.coefficient(Variable(i))!=0)
			crs.insert(i);
	}
	return crs;
}

clock_ref_set get_positive_vars(const Linear_Expression& c) {
	clock_ref_set crs;
	for (dimension_type i=0; i<c.space_dimension(); ++i) {
		if (c.coefficient(Variable(i))>0)
			crs.insert(i);
	}
	return crs;
}

/*
Constraint_List get_nonconsts_in_rel(const Constraint_List& cl) {
	Constraint_List ncl;
	clock_ref_set crs;
	clock_ref r1, r2;
	bool is_const;
	Constraint_List::const_iterator it=cl.begin();
	while (it!=cl.end()) {
		assert(cl.space_dimension()%2==0);
		is_const=false;
		crs=ppl_polyhedron::get_nonzero_vars(*it);
		if (crs.size()==2) {
			r1=*crs.begin();
			// get the second clock
			for (clock_ref_set::const_iterator ii=crs.begin(); ii!=crs.end(); ++ii)
				r2=*ii;
			if (r1+it->space_dimension()/2==r2 && it->coefficient(Variable(r1))
					==-it->coefficient(Variable(r2))) {
				// this constraint denotes that variavle r1 remains constant in the relation
				is_const=true;
			}
		}
		if (!is_const)
			ncl.push_back(*it);
		++it;
	}
	return ncl;
}
*/

// --------------------------------------------------------------
// from print.hh
// --------------------------------------------------------------

void print_constraint(const Constraint& c, const string& intro, ostream& s) {
	if (!intro.empty())
		s << intro << endl;
	s << c << endl;
}

void print_constraints(const Constraint_System& cs, const string& intro,
		ostream& s) {
	if (!intro.empty())
		s << intro << endl;
	Constraint_System::const_iterator i = cs.begin();
	Constraint_System::const_iterator cs_end = cs.end();
	bool printed_something = i != cs_end;
	while (i != cs_end) {
		s << *i++;
		if (i != cs_end)
			//      s << "," << endl; gf, 24.11.02
			s << " & "; // << endl;
	}
	s << (printed_something ? "." : "true.") << endl;
}

void print_constraints(const Polyhedron& ph, const string& intro, ostream& s) {
	ppl_polyhedron::print_constraints(ph.constraints(), intro, s);
}

// --------------------------------------------------------------


dimension_type Constraint_List::get_max_space_dimension() const {
	dimension_type d=0;
	for (Constraint_List::const_iterator i=begin(); i!=end(); ++i) {
		if (i->space_dimension() > d)
			d=i->space_dimension();
	}
	return d;
}

void Constraint_List::add_space_dimensions_and_embed(dimension_type m) {
	// rebuild the constraints
	for (Constraint_List::iterator i=begin(); i!=end(); ++i) {
		Constraint_add_space_dimensions_and_embed(*i, m);
	}
}

void Constraint_List::add_space_dimensions_and_embed_to_dim(dimension_type m)
// add dimensions to constraints such that all have dimension m
{
	// rebuild the constraints
	for (Constraint_List::iterator i=begin(); i!=end(); ++i) {
		Constraint_add_space_dimensions_and_embed(*i, m-i->space_dimension());
	}
}

void Constraint_List::dimension_move_assign(dimension_type x1,
		dimension_type x2, dimension_type y) {
	index_to_index_bimap pfunc;
	pfunc.move_assign(x1, x2, y, space_dimension());
	map_space_dimensions(pfunc);
}

// void
// Constraint_List::intersection_assign(const Constraint_List& cl)
// {
// 	for (Constraint_List::const_iterator i=cl.begin(); i!=cl.end(); ++i)
// 	{
// 		push_back(*i);
// 	}
// }

void Constraint_List::union_assign(const Constraint_List& cl) {
	for (Constraint_List::const_iterator i=cl.begin(); i!=cl.end(); ++i) {
		push_back(*i);
	}
}

void Constraint_List::reverse_first_half() {
	Constraint_List result(mydim);
	Linear_Expression e;

	// for reversing time: change the sign on the first dim/2 variables
	for (Constraint_List::iterator it=begin(); it!=end(); ++it) {
		e=0*Variable(it->space_dimension()-1);
		for (int i = it->space_dimension() - 1; i >= 0; i--) {
			if ((unsigned int)2*i>=it->space_dimension() )
				e += it->coefficient(Variable(i)) * Variable(i);
			else
				e -= it->coefficient(Variable(i)) * Variable(i);
		}
		e += it->inhomogeneous_term();
		if (it->is_equality()) {
			result.push_back(e==0);
		} else {
			if (it->is_strict_inequality())
				result.push_back(e > 0);
			else
				result.push_back(e >= 0);
		}
	}
	*this=result;

}

void Constraint_List::strictify() {
	// turn nonstrict into strict inequalities; equalities result in contradicting constraints
	for (Constraint_List::iterator it=begin(); it!=end(); ++it) {
		if (it->is_equality()) {
			*it=(Constraint2Linear_Expression(&*it)>0);
			push_front(Constraint2Linear_Expression(&*it)<0); // insert it in the beginning so it won't be strictified (although it wouldn't hurt)
		} else {
			if (!it->is_strict_inequality())
				*it=(Constraint2Linear_Expression(&*it)>0);
		}
	}
}

std::ostream& operator<<(std::ostream& os, const Constraint_List &cl) {
	for (Constraint_List::const_iterator i=cl.begin(); i!=cl.end(); ++i) {
		if (i!=cl.begin())
			os << ", ";
		os << *i;
	}
	return os;
}

void Constraint_add_space_dimensions_and_embed(Constraint& c, dimension_type m) {
	//	Linear_Expression e=0*Variable(c.space_dimension()-1+m)+Linear_Expression(c); // make sure it has the same dimension
	Linear_Expression e=0*Variable(c.space_dimension()-1+m)
			+Constraint2Linear_Expression(&c); // make sure it has the same dimension

	if (c.is_equality()) {
		c=(e==0);
	} else {
		if (c.is_strict_inequality())
			c=(e > 0);
		else
			c=(e >= 0);
	}
}

Linear_Expression Constraint2Linear_Expression(const Constraint* pcon,
		dimension_type start_dim) {
	// utility function: get the Linear_Expression part of a constraint
	// code taken from example in ppl documentation
	// constraint is e > 0 or e >= 0
	Linear_Expression e;

	// dimension_type is always >0
	// 	if (start_dim<0)
	//	throw_error("Constraint2Linear_Expression: start_dim < 0");

	for (int i = pcon->space_dimension() - 1; i >= (int)start_dim; i--)
		e += pcon->coefficient(Variable(i)) * Variable(i);
	e += pcon->inhomogeneous_term();
	return e;
}

Linear_Expression Constraint2Linear_Expression_moved_down(
		const Constraint* pcon, dimension_type start_dim) {
	// utility function: get the Linear_Expression part of a constraint
	// note: result includes inhomogeneous term!
	// code taken from example in ppl documentation
	// constraint is e > 0 or e >= 0
	Linear_Expression e;

	// dimension_type is always >0
	//	if (start_dim<0)
	//		throw_error("Constraint2Linear_Expression: start_dim < 0");

	for (int i = pcon->space_dimension() - 1; i >= (int)start_dim; i--)
		e += pcon->coefficient(Variable(i)) * Variable(i-start_dim);
	e += pcon->inhomogeneous_term();
	return e;
}

Linear_Expression Constraint2Linear_Expression_up_to(const Constraint* pcon,
		dimension_type stop_dim) {
	// utility function: get the Linear_Expression part of a constraint
	// up to Variable(stop_dim).
	// Attention: Linear_Expression is then of dimension (stop_dim+1) !!!

	if (stop_dim>pcon->space_dimension() - 1)
		throw_error("Constraint2Linear_Expression_up_to: stop_dim > max");

	Linear_Expression e;
	for (int i = stop_dim; i >= 0; i--)
		e += pcon->coefficient(Variable(i)) * Variable(i);
	e += pcon->inhomogeneous_term();
	return e;
}

Linear_Expression ConstraintHom2Linear_Expression(const Constraint* pcon) {
	// utility function: get the Linear_Expression part of a constraint
	// inhomogenous part is dropped
	Linear_Expression e;
	for (int i = pcon->space_dimension() - 1; i >= 0; i--)
		e += pcon->coefficient(Variable(i)) * Variable(i);
	return e;
}

Linear_Expression ConstraintHom2Linear_Expression_up_to(const Constraint* pcon,
		dimension_type stop_dim) {
	// utility function: get the Linear_Expression part of a constraint
	// inhomogenous part is dropped

	if (stop_dim>pcon->space_dimension() - 1)
		throw_error("ConstraintHom2Linear_Expression_up_to: stop_dim > max");

	Linear_Expression e;
	for (int i = stop_dim; i >= 0; i--)
		e += pcon->coefficient(Variable(i)) * Variable(i);
	return e;
}

Constraint_List Constraint_System_convert_equalities(const Constraint_System& cs) {
	// utility function: replace each equality by two inequalities
	// simply copy other constraints
	// Comment: Constraint_List type has to be used because a Constraint_System would just
	//          be minimized to equalities again
	Constraint_List csl;
	Constraint_System::const_iterator i;
	for (i = cs.begin(); i != cs.end(); i++) {
		if (i->is_equality()) {
			Linear_Expression e = Constraint2Linear_Expression(&(*i));
			csl.push_back(e<=0);
			csl.push_back(e>=0);
		} else {
			csl.push_back(*i);
		}
	}
	return csl;
}

Constraint_List Constraint_System_no_equalities(const Constraint_System& cs) {
	// utility function: copy constraints that are not equalities
	Constraint_List csl;
	Constraint_System::const_iterator i;
	for (i = cs.begin(); i != cs.end(); i++) {
		if (!(i->is_equality())) {
			csl.push_back(*i);
		}
	}
	return csl;
}

Constraint closed_inequality_complement(const Constraint& cs) {
	// utility function: return the closure of the complement of cs
	// Comment: only valid if cs is an inequality

	if (cs.is_equality()) {
		throw_error("attempting to complement equality constraint in inequality_complement");
		return Constraint::zero_dim_false();
	} else {
		Linear_Expression e = Constraint2Linear_Expression(&cs);
		return Constraint(e <= 0);
	}
}

Constraint inequality_complement(const Constraint& cs) {
	// utility function: replace the inequality by its complement
	// Comment: only valid if cs is an inequality

	if (cs.is_equality()) {
		throw_error("attempting to complement equality constraint in inequality_complement");
		return Constraint::zero_dim_false();
	} else {
		Linear_Expression e = Constraint2Linear_Expression(&cs);
		if (cs.is_strict_inequality())
			return Constraint(e <= 0);
		else
			return Constraint(e < 0);
	}
}

Constraint equality_constraint(const Constraint& cs) {
	// utility function: turn cs into an equality constraint

	if (cs.is_equality()) {
		return cs;
	} else {
		Linear_Expression e = Constraint2Linear_Expression(&cs);
		return Constraint(e == 0);
	}
}

Constraint_List complement(const Constraint& cs) {
	// utility function: replace the (in)equality by its complement
	// Comment: This results in a list of constraints that's not a Constraint_System.
	Constraint_List csl;
	Linear_Expression e = Constraint2Linear_Expression(&cs);
	if (cs.is_equality()) {
		csl.push_back(e<0);
		csl.push_back(e>0);
	} else {
		if (cs.is_strict_inequality())
			csl.push_back(e <= 0);
		else
			csl.push_back(e < 0);
	}
	return csl;
}

Constraint_List complement(const Constraint_System& cs) {
	// utility function: replace each (in)equality by its complement
	// Comment: This results in a list of constraints that's not a Constraint_System.
	Constraint_List csl;
	Constraint_System::const_iterator i;
	for (i = cs.begin(); i != cs.end(); i++) {
		Linear_Expression e = Constraint2Linear_Expression(&(*i));
		if (i->is_equality()) {
			csl.push_back(e<0);
			csl.push_back(e>0);
		} else {
			if (i->is_strict_inequality())
				csl.push_back(e <= 0);
			else
				csl.push_back(e < 0);
		}
	}
	return csl;
}

Constraint constraint_homogenous_part(const Constraint* pcon) {
	Constraint c(Variable(0)==0); // initialize with some arbitrary value
	Linear_Expression e = Constraint2Linear_Expression(pcon);
	if (pcon->is_equality()) {
		c=(e-pcon->inhomogeneous_term()==0);
	} else {
		if (pcon->is_strict_inequality())
			c=(e-pcon->inhomogeneous_term()>0);
		else
			c=(e-pcon->inhomogeneous_term()>=0);
	}
	return c;
}

Constraint constraint_homogenous_part(const Constraint& con) {
	Constraint c(Variable(0)==0); // initialize with some arbitrary value
	Linear_Expression e = Constraint2Linear_Expression(&con);
	if (con.is_equality()) {
		c=(e-con.inhomogeneous_term()==0);
	} else {
		if (con.is_strict_inequality())
			c=(e-con.inhomogeneous_term()>0);
		else
			c=(e-con.inhomogeneous_term()>=0);
	}
	return c;
}

Constraint constraint_to_equality(const Constraint& cs) {
	Constraint c(Variable(0)==0); // initialize with some arbitrary value
	Linear_Expression e = Constraint2Linear_Expression(&cs);
	c=(e == 0);
	return c;
}

Constraint constraint_to_nonstrict_inequality(const Constraint& cs) {
	Constraint c(Variable(0)==0); // initialize with some arbitrary value
	Linear_Expression e = Constraint2Linear_Expression(&cs);
	c=(e >= 0);
	return c;
}

Constraint constraint_closure(const Constraint& cs) {
	// return the constraint closed
	if (cs.is_strict_inequality()) {
		Linear_Expression e = Constraint2Linear_Expression(&cs);
		return Constraint(e>=0);
	} else
		return cs;
}

int Constraint_max_bit_size(const Constraint& c) {
	dimension_type dum;
	dimension_type m=0;

	for (dimension_type j=0; j<c.space_dimension(); ++j) {
		dum=mpz_sizeinbase(c.coefficient(Variable(j)).get_mpz_t(),2);
		if (dum>m)
			m=dum;
	}
	return (int)m;
}

Constraint Constraint_restrict_to_dim(const Constraint& cs,
		dimension_type stop_dim) {
	// restrict the constraint to stop_dim dimensions.
	Constraint c(Variable(0)==0); // initialize with some arbitrary value
	Linear_Expression e = ppl_polyhedron::Constraint2Linear_Expression_up_to(&cs, stop_dim-1);
	if (cs.is_equality()) {
		c=(e==0);
	} else {
		if (cs.is_strict_inequality())
			c=(e>0);
		else
			c=(e>=0);
	}
	return c;
}

int Constraint_count_coefficients(const Constraint& c) {
	// return the number of nonzero inhomogeneous terms
	int n=0;
	for (dimension_type j=0; j<c.space_dimension(); ++j) {
		if (c.coefficient(Variable(j))!=Integer(0)) {
			++n;
		}
	}
	return n;
}

bool Constraint_is_carthesian(const Constraint& c) {
	// return if the constraint has at most one homogeneous coefficient
	bool found_one = false;
	for (dimension_type j=0; j<c.space_dimension(); ++j) {
		if (c.coefficient(Variable(j))!=Integer(0)) {
			if (found_one)
				return false;
			else
				found_one=true;
		}
	}
	return true;
}

bool Constraint_is_octagonal(const Constraint& c) {
	// return if the constraint has at most one homogeneous coefficient
	if (ppl_polyhedron::Constraint_count_coefficients(c)==2) {
		dimension_type first=c.space_dimension();
		for (dimension_type j=0; j<c.space_dimension(); ++j) {
			if (c.coefficient(Variable(j))!=Integer(0)) {
				if (first==c.space_dimension())
					first=j;
				else if (c.coefficient(Variable(first))
						==c.coefficient(Variable(j))
						|| c.coefficient(Variable(first))
								==-c.coefficient(Variable(j)))
					return true;
				else
					return false;
			}
		}
		return false;
	} else
		return false;
}

void print_constraints(const Constraint_List& C) {
	Constraint_List::const_iterator i;
	int counter=0;
	for (i = C.begin(); i != C.end(); ++i) {
		cout << ++counter <<": ";
		ppl_polyhedron::print_constraint(*i, "", cout);
		// new  cout << (*i);
	}
}

// ---------------------------------------------------------
// ---------------------------------------------------------
// misc. polyhedral functions
// ---------------------------------------------------------
// ---------------------------------------------------------

dimension_type get_true_dimension(const Linear_Expression& linex) { // Return the highest non-zero coefficient + 1
//cout << linex << endl;
	if (linex.space_dimension()>0) {
		int i=linex.space_dimension()-1;
		while ((i>=0) && (linex.coefficient(Variable((dimension_type)i))==0))
			--i;
		//cout << ":" << (dimension_type)i;
		return (dimension_type)(i+1);
	} else
		return 0;
}

dimension_type get_true_dimension(const Constraint& linex) { // Return the highest non-zero coefficient + 1
	if (linex.space_dimension()>0) {
		int i=linex.space_dimension()-1;
		while ((i>=0) && (linex.coefficient(Variable((dimension_type)i))==0))
			--i;
		return (dimension_type)(i+1);
	} else
		return 0;
}

bool attempt_multiply(Linear_Expression& l1, Linear_Expression& l2) { // multiply if l1 or l2 is a scalar and put the result in l1. Return true if successful.

	if (ppl_polyhedron::get_true_dimension(l1)==0) {
		l1=l1.inhomogeneous_term()*l2;
		return true;
	} else if (ppl_polyhedron::get_true_dimension(l2)==0) {
		l2=l2.inhomogeneous_term()*l1;
		return true;
	} else
		return false;
}

}

