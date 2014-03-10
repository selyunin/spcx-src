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
#ifndef GUARD_extended_ppl_h
#define GUARD_extended_ppl_h

//#include <stdio.h>
//#include <istd::ostream>
//#include <vector>
#include <list>
//#include <set>
//#include <map>
//#include <algorithm>
//#include <stdexcept>
//#include <stdlib.h>

#include <ppl.hh>

#include "core/continuous/polyhedra/ppl_polyhedron/general.h"
#include "core/continuous/polyhedra/ppl_polyhedron/variable_id_map.h"
//#include "myPFunction.h"



//using namespace Parma_Polyhedra_Library;
//using namespace Parma_Polyhedra_Library::IO_Operators;
//using namespace std;

namespace ppl_polyhedron {

Parma_Polyhedra_Library::Linear_Expression Linear_Expression_subset(
		const Parma_Polyhedra_Library::Linear_Expression& o,
		const clock_ref_set& crs);
void print_Linear_Expression_hom(std::ostream& s,
		const Parma_Polyhedra_Library::Linear_Expression& c,
		const variable_id_map& vnvec, bool& first);
void print_relation_symbol(std::ostream& s,
		Parma_Polyhedra_Library::Constraint::Type t, bool reverse_sign);
void print_inhomogeneous_term(std::ostream& s,
		const Parma_Polyhedra_Library::Constraint& c, bool reverse_sign,
		bool& first);

// --------------------------------------------------------------
// from print.hh
// --------------------------------------------------------------
void print_constraint(const Parma_Polyhedra_Library::Constraint& c,
		const std::string& intro, std::ostream& s);

void print_constraints(const Parma_Polyhedra_Library::Constraint_System& cs,
		const std::string& intro, std::ostream& s);

void print_constraints(const Parma_Polyhedra_Library::Polyhedron& ph,
		const std::string& intro, std::ostream& s);

// --------------------------------------------------------------


Parma_Polyhedra_Library::Linear_Expression Constraint2Linear_Expression(
		const Parma_Polyhedra_Library::Constraint* pcon,
		Parma_Polyhedra_Library::dimension_type start_dim = 0);
Parma_Polyhedra_Library::Linear_Expression
		Constraint2Linear_Expression_moved_down(
				const Parma_Polyhedra_Library::Constraint* pcon,
				Parma_Polyhedra_Library::dimension_type start_dim);
Parma_Polyhedra_Library::Linear_Expression Constraint2Linear_Expression_up_to(
		const Parma_Polyhedra_Library::Constraint* pcon,
		Parma_Polyhedra_Library::dimension_type stop_dim);
Parma_Polyhedra_Library::Linear_Expression ConstraintHom2Linear_Expression(
		const Parma_Polyhedra_Library::Constraint* pcon);
Parma_Polyhedra_Library::Linear_Expression
		ConstraintHom2Linear_Expression_up_to(
				const Parma_Polyhedra_Library::Constraint* pcon,
				Parma_Polyhedra_Library::dimension_type stop_dim);

void print_constraint(std::ostream& s, const Parma_Polyhedra_Library::Constraint& c,
		const variable_id_map& vnvec);

// -------------------------------------------------------------------------------
// Constraint List
// -------------------------------------------------------------------------------
// - necessary because Constraint_System is automatically simplified,
//   so inequalities can turn into equalities etc.

class Constraint_List : public std::list<Parma_Polyhedra_Library::Constraint> {
public:
	Constraint_List() :
		mydim(0) {
	}
	;
	Constraint_List(Parma_Polyhedra_Library::dimension_type dim) :
		mydim(dim) {
	}
	;

	bool is_empty() const {
		return empty();
	}
	;
	Parma_Polyhedra_Library::dimension_type get_max_space_dimension() const;
	Parma_Polyhedra_Library::dimension_type space_dimension() const {
		return get_max_space_dimension();
	}
	;

	void add_space_dimensions_and_embed(
			Parma_Polyhedra_Library::dimension_type m);
	void add_space_dimensions_and_embed_to_dim(
			Parma_Polyhedra_Library::dimension_type m);
	void dimension_move_assign(Parma_Polyhedra_Library::dimension_type x1,
			Parma_Polyhedra_Library::dimension_type x2,
			Parma_Polyhedra_Library::dimension_type y);

	//	void intersection_assign(const Constraint_List& cl);
	void union_assign(const Constraint_List& cl);

	template<typename PartialFunction> void map_space_dimensions(
			PartialFunction pfunc) {
		for (Constraint_List::iterator i=begin(); i!=end(); ++i) {
			Constraint_map_space_dimensions(*i, pfunc);
		}
	}

	void constraint_to_equality(Constraint_List::iterator i) {
		push_back(Constraint2Linear_Expression(&(*i))<=0); // make copy
		*i=(Constraint2Linear_Expression(&(*i))>=0);
	}
	;

	void reverse_first_half();
	void strictify();

	void print_phaver(std::ostream& s, const variable_id_map& vnvec) const;
	int get_memory() const;

private:
	Parma_Polyhedra_Library::dimension_type mydim; // CAREFUL!!!!!! mydim is not maintained, don't use it
};

std::ostream& operator<<(std::ostream& os, const Constraint_List &cl);

void Constraint_add_space_dimensions_and_embed(
		Parma_Polyhedra_Library::Constraint& c,
		Parma_Polyhedra_Library::dimension_type m);

template<typename PartialFunction> void Constraint_map_space_dimensions(
		Parma_Polyhedra_Library::Constraint& c, PartialFunction pfunc) {
	// rebuild the constraint
	Parma_Polyhedra_Library::Linear_Expression e=0
			*Parma_Polyhedra_Library::Variable(c.space_dimension()-1); // make sure it has the same dimension

	for (int i = c.space_dimension() - 1; i >= 0; i--)
		e += c.coefficient(Parma_Polyhedra_Library::Variable(i))
				* Parma_Polyhedra_Library::Variable(pfunc.get_map(i));
	e += c.inhomogeneous_term();

	if (c.is_equality()) {
		c=(e==0);
	} else {
		if (c.is_strict_inequality())
			c=(e > 0);
		else
			c=(e >= 0);
	}
}
;

Constraint_List Constraint_System_convert_equalities(
		const Parma_Polyhedra_Library::Constraint_System& cs);

Constraint_List Constraint_System_no_equalities(
		const Parma_Polyhedra_Library::Constraint_System& cs);

Parma_Polyhedra_Library::Constraint closed_inequality_complement(
		const Parma_Polyhedra_Library::Constraint& cs);

Parma_Polyhedra_Library::Constraint inequality_complement(
		const Parma_Polyhedra_Library::Constraint& cs);

Parma_Polyhedra_Library::Constraint equality_constraint(
		const Parma_Polyhedra_Library::Constraint& cs);

Constraint_List complement(const Parma_Polyhedra_Library::Constraint& cs);

Constraint_List
		complement(const Parma_Polyhedra_Library::Constraint_System& cs);

int Constraint_count_coefficients(const Parma_Polyhedra_Library::Constraint& c);

bool Constraint_is_carthesian(const Parma_Polyhedra_Library::Constraint& c);

bool Constraint_is_octagonal(const Parma_Polyhedra_Library::Constraint& c);

void print_constraints(const Constraint_List& C);

// -------------------------------------------------------------------------------
// Generator List
// -------------------------------------------------------------------------------
// - necessary because Generator_System is automatically simplified,
//   so interior points can disappear

typedef std::list<Parma_Polyhedra_Library::Generator> Generator_List;

// -------------------------------------------------------------------------------

Parma_Polyhedra_Library::Constraint constraint_homogenous_part(
		const Parma_Polyhedra_Library::Constraint* pcon);

Parma_Polyhedra_Library::Constraint constraint_homogenous_part(
		const Parma_Polyhedra_Library::Constraint& con);

Parma_Polyhedra_Library::Constraint constraint_to_equality(
		const Parma_Polyhedra_Library::Constraint& cs);

Parma_Polyhedra_Library::Constraint constraint_to_nonstrict_inequality(
		const Parma_Polyhedra_Library::Constraint& cs);

Parma_Polyhedra_Library::Constraint constraint_closure(
		const Parma_Polyhedra_Library::Constraint& cs);

int Constraint_max_bit_size(const Parma_Polyhedra_Library::Constraint& c);

Parma_Polyhedra_Library::Constraint Constraint_restrict_to_dim(
		const Parma_Polyhedra_Library::Constraint& cs,
		Parma_Polyhedra_Library::dimension_type stop_dim);

// ---------------------------------------------------------
// ---------------------------------------------------------
// misc. polyhedral functions
// ---------------------------------------------------------
// ---------------------------------------------------------

Parma_Polyhedra_Library::dimension_type get_true_dimension(
		const Parma_Polyhedra_Library::Linear_Expression& linex);

Parma_Polyhedra_Library::dimension_type get_true_dimension(
		const Parma_Polyhedra_Library::Constraint& linex);

bool attempt_multiply(Parma_Polyhedra_Library::Linear_Expression& l1,
		Parma_Polyhedra_Library::Linear_Expression& l2);

// -------------------------------------------------------------------------
// Enhancing Partial Functions
// -------------------------------------------------------------------------

template<typename PartialFunction> void PartialFunction_Double(
		PartialFunction& pfunc, Parma_Polyhedra_Library::dimension_type tdim) {
	// Repeat the partial function for higher dimensions, n=tdim
	// i.e. 0->pfunc(0), ..., (n-1)->pfunc(n-1)
	// turns to
	//      0->pfunc(0), ..., (n-1)->pfunc(n-1),
	//      n->pfunc(0)+n, ..., (2n-1)->pfunc(n-1)+n

	// tdim is the dimension of the domain of pfunc!

	dimension_t ndim=pfunc.max_in_codomain()+1;
	for (Parma_Polyhedra_Library::dimension_type i=0; i<tdim; ++i) {
		if (pfunc.in_domain(i))
			pfunc.insert(i+tdim, pfunc.get_map(i)+ndim);
	};
}
;

// -------------------------------------------------------------------------

clock_ref_set get_nonzero_vars(const Parma_Polyhedra_Library::Constraint& c);
clock_ref_set get_positive_vars(const Parma_Polyhedra_Library::Constraint& c);
clock_ref_set get_nonzero_vars(
		const Parma_Polyhedra_Library::Linear_Expression& c);
clock_ref_set get_positive_vars(
		const Parma_Polyhedra_Library::Linear_Expression& c);
Constraint_List get_nonconsts_in_rel(const Constraint_List& cl);

}

#endif
