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

#include "core/continuous/polyhedra/ppl_polyhedron/fp_interface.h"

#include <vector>
#include <list>
#include <functional>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include <stdexcept>

#include <ppl.hh>
#include <time.h>

namespace ppl_polyhedron {

using namespace std;
using namespace Parma_Polyhedra_Library;

const unsigned int GENERATOR_TO_DOUBLE_PRECISION=10;

double Round(double Zahl, const unsigned int& Stellen) {
	// if (Zahl-floor(Zahl * pow((double)10, Stellen) + 0.5) * pow((double)10, -Stellen)!=0)
	//	cout << Zahl << "->" << floor(Zahl * pow((double)10, Stellen) + 0.5) * pow((double)10, -Stellen) << endl;

	if (Stellen > 0)
		return floor(Zahl * pow((double)10, (double)Stellen) + 0.5) * pow(
				(double)10, -(double)Stellen);
	else
		return Zahl;
}

double generator_to_double(Generator& g, dimension_type pos) {
	mpq_class q(g.coefficient(Variable(pos)), g.divisor());
	return Round(q.get_d(), GENERATOR_TO_DOUBLE_PRECISION);
}

void generator_to_double_point(const Generator& g, double_point& p) {
	if (p.size()!=g.space_dimension())
		p=double_point(g.space_dimension());
	if (g.space_dimension()>0) {
		/*   mpq_class q(g.coefficient(Variable(0)),g.divisor());
		 p[0]=q.get_d();
		 for (int i=1;i<g.space_dimension();++i)
		 {
		 q.get_num()=g.coefficient(Variable(i));
		 p[i]=q.get_d();
		 }*/
		mpq_class q(0);
		for (int i=0; i<(int)g.space_dimension(); ++i) {
			q=mpq_class(g.coefficient(Variable(i)), g.divisor());
			p[i]=Round(q.get_d(), GENERATOR_TO_DOUBLE_PRECISION);
		}
	}
}

void get_max_coeff(const Linear_Expression& g, Integer& d) {
	// get maximum coefficient of linear expression

	d=Integer(0);
	for (int i=0; i<(int)g.space_dimension(); ++i) {
		if ((g.coefficient(Variable(i))>d) || (g.coefficient(Variable(i))<-d)) {
			d=abs(g.coefficient(Variable(i)));
		}
	}

	if ((g.inhomogeneous_term()>d) || (g.inhomogeneous_term()<-d)) {
		d=abs(g.inhomogeneous_term());
	}
}

void get_max_coeff(const Constraint& g, Integer& d) {
	// get maximum coefficient of constraint

	d=Integer(0);
	for (int i=0; i<(int)g.space_dimension(); ++i) {
		if ((g.coefficient(Variable(i))>d) || (g.coefficient(Variable(i))<-d)) {
			d=abs(g.coefficient(Variable(i)));
		}
	}

	if ((g.inhomogeneous_term()>d) || (g.inhomogeneous_term()<-d)) {
		d=abs(g.inhomogeneous_term());
	}
}

void hom_linexpression_to_double_point(const Linear_Expression& g,
		double_point& p) {
	// return the homogeneous art of a linear expression as a double_point

	if (p.size()!=g.space_dimension())
		p=double_point(g.space_dimension());

	if (g.space_dimension()>0) {

		// get maximum coefficient so the numbers don't get too big for double
		Integer d(1);
		// get maximum homoegeneous coefficient > 1 of linear expression
		for (int i=0; i<(int)g.space_dimension(); ++i) {
			if ((g.coefficient(Variable(i))>d) || (g.coefficient(Variable(i))
					<-d)) {
				d=abs(g.coefficient(Variable(i)));
			}
		}

		mpq_class q(0);
		for (int i=0; i<(int)g.space_dimension(); ++i) {
			q=mpq_class(g.coefficient(Variable(i)), d);
			p[i]=Round(q.get_d(), GENERATOR_TO_DOUBLE_PRECISION);
		}
	}
}

void constraint_to_double_point(const Constraint& c, double_point& dx, double& b) {
	// Returns the constraint as a normal vector in dx and b separately
	// the 2-norm of dx is 1 if dx!=0.
	if (dx.size()!=c.space_dimension()) {
		dx=double_point(c.space_dimension());
		for (unsigned int i=0; i<dx.size(); ++i) {
			dx[i]=0.0;
		}
	}

	if (c.space_dimension()>0) {
		// get maximum coefficient so the numbers don't get too big for double
		Integer d(1);
		get_max_coeff(c, d); // get max coefficient
		if (d==Integer(0))
			d=Integer(1);

		mpq_class q(0);
		for (unsigned int i=0; i<c.space_dimension(); ++i) {
			q.get_num()=c.coefficient(Variable(i));
			q.get_den()=d;
			//			dx[i]=Round(q.get_d(),GENERATOR_TO_DOUBLE_PRECISION); //q.get_d();
			dx[i]=q.get_d();
		}
		// inhomogenous term
		q.get_num()=c.inhomogeneous_term();
		q.get_den()=d;
		//		b=Round(q.get_d(),GENERATOR_TO_DOUBLE_PRECISION); //q.get_d();
		b=q.get_d();

		// norm the vector
		double length(sqrt(norm2(dx))); // get length of vector
		if (length>0) {
			dx/=length;
			b/=length;
		}
	}
}

void constraint_to_double_point(const Constraint& c, double_point& dx) {
	// Returns the constraint as a normal vector in dx=[a_1 ... a_n b]
	if (dx.size()!=c.space_dimension()+1)
		dx=double_point(c.space_dimension()+1);
		for (unsigned int i=0; i<dx.size(); ++i) {
			dx[i]=0.0;
		}

	if (c.space_dimension()>0) {
		// get maximum coefficient so the numbers don't get too big for double
		Integer d(1);
		get_max_coeff(c, d); // get max coefficient
		if (d==Integer(0))
			d=Integer(1);

		mpq_class q(0);
		for (unsigned int i=0; i<c.space_dimension(); ++i) {
			q.get_num()=c.coefficient(Variable(i));
			q.get_den()=d;
			//			dx[i]=Round(q.get_d(),GENERATOR_TO_DOUBLE_PRECISION); //q.get_d();
			dx[i]=q.get_d();
		}
		// normalize dx before adding the inhomogeneous coefficient
		dx[c.space_dimension()]=0.0;
		// norm the vector
		double length(sqrt(norm2(dx))); // get length of vector
		if (length>0) {
			dx/=length;
		}

		// inhomogenous term
		q.get_num()=c.inhomogeneous_term();
		q.get_den()=d;
		//		dx[c.space_dimension()]=Round(q.get_d(),GENERATOR_TO_DOUBLE_PRECISION); //q.get_d();
		dx[c.space_dimension()]=q.get_d();
		if (length>0) {
			dx[c.space_dimension()]/=length;
		}
	}
}

void double_point_to_Linear_Expression(const double_point& p,
		Linear_Expression& lin, Integer& d) {
	// convert double vector p to linear expression and denominator d
	mpq_class q(0);

	vector <Integer> num(p.size(), 0);
	vector <Integer> den(p.size(), 0);
	Integer gcd=1;
	lin=Linear_Expression(0*Variable(p.size()-1));
	d=Integer(1);

	if (p.size()>0) {
		// extract numerators and common denominator
		for (unsigned int i=0; i<p.size(); ++i) {
			// For each dimension
			q=mpq_class(p[i]);
			// expand the denominator
			d*=q.get_den();
			num[i]=q.get_num();
			den[i]=q.get_den();
		}

		// multiply all numerators with the expansion = the rest of the denominators
		for (unsigned int i=0; i<p.size(); ++i) {
			for (unsigned int j=0; j<p.size(); ++j) {
				if (i!=j)
					num[i]*=den[j];
			}
		}

		// could this come before the multiplication?
		// find the greatest common divisor of the numerator
		/*    gcd=num[0];
		 for (unsigned int i=1;i<p.size();++i)
		 {
		 gcd_assign(gcd,num[i]);
		 }
		 gcd_assign(gcd,d); // now it's the gcd between the numerators and the total den

		 if (gcd>1)
		 {
		 // divide all numerators by gcd
		 for (unsigned int i=0;i<p.size();++i)
		 {
		 num[i]/=gcd;
		 }
		 d/=gcd;
		 }*/

		// assign to linear expression
		for (unsigned int i=0; i<p.size(); ++i) {
			lin+=num[i]*Variable(i);
		}

	}
}

void double_point_to_generator(const double_point& p, Generator& g) {
	Linear_Expression lin;
	Integer d=1;
	double_point_to_Linear_Expression(p, lin, d);

	mpq_class q(0);

	vector <Integer> num(p.size(), 0);
	vector <Integer> den(p.size(), 0);

	// final assignment:
	// if dimension is zero, then lin=0, so result is correct
	g=point(lin, d);

	// Quality Test:
	/*double_point p2;
	 generator_to_double_point(g, p2);
	 double a=norm2(p2-p);
	 if (true) {
	 cout << "Error: "<<a<<" from ";
	 cout << p << " to " << p2;}*/
}

void add_Generator_List_to_double_point_list(const Generator_List& gl,
		double_point_list& pl) {
	// copy elements
	// ATTENTION: pl is not cleared, points will be added!

	double_point p(0);

	for (Generator_List::const_iterator i=gl.begin(); i!=gl.end(); ++i) {
		if (i->is_point() || i->is_closure_point()) {
			generator_to_double_point(*i, p);
			pl.push_back(p);
		}
	}
}

void get_constraint_through_double_point(Linear_Expression linex,
		double_point& p, Constraint& c, bool strict) {
	// Returns in c the constraint given by linex (sign) linex*p, i.e. the plane defined by linex that
	// goes through the point p.
	// The (sign) is "<" is strict==true, and "<=" othewise.
	// if linex=sum(a_i*x_i), then c:= sum(d*a_i*x_i)-sum(a_i*c_i)==0.
	// Attention: dimension is taken from linex!!!

	// First, convert p into a generator to get at the integer coefficients
	// todo: this can probably be done more efficiently
	Generator g=point(0*Variable(0), 1); // dummy
	double_point_to_generator(p, g);

	Integer scale(1); //("1000000000000000000000000000");

	// get sum(a_i*c_i)
	Integer b=0;
	for (unsigned int i=0; i<linex.space_dimension(); ++i) {
		b+=linex.coefficient(Variable(i))*g.coefficient(Variable(i));
	}
	//	while (b/scale>1024*1024)
	// 	while (b>scale*1024*1024 || b<-scale*1024*1024)
	// 		scale=scale*10;
	// 	b=b/scale;

	// linex'=sum(d*a_i*x_i)
	linex*=g.divisor()/scale;

	if (strict)
		c=Constraint(linex-b > 0);
	else
		c=Constraint(linex-b >= 0);
}

void get_constraint_through_double_point(double_point& dx, double_point& p,
		Constraint& c, bool strict) {
	// returns the constraint given by sum(dx[i]*x_i)+b with b such that the plane goes through p.
	// Attention: dimension is taken from dx!!!

	// convert dx into a generator to get at the integer coefficients
	Generator gd=point(0*Variable(0), 1); // dummy
	double_point_to_generator(dx, gd);

	Linear_Expression linex;
	for (unsigned int i=0; i<dx.size(); ++i) {
		linex+=gd.coefficient(Variable(i))*Variable(i);
	}

	get_constraint_through_double_point(linex, p, c, strict);
}

void get_constraint_through_double(double_point& dx, double q, Constraint& c,
		bool strict) {
	// returns the constraint given by sum(dx[i]*x_i)+b with b=q.
	// Attention: dimension is taken from dx!!!

	// find a variable with non-zero coefficient to determine b
	unsigned int i=0;
	while (dx[i]==0.0 && i<dx.size())
		++i;

	if (i<dx.size() && dx[i]!=0.0) {
		// dx[i]*x=q	-> x=q/dx[i]
		double_point p(dx.size());

		p[i]=-q/dx[i];
		get_constraint_through_double_point(dx, p, c, strict);
	}
}

// Computation Methods


void get_double_point_list_center(double_point_list& pl, double_point& p) {
	// Computes the arithmetic center of the points in pl.
	// Returns a zero point if pl is empty.
	// ATTENTION: the dimension of p must be equal to the dimension of the points in pl

	p=double_point(p.size());
	for (unsigned int i=0; i<p.size(); ++i) {  // initialize p to zero
		p[i]=0.0;
	}

	if (pl.size()>0) {
		for (double_point_list::const_iterator i=pl.begin(); i!=pl.end(); ++i)
			p+=*i;
		for (unsigned int i=0; i<p.size(); ++i)
			p[i]=p[i]/pl.size();
	}
}

void get_double_point_list_min_max(const double_point_list& pl,
		const double_point& dx, double& min_dx, double& max_dx,
		double_point& x_min, double_point& x_max) {
	// find the min and max of pl[i]*dx
	double v=0;
	bool no_min=true;
	bool no_max=true;

	for (double_point_list::const_iterator i=pl.begin(); i!=pl.end(); ++i) {
		v=scalar_product(*i, dx);
		if (v<min_dx || no_min) {
			min_dx=v;
			no_min=false;
			x_min=*i;
		}
		if (v>max_dx || no_max) {
			max_dx=v;
			no_max=false;
			x_max=*i;
		}
	}
}

double get_double_point_list_angle(double_point_list& pl) {
	// Computes the minimum cos(theta) between any two points in pl
	// Returns zero if pl is empty.
	// ATTENTION: the dimension of p must be equal to the dimension of the points in pl

	double p, n;
	double pmin=1;

	double_point_list::const_iterator j;
	if (pl.size()>0) {
		for (double_point_list::const_iterator i=pl.begin(); i!=pl.end(); ++i) {
			j=i;
			++j;
			while (j!=pl.end()) {
				// angle is equal to a^T b / |a||b|
				n=sqrt(norm2(*i))*sqrt(norm2(*j));
				if (n!=0)
					p=scalar_product(*i, *j)/n;
				else
					p=1;

				if (p<pmin)
					pmin=p;

				++j;
			}
		}
		//cout << pmin << endl;
		return pmin;
	}

	return 0;
}

double get_double_point_list_angle_vecs(double_point_list& pl,
		double_point& g1, double_point& g2) {
	// Computes the minimum cos(theta) between any two points in pl
	// Returns the two points that have the maximal angle
	// Returns one if pl is empty.
	// ATTENTION: the dimension of p must be equal to the dimension of the points in pl

	double p, n;
	double pmin=1;

	if (pl.size()>1) {
		g1=*pl.begin();
		g2=*pl.begin();
		double_point_list::const_iterator j;
		for (double_point_list::const_iterator i=pl.begin(); i!=pl.end(); ++i) {
			j=i;
			++j;
			while (j!=pl.end()) {
				// angle is equal to a^T b / |a||b|
				n=sqrt(norm2(*i))*sqrt(norm2(*j));
				if (n!=0)
					p=scalar_product(*i, *j)/n;
				else
					p=1;

				if (p<pmin) {
					g1=*i;
					g2=*j;
					pmin=p;
				}

				++j;
			}
		}
		//cout << pmin << endl;
		return pmin;
	}

	return 1;
}

double get_double_point_list_variance(double_point_list& pl) {
	// Computes the arithmetic center of the points in pl.
	// Returns a zero point if pl is empty.
	// ATTENTION: the dimension of p must be equal to the dimension of the points in pl

	double_point p(pl.begin()->size()); // initialize p to zero
	for (unsigned int i=0; i<p.size(); ++i) {  // initialize p to zero
		p[i]=0.0;
	}

	get_double_point_list_center(pl, p);

	double v=0;

	if (pl.size()>0) {
		for (double_point_list::const_iterator i=pl.begin(); i!=pl.end(); ++i) {
			v+=norm2<double>(p-(*i));
		}
	}
	//cout << v;

	return v;
}

void print_generator_fp_raw(ostream& s, const Generator& g) {
	const int num_variables = g.space_dimension();
	//  if (vnvec.size()<num_variables)
	//    throw_error("Not enough variable names for constraint.")
	//  bool first = true;

	// for now: only print points and closure points
	// todo: do something about rays and lines
	if (g.is_point() || g.is_closure_point()) {
		Integer den = g.divisor();
		Integer cv;

		//  double mpq_get_d (mpq_t OP)
		///  void mpq_div (mpq_t QUOTIENT, mpq_t DIVIDEND, mpq_t
		//          DIVISOR)

		// - Function: void mpq_set_num (mpq_t RATIONAL, mpz_t NUMERATOR)
		// - Function: void mpq_set_den (mpq_t RATIONAL, mpz_t DENOMINATOR)
		mpq_class q(0);
		//  q.set_den(g.divisor());
		q.get_den()=g.divisor();
		for (int v = 0; v < num_variables; ++v) {
			//    q.set_num(g.coefficient(Variable(v)));
			q.get_num()=g.coefficient(Variable(v));
			s << " " << q.get_d();
			//cout << "("<<g.coefficient(Variable(v))<<"/"<<g.divisor()<<")";
		}
		s << endl;
	}
}

/*
 void print_generator_fp_raw(ostream& s, const Generator& g)
 {
 const int num_variables = g.space_dimension();
 vector <double> vals(num_variables);
 //  if (vnvec.size()<num_variables)
 //    throw_error("Not enough variable names for constraint.")
 //  bool first = true;

 // for now: only print points and closure points
 // todo: do something about rays and lines
 if (g.is_point() || g.is_closure_point())
 {
 Integer den = g.divisor();
 Integer cv;

 //  double mpq_get_d (mpq_t OP)
 ///  void mpq_div (mpq_t QUOTIENT, mpq_t DIVIDEND, mpq_t
 //          DIVISOR)

 // - Function: void mpq_set_num (mpq_t RATIONAL, mpz_t NUMERATOR)
 // - Function: void mpq_set_den (mpq_t RATIONAL, mpz_t DENOMINATOR)
 mpq_class q(0);
 //  q.set_den(g.divisor());
 q.get_den()=g.divisor();
 for (int v = 0; v < num_variables; ++v)
 {
 //    q.set_num(g.coefficient(Variable(v)));
 q.get_num()=g.coefficient(Variable(v));
 vals[v]=q.get_d();
 }

 if (num_variables!=3)
 {
 for (int v = 0; v < num_variables; ++v)
 {
 s << " " << vals[v];
 }
 }
 else
 {
 s << " " << vals[0]+vals[2]*0.866;
 s << " " << vals[1]+vals[2]*0.5;
 }

 s << endl;
 }
 }
 */

void add_ccvs_to_double_point_list(const Parma_Polyhedra_Library::NNC_Polyhedron& ccvs,double_point_list& pl)
{
  // add the vertices of ccvs to pl
  // ATTENTION: pl is not cleared, points will be added!

 Generator_List gen_list;
 get_and_add_generators(ccvs,gen_list);
//for (Generator_List::const_iterator i=gen_list.begin();i!=gen_list.end();++i) cout << *i <<";";
 add_Generator_List_to_double_point_list(gen_list,pl);
}

void double_point_list_to_ccvs(double_point_list& pl, Parma_Polyhedra_Library::NNC_Polyhedron& ccvs, unsigned int dim)
{
  // set ccvs to the convex hull of the points in pl

  ccvs=Parma_Polyhedra_Library::NNC_Polyhedron(dim,EMPTY);
  Generator g=point(0*Variable(0),1); // dummy
  //double_point p(dim);

  for (double_point_list::iterator i=pl.begin();i!=pl.end();++i)
  {
    double_point_to_generator(*i,g);
    ccvs.add_generator(g);
  }
}

void add_ccvs_consys_to_double_point_list(const Parma_Polyhedra_Library::NNC_Polyhedron& ccvs,double_point_list& pl)
{
  // add the Constraints of ccvs to pl
  // ATTENTION: pl is not cleared, points will be added!
  // the points are of the form [a_1 ... a_n b]
  double_point dx_zero(ccvs.space_dimension()+1);
	for (unsigned int i=0; i<dx_zero.size(); ++i) {  // initialize p to zero
		dx_zero[i]=0.0;
	}

	Constraint_System cs=ccvs.minimized_constraints();
	double_point dx(ccvs.space_dimension()+1);

	for (Constraint_System::const_iterator i=cs.begin();i!=cs.end();++i)
	{
		constraint_to_double_point(*i, dx);
		pl.push_back(dx);
	  if (i->is_equality())
	    {
	      // also add complement
	      pl.push_back(-dx);
	    }
	}
}

void get_and_add_generators(const Parma_Polyhedra_Library::NNC_Polyhedron& ccvs, Generator_List& l)
{
Generator_System gs;
gs=ccvs.minimized_generators();
for (Generator_System::const_iterator i=gs.begin();i!=gs.end();++i)
{
  l.push_back(*i);
}
}

}
