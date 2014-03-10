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
#ifndef GUARD_fp_interface_h
#define GUARD_fp_interface_h

//#include <vector>
//#include <list>
//#include <functional>
//#include <algorithm>
//#include <map>
//#include <stdlib.h>
//#include <stdexcept>
//
#include <ppl.hh>
//#include <time.h>

#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include "core/continuous/polyhedra/ppl_polyhedron/general.h"
//#include "convex_clock_val_set.h"
//#include "clock_val_set.h"
#include "core/continuous/polyhedra/ppl_polyhedron/extended_ppl.h"

#include "math/scalar_types/rational.h"

namespace ppl_polyhedron {

typedef boost::numeric::ublas::vector<double> double_point;
typedef std::list <double_point> double_point_list;

double Round(double Zahl, int Stellen);

// Convex Sorting

template <class T> bool is_counterclockwise_order(const boost::numeric::ublas::vector<T> &A,
		const boost::numeric::ublas::vector<T> &B, const boost::numeric::ublas::vector<T> &C) {
	// Returns whether the points given by 2-dimensional arrays A,B,C in counterclockwise order
	if (A.size()==2 && B.size()==2 && C.size()==2) {
		//		return (C[1]-A[1])*(B[0]-A[0])>=(B[1]-A[1])*(C[0]-A[0]);
		return A[0]*B[1]-A[1]*B[0]+A[1]*C[0]-A[0]*C[1]+B[0]*C[1]-C[0]*B[1]>=0;
	} else
		throw std::runtime_error("wrong dimension in is_counterclockwise_order");
	return false;
};

template <class T> bool compare_second_coord(const boost::numeric::ublas::vector <T> &p1,
		const boost::numeric::ublas::vector <T> &p2) {
	return p1[1]<=p2[1];
};

template <class T> void sort_second_coord(std::list <boost::numeric::ublas::vector <T> > &l) {
	if (l.size()>1) {
		typename std::list <boost::numeric::ublas::vector <T> >::iterator i=l.begin();
		typename std::list <boost::numeric::ublas::vector <T> >::iterator j=l.begin();
		typename std::list <boost::numeric::ublas::vector <T> >::iterator m=l.end();
		--m;
		boost::numeric::ublas::vector <T> dum;

		/*  int m, i;
		 for (m = a.length - 1; m > 0; m--) {  // counting down
		 for (i = 0; i < m; i++) {       // bubbling up
		 if (a[i] > a[i + 1]) { // if out of order...
		 int temp = a[i];          // ...then swap
		 a[i] = a[i + 1];
		 a[i + 1] = temp;
		 }
		 }
		 } */

		/*  for (m = (array_size - 1); m >= 0; m--)
		 {
		 for (i = 0; i < m; i++)
		 {
		 if (numbers[i] > numbers[i+1])
		 {
		 temp = numbers[i];
		 numbers[i] = numbers[i+1];
		 numbers[i+1] = temp;
		 }
		 }
		 }*/

		while (m!=l.begin()) {
			for (i=l.begin(); i!=m; ++i) {
				j=i;
				++j;
				if (!compare_second_coord(*i, *j)) {
					// swap
					dum=(*i);
					*i=(*j);
					*j=dum;
				}
			}
			--m;
		}

		// Eliminate doubles
		i=l.begin();
		while (i!=l.end()) {
			j=i;
			while (j!=l.end()) {
				if (i!=j && *i==*j) {
					j=l.erase(j);
				} else
					++j;
			}
			++i;
		}

		//		for (i=l.begin(); i!=l.end(); ++i)
		//			cout << (*i)[0] <<"," << (*i)[1] << endl;
		//		cout << ";" << endl;
	}
};

template <class T> void sort_counterclockwise(std::list <boost::numeric::ublas::vector <T> > &l) {
	// ATTENTION: there must not be doubles!
	if (l.size()>2) {
		sort_second_coord(l);
		typename std::list <boost::numeric::ublas::vector <T> >::iterator i=l.begin();
		typename std::list <boost::numeric::ublas::vector <T> >::iterator i2=l.begin();
		typename std::list <boost::numeric::ublas::vector <T> >::iterator j=l.begin();
		++j;
		typename std::list <boost::numeric::ublas::vector <T> >::iterator m=l.end();
		--m;
		--m;
		boost::numeric::ublas::vector <T> dum;

		while (m!=l.begin()) {
			for (i=l.begin(); i!=m; ++i) {
				j=i;
				++j;
				i2=j;
				++i2;
				if (!is_counterclockwise_order(*i, *j, *i2)) {
					// swap
					dum=(*i);
					*i=(*j);
					*j=dum;
				}
			}
			--m;
		}
	}
};

/*
 // preprocess so that p[1] has smallest y-coordinate
 // sort by angle with p[1]
 points[0] = points[N]; // sentinel
 int M = 3;
 for (int i = 4; i <= N; i++) {
 while (Point.ccw(p[M], p[M-1], p[i]) >= 0) {
 M--; // back up to include i on hull
 }
 M++;
 swap(points, M, i); // add i to putative hull
 }
 */

template <class T> void sort_counterclockwise2(std::list <boost::numeric::ublas::vector <T> > &l) {
	if (l.size()>=3) {
		//		l.sort(compare_first_coord);
		//		l.sort();
		sort_second_coord(l);

		typename std::list <boost::numeric::ublas::vector <T> >::iterator i=l.end();
		l.push_front(*(--i)); // sentinel

		typename std::list <boost::numeric::ublas::vector <T> >::iterator m=++(++(++(l.begin()))); // M=3
		typename std::list <boost::numeric::ublas::vector <T> >::iterator m1;

		if (l.begin()->dim()==2) {
			boost::numeric::ublas::vector <T> dum;

			for (i=++(++(++(++(l.begin())))); i!=l.end(); ++i) {
				m1=m;
				--m1;
				while (is_counterclockwise_order(*m, *m1, *i)) {
					--m;
				}

				++m;

				// swap i and m
				dum=(*i).copy();
				*i=(*m).copy();
				*m=dum.copy();
			}
		}
	}
};

/*
 template <class T>
 void sort_counterclockwise(list <boost::numeric::ublas::vector <T> > &l)
 {
 // Todo: This is an inefficient implementation.
 if (l.size()>=3)
 {
 if (l.begin()->dim()==2)
 {
 typename list <boost::numeric::ublas::vector <T> >::iterator j,m;
 typename list <boost::numeric::ublas::vector <T> >::iterator i=l.begin();
 j=i;
 ++j;
 m=j;
 ++m;
 boost::numeric::ublas::vector <T> dum;

 while (i!=l.end())
 {
 if (!is_counterclockwise_order(*i,*j,*m))
 {
 // swap j and m
 dum=(*j).copy();
 *j=(*m).copy();
 *m=dum.copy();
 // decrease i so that i and m are compared next
 //					i=l.begin();
 if (i!=l.begin())
 --i;
 else
 ++i;
 if (m==l.begin() || j==l.begin()) // start over if end was changed
 i=l.begin();
 }
 else
 ++i;
 j=i;
 ++j;
 if (j==l.end())
 {
 j=l.begin();
 }
 m=j;
 ++m;
 if (m==l.end())
 {
 m=l.begin();
 }
 }
 }
 }
 }
 */

template <class T> void print_fp_raw(std::ostream& o, const boost::numeric::ublas::vector<T> &A) {
	for (unsigned int i=0; i<A.size(); ++i)
	{
		if (i>0) o << " ";
		o << A[i];
	}
	o << std::endl;
};

template <class T> void print_fp_raw(std::ostream& o, std::list<boost::numeric::ublas::vector<T> > &l) {
	// if it is two-dimensional, sort it so it's counter-clockwise
	if (!l.empty() && l.begin()->size()==2) {
		sort_counterclockwise(l);
		//		sort_counterclockwise(l);

		// close the curve
		if (l.size()>2)
			l.push_back(*l.begin());
	}
	for (typename std::list <boost::numeric::ublas::vector <T> >::const_iterator i=l.begin(); i
			!=l.end(); ++i)
		print_fp_raw(o, *i);
};

// Matrix Operations

template <class T> T scalar_product(const boost::numeric::ublas::vector<T> &A,
		const boost::numeric::ublas::vector<T> &B) {
	// Scalar Prodcut of 2 Vectors, i.e., sum_i(A[i]*B[i])
	int n = A.size();
	int m = B.size();

	if (n != m || n<=0)
		return T(0);

	else {
		T c;

		c=A[0]*B[0];
		for (int i=1; i<n; i++) {
			c+=A[i]*B[i];
		}
		return c;
	}
};

template <class T> boost::numeric::ublas::vector<T> operator*(const T &B,
		const boost::numeric::ublas::vector<T> &A) {
	// Scalar*Vector

	int n = A.size();
	boost::numeric::ublas::vector<T> C(n);

	for (int i=0; i<n; i++) {
		C[i] = B*A[i];
	}
	return C;
};

template <class T> boost::numeric::ublas::vector<T> operator*(const boost::numeric::ublas::matrix<T> &A,
		const boost::numeric::ublas::vector<T> &B) {
	// Matrix * Vector
	int n = A.size1();
	int m = A.size2();

	if (B.size() != m || n<=0 || m<=0)
		return boost::numeric::ublas::vector<T>();

		else
		{
			boost::numeric::ublas::vector<T> C(n);

			for (int i=0; i<n; i++)
			{
				C[i]=A[i][0]*B[0];
				for (int j=1; j<m; j++)
				C[i] += A[i][j] * B[j];
			}
			return C;
		}
	};

template <class T> bool operator==(const boost::numeric::ublas::vector<T> &A,
		const boost::numeric::ublas::vector<T> &B) {
	// Scalar Prodcut of 2 Vectors, i.e., sum_i(A[i]*B[i])
	int n = A.size();
	int m = B.size();

	if (n != m || n<=0)
		return false;

	else {
		for (int i=0; i<n; i++) {
			if (!(A[i]==B[i])) {
				return false;
			}
		}
		return true;
	}
};

template <class T> T norm2(const boost::numeric::ublas::vector<T> &A) {
	// return |A^T*A|^2=A^T*A
	int n = A.size();

	T c=0;

	if (n>0) {
		for (int i=0; i<n; i++) {
			c += A[i]*A[i];
		}
	}
	return c;
};

// Conversion Methods

double generator_to_double(Parma_Polyhedra_Library::Generator& g, Parma_Polyhedra_Library::dimension_type pos);

void generator_to_double_point(const Parma_Polyhedra_Library::Generator& g, double_point& p);

void get_max_coeff(const Parma_Polyhedra_Library::Linear_Expression& g, Integer& d);

void get_max_coeff(const Parma_Polyhedra_Library::Constraint& g, Integer& d);

void hom_linexpression_to_double_point(const Parma_Polyhedra_Library::Linear_Expression& g,
		double_point& p);

void double_point_to_Linear_Expression(const double_point& p,
		Parma_Polyhedra_Library::Linear_Expression& lin, Integer& d);

void double_point_to_generator(const double_point& p, Parma_Polyhedra_Library::Generator& g);

void
		constraint_to_double_point(const Parma_Polyhedra_Library::Constraint& c, double_point& dx,
				double& b);

void constraint_to_double_point(const Parma_Polyhedra_Library::Constraint& c, double_point& dx);

void add_Generator_List_to_double_point_list(const Generator_List& gl,
		double_point_list& pl);

// Computation Methods

void get_constraint_through_double_point(Parma_Polyhedra_Library::Linear_Expression linex,
		double_point& p, Parma_Polyhedra_Library::Constraint& c, bool strict=false);

void get_constraint_through_double_point(double_point& dx, double_point& p,
		Parma_Polyhedra_Library::Constraint& c, bool strict=false);

void get_constraint_through_double(double_point& dx, double q, Parma_Polyhedra_Library::Constraint& c,
		bool strict=false);

void get_double_point_list_center(double_point_list& pl, double_point& p);

//void get_double_point_list_min_max(const double_point_list& pl,const double_point& dx,double& min_dx,double& max_dx);
void get_double_point_list_min_max(const double_point_list& pl,
		const double_point& dx, double& min_dx, double& max_dx,
		double_point& x_min, double_point& x_max);

double get_double_point_list_angle(double_point_list& pl);

double get_double_point_list_angle_vecs(double_point_list& pl,
		double_point& g1, double_point& g2);

double get_double_point_list_variance(double_point_list& pl);

void print_generator_fp_raw(std::ostream& s, const Parma_Polyhedra_Library::Generator& g);


void add_ccvs_to_double_point_list(const Parma_Polyhedra_Library::NNC_Polyhedron& ccvs,double_point_list& pl);

void double_point_list_to_ccvs(double_point_list& pl, Parma_Polyhedra_Library::NNC_Polyhedron& ccvs, unsigned int dim);

void add_ccvs_consys_to_double_point_list(const Parma_Polyhedra_Library::NNC_Polyhedron& ccvs,double_point_list& pl);

void get_and_add_generators(const Parma_Polyhedra_Library::NNC_Polyhedron& ccvs, Generator_List& l);

}

#endif
