/*
 * piecewise_linear_function.h
 *
 *  Created on: Oct 18, 2012
 *      Author: kateja
 */

#ifndef PIECEWISE_LINEAR_FUNCTION_H_
#define PIECEWISE_LINEAR_FUNCTION_H_

#include <iostream>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "interval.h"
#include "breakpoint.h"

namespace plif {

/**
 * These two are symbols to be used for positive and negative infinity, using the definitions in C++ standard library
 */
#define NEG_INFTY (-1*std::numeric_limits<precision_type>::infinity())
#define POS_INFTY (std::numeric_limits<precision_type>::infinity())

/**
 * Class for a piecewise Linear funcion.
 *
 * Data Memebers are domain, list of points, boolean left and right unbounded, and left and right slope.
 * Left and Right slope are meaningful only when function is left unbounded and right unbounded respectively.
 * On bounded domains, the end points will not have valid left and right limits. Have to implement that check.
 */
class piecewise_linear_function
{
public:
	typedef std::vector<breakpoint> list_of_points_type;
private:
	interval domain;						/** Domain of the PLF */
	list_of_points_type list_of_points;		/** Breakpoint list */
	bool left_unbounded;					/** Is the function unbounded on the left */
	bool right_unbounded;					/** Is the function unbounded on the right */
	bool left_closed;				/** Is the leftmost point in the list of breakpoints a part of the domain, is the domain left closed*/
	bool right_closed;				/** Is the rightmost point in the list of breakpoints a part of the domain, is the domain right closed*/
	precision_type left_slope;				/** If left_unbounded, left_slope value */
	precision_type right_slope;				/** If left_unbounded, left_slope value */
	computation_parameters comp_params; /** Computation parameters like rounding */
public:
	/** Construct a closed piecewise linear function from a list of breakpoints
	 *
	 * If no list is passed, there are no breakpoints and the domain is the empty interval. */
	static piecewise_linear_function create_closed(const list_of_points_type& breakpts=list_of_points_type(), const computation_parameters& params=computation_parameters()) {
		piecewise_linear_function f;
		f.set_list(breakpts);
		f.set_right_closed(true);
		f.set_left_closed(true);
		f.set_left_bounded();
		f.set_right_bounded();
		f.set_computation_parameters(params);

		// find the interval bounds
		interval dom=interval::empty(); // start with the empty interval
		if (breakpts.size()>0) {
			// sort the breakpoints
			std::sort(f.list_of_points.begin(),f.list_of_points.end());
// GF: do the following if you don't want to sort the list
//			precision_type xmin=*std::min_element(f.list_of_points.begin(),f.list_of_points.end());
//			precision_type xmax=*std::max_element(f.list_of_points.begin(),f.list_of_points.end());
			precision_type xmin=f.list_of_points.begin()->get_x(); // first element in sorted list
			precision_type xmax=f.list_of_points.rbegin()->get_x(); // last element in sorted list
			dom = interval(xmin,xmax);
		}
		f.set_domain(dom);

		return f;
	};

	/** Construct a closed piecewise linear function from a list of breakpoints, adopting the list
	 *
	 * Assumes that the list is sorted. */
	static piecewise_linear_function create_closed_from_temp(list_of_points_type& breakpts, const computation_parameters& params=computation_parameters()) {
		piecewise_linear_function f;
		f.list_of_points.swap(breakpts);
		f.set_right_closed(true);
		f.set_left_closed(true);
		f.set_left_bounded();
		f.set_right_bounded();
		f.set_computation_parameters(params);

		// find the interval bounds
		interval dom=interval::empty(); // start with the empty interval
		if (f.list_of_points.size()>0) {
			precision_type xmin=f.list_of_points.begin()->get_x(); // first element in sorted list
			precision_type xmax=f.list_of_points.rbegin()->get_x(); // last element in sorted list
			dom = interval(xmin,xmax);
		}
		f.set_domain(dom);

		return f;
	};
	/** Get the number of breakpoints
	 */
	size_t size() const {
		return list_of_points.size();
	}
	/**
	 *Functions to set and get values of the memeber data fields
	 */
	void set_domain(interval I)
		{domain = I;}
	const interval& get_domain() const
		{return domain;}
	void set_left_unbounded(precision_type slope)
	{
		left_slope = slope;
		left_closed = true;
		left_unbounded = true;
	}
	void set_left_bounded()
		{left_unbounded = false;}
	const bool& is_left_unbounded() const
		{return left_unbounded;}
	void set_right_unbounded(precision_type slope)
	{
		right_slope = slope;
		right_closed = true;
		right_unbounded = true;
	}
	void set_right_bounded()
		{right_unbounded = false;}
	bool is_right_unbounded() const
		{return right_unbounded;}
	const precision_type& get_left_slope() const
	{
		assert(left_unbounded);
		return left_slope;
	}
	const precision_type& get_right_slope() const
	{
		assert(right_unbounded);
		return right_slope;
	}
	void set_left_closed(bool b)
		{left_closed = b;}
	const bool& is_left_closed() const
		{return left_closed;}
	void set_right_closed(bool b)
		{right_closed = b;}
	const bool& is_right_closed() const
		{return right_closed;}
	/** Set computation parameters */
	void set_computation_parameters(const computation_parameters& params) {
		comp_params = params;
	};
	/** Returns computation parameters */
	const computation_parameters& get_computation_parameters() const {
		return comp_params;
	};
	/** Sets computations to round up */
	void set_rounding_up() {
		comp_params.rounding = computation_parameters::up;
	};
	/** Sets computations to round up */
	void set_rounding_down() {
		comp_params.rounding = computation_parameters::down;
	};
	void insert_point(breakpoint point_to_insert)
	{
		list_of_points_type::iterator iter;
		bool present = false;
		for(iter=list_of_points.begin();iter!=list_of_points.end();iter++)
		{
			if((*iter).get_x() == point_to_insert.get_x())
			{
				present = true;
				break;
			}
			if((*iter).get_x() > point_to_insert.get_x())
			{
				list_of_points.insert(iter, point_to_insert);
				break;
			}
		}
		if(!present && iter == list_of_points.end())
			list_of_points.push_back(point_to_insert);

		// adapt domain to include point_to_insert
		if (empty(domain)) {
			domain = interval(point_to_insert.get_x(), point_to_insert.get_x());
		} else {
			if (!is_left_unbounded() && domain.lower()
					> point_to_insert.get_x())
				domain = interval(point_to_insert.get_x(), domain.upper());
			if (!is_right_unbounded() && domain.upper()
					< point_to_insert.get_x())
				domain = interval(domain.lower(), point_to_insert.get_x());
		}
	}
	/** @author Goran Frehse */
	void push_back(breakpoint point_to_insert)
	{
		list_of_points.push_back(point_to_insert);

		// adapt domain to include point_to_insert
		if (empty(domain)) {
			domain = interval(point_to_insert.get_x(), point_to_insert.get_x());
		} else {
			if (!is_right_unbounded() && domain.upper()
					< point_to_insert.get_x())
				domain = interval(domain.lower(), point_to_insert.get_x());
		}
	}
	void insert_point(precision_type a, precision_type b, precision_type c, precision_type d)
	{
		breakpoint p(a, b, c, d);
		insert_point(p);
	}
	/** Insert a point (x,y) at which the function is continuous. */
	void insert_continuous_point(precision_type x, precision_type y)
	{
		breakpoint p(x, y);
		insert_point(p);
	}
	/** Sets the list of breakpoints. */
	void set_list(const list_of_points_type& list) {
		list_of_points = list;
	}
	/** Swaps the breakpoints of *this with list */
	void set_list_from_temp(list_of_points_type& list) {
		list_of_points.swap(list);
	}
	/** Returns a const reference to the list of breakpoints. */
	const list_of_points_type& get_list() const {
		return list_of_points;
	}
	/** Returns a const reference to the list of breakpoints. */
	list_of_points_type& get_list() {
		return list_of_points;
	}
	/**
	 * Displays the function fields, to be used mainly for debugging purposes.
	 *
	 * If no stream is passed as an argument, the output is sent to std::cout.
	 */
	void display(std::ostream& cout=std::cout) const
	{
		cout << "Domain: [" << domain.lower() << ", " << domain.upper() << "]" << std::endl;
		cout << "left_unbounded " << left_unbounded << std::endl;
		cout << "right_unbounded " << right_unbounded << std::endl;
		cout << "Left Slope: " << left_slope << std::endl;
		cout << "Right Slope: " << right_slope << std::endl;
		cout << "Left Closed: " << left_closed << std::endl;
		cout << "Right Closed: " << right_closed << std::endl;
		cout << "List of points" << std::endl;
		int i = 0;
		for(i=0;i<list_of_points.size();i++)
		{
			cout <<	"\t" << list_of_points[i].get_x() << "  " << list_of_points[i].get_y_left() << "  " << list_of_points[i].get_y() << "  " << list_of_points[i].get_y_right() << std::endl;
		}

	}
	/**
	 *Computes the value of the function at the provided x
	 */
	precision_type value_at(precision_type x) const
	{
		breakpoint b = at(x);
		return b.get_y();
	}
	/**
	 *Computes the left limit of the function at the provided x. Implementation is same as that of value at function
	 */
	precision_type left_lim_at(precision_type x) const
	{
		breakpoint b = at(x);
		return b.get_y_left();
	}
	/**
	 *Computes the right limit of the function at the provided x. Implementation is same as that of value at function
	 */
	precision_type right_lim_at(precision_type x) const
	{
		breakpoint b = at(x);
		return b.get_y_right();
	}
	/** Returns a breakpoint at position x
	 *
	 * If necessary, this is an interpolation of the surrounding points. */
	breakpoint at(const precision_type& x) const {
		const breakpoint& leftmost = *list_of_points.begin();
		const breakpoint& rightmost = *list_of_points.rbegin();
		if(x < leftmost.get_x()) {
			if (left_unbounded) {
				return breakpoint(x,leftmost.get_y_left() + left_slope*(leftmost.get_x()-x));
			}
			else {
				throw std::runtime_error("x below domain");
			}
		}
		if(x > rightmost.get_x()) {
			if (left_unbounded) {
				return breakpoint(x,rightmost.get_y_right() + right_slope*(x-rightmost.get_x()));
			}
			else {
				throw std::runtime_error("x above domain");
			}
		}

		const list_of_points_type& list=get_list();
		breakpoint b(x,precision_type(0));
		// first that is not less than x
		list_of_points_type::const_iterator bit=std::lower_bound(list.begin(),list.end(),b);
		// bit->x >= x
		if (bit->get_x() == x) {
			return *bit;
		} else {
			// bit->x > x
			// first that is greater than x
			list_of_points_type::const_iterator eit=bit;
			// move bit one back so its x is less than x
			--bit;
			// bit->x < x
			return interpolate(*bit,*eit,x,get_computation_parameters());
		}
	}
	/**
	 * Copies from the PLF f to this PLF on the domain identified by (star, end).
	 * left_closed and right_closed are to be set outside this function.
	 * same with left_unbouned, right_unbounded, left_slope and right_slope if applicable
	 * so it basically copies only the list of breakpoints in the domain
	 *
	 * @attention requires the breakpoints to be sorted
	 *
	 * @author complete reimplementation by Goran Frehse
	 */
	void copy_on_domain(const piecewise_linear_function& f, precision_type start, precision_type end)
	{
		const list_of_points_type& list=f.get_list();
		domain.set(start, end);
		// find lowest element in domain
		breakpoint bsta(start,precision_type(0));
		breakpoint bend(end,precision_type(0));
		// first one >=
		list_of_points_type::const_iterator bit=std::lower_bound(list.begin(),list.end(),bsta);
		// first one >
		list_of_points_type::const_iterator eit=std::upper_bound(bit,list.end(),bend);
		// point to the last element in the range
		list_of_points_type::const_iterator rit=eit;
		if (rit!=list.begin()) {
			--rit;
		}
		// need additional start and end points if not exactly on bit and eit
		unsigned int s = eit-bit;
		bool insert_first = bit == list.end() || bit->get_x()>start;
		bool insert_last = rit == list.end() || rit->get_x()>end || rit->get_x()<end;
		if (insert_first) {
			++s;
		}
		if (insert_last) {
			++s;
		}
		list_of_points=list_of_points_type(s);
		list_of_points_type::iterator it = list_of_points.begin();
		if (insert_first) {
			*it = f.at(start);
			++it;
		}
		// copy the breakpoints from f
		it = std::copy(bit,eit,it);
		if (insert_last) {
			*it = f.at(end);
			++it;
		}
		// check
		if (it!=list_of_points.end()) {
			std::cerr << "orig: ";
			f.display(std::cerr);
			std::cerr << "copied: ";
			display(std::cerr);
			std::cerr << insert_first << insert_last << "s: " << s << ", it: " << it-bit << ", copied: " << eit-bit+1 << std::endl;
			throw std::runtime_error("messed up copy_on_domain");
		}
	}
	/** Equality check
	 *
	 * Does not take into account left_slope and right_slope if the domain is bounded on that side.
	 * */
	bool operator== (const piecewise_linear_function& g) const
	{
		/** Note that comparing intervals for equality does not have the usual meaning with boost!
		 * So we compare the interval bounds here explicitly.
		 */
		return ((empty(domain) && empty(g.domain)) ||
					(domain.lower()==g.domain.lower() && domain.upper()==g.domain.upper())) &&
				list_of_points==g.list_of_points
				&& (!left_unbounded || (left_slope==g.left_slope && left_closed==g.left_closed ))
				&& (!right_unbounded || (right_slope==g.right_slope && right_closed==g.right_closed ));
	}

	/** Scalar Multiplication */
	piecewise_linear_function& operator*=(precision_type c) {
		if (is_left_unbounded()) {
			left_slope *= c;
		}
		if (is_right_unbounded()) {
			right_slope *= c;
		}
		// multiply breakpoints
		for (list_of_points_type::iterator it = list_of_points.begin();
				it != list_of_points.end(); ++it) {
			it->operator *=(c);
		}
		return *this;
	}

	/** Scalar Addition */
	piecewise_linear_function& operator+=(precision_type c) {
		// add to breakpoints
		for (list_of_points_type::iterator it = list_of_points.begin();
				it != list_of_points.end(); ++it) {
			it->operator+=(c);
		}
		return *this;
	}


	friend piecewise_linear_function operator*(const piecewise_linear_function& f, const precision_type& lambda);
};

}

#endif /* PIECEWISE_LINEAR_FUNCTION_H_ */
