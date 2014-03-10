/*
 * piecewise_linear_function_paths.h
 *
 *  Created on: Oct 18, 2012
 *      Author: kateja
 */

#ifndef PIECEWISE_LINEAR_FUNCTION_PATHS_H_
#define PIECEWISE_LINEAR_FUNCTION_PATHS_H_

#include "piecewise_linear_function_operators.h"

#include <list>
#include <deque>
#include <sstream>

#include "interval.h"
#include "piecewise_linear_interval_function.h"
#include "piecewise_linear_interval_function_utility.h"

namespace plif {

/** 
 * Checks if B is visible to A, i.e. is there a straight line joining A to B lying above LOWER and below UPPER
 * for now, in the case when a.x() == b.x(), the algo might run into bugs
 * @todo: fix the case when x1 == x2, require multiple checks for values since left limits and right limits are to be considered and x1 may be
 *	 	 an endpoint. for now, an assertion x1!=x2 has been placed.
 */
inline
bool is_visible(piecewise_linear_function lower, piecewise_linear_function upper, breakpoint a, breakpoint b)
{
	
	const std::vector<breakpoint>& list_lower = lower.get_list();
	const std::vector<breakpoint>& list_upper = upper.get_list();
	int i;
	precision_type x1;
	precision_type y1;
	precision_type x2;
	precision_type y2;
	if(a.get_x() < b.get_x())
	{
		x1 = a.get_x();
		y1 = a.get_y();
		x2 = b.get_x();
		y2 = b.get_y();
	}
	else 
	{
		x1 = b.get_x();
		y1 = b.get_y();
		x2 = a.get_x();
		y2 = a.get_y();
	}
	assert(x1!=x2);
	
	if(x1!=lower.get_domain().lower() || (x1==lower.get_domain().lower() && lower.is_left_closed()))
	{
		if(lower.value_at(x1) > y1)
			return false;
	}		
	if(lower.right_lim_at(x1) > y1)
		return false;	
	for(i=0;i<list_lower.size();i++)
	{
		if(list_lower[i].get_x() <= x1)
			continue;
		else if(list_lower[i].get_x() >= x2)
			break;
		else
		{	
			precision_type x = list_lower[i].get_x();
			precision_type y;
			if(i==0)
			{
				y = list_lower[i].get_y_right();
				if(y > y1+(y2-y1)/(x2-x1)*(x-x1))
					return false;
				if(lower.is_left_closed())
				{
					y = list_lower[i].get_y();
					if(y > y1+(y2-y1)/(x2-x1)*(x-x1))
						return false;
				}	
			}
			else if(i==list_lower.size()-1)
			{
				y = list_lower[i].get_y_left();
				if(y > y1+(y2-y1)/(x2-x1)*(x-x1))
					return false;
				if(lower.is_right_closed())
				{
					y = list_lower[i].get_y();
					if(y > y1+(y2-y1)/(x2-x1)*(x-x1))
						return false;
				}	
			}
			else
			{
				y = list_lower[i].get_y_left();	
				if(y > y1+(y2-y1)/(x2-x1)*(x-x1))
						return false;
				y = list_lower[i].get_y();
				if(y > y1+(y2-y1)/(x2-x1)*(x-x1))
					return false;
				y = list_lower[i].get_y_right();
				if(y > y1+(y2-y1)/(x2-x1)*(x-x1))
					return false;		
			}		
		}		
	}
	if(x2!=lower.get_domain().upper() || (x2==lower.get_domain().upper() && lower.is_right_closed()))
	{
		if(lower.value_at(x2) > y2)
			return false;
	}
	if(lower.left_lim_at(x2) > y2)
		return false;	


	if(x1!=upper.get_domain().lower() || (x1==upper.get_domain().lower() && upper.is_left_closed()))
	{
		if(upper.value_at(x1) < y1)
			return false;
	}		
	if(upper.right_lim_at(x1) < y1)
		return false;	
	for(i=0;i<list_upper.size();i++)
	{
		if(list_upper[i].get_x() <= x1)
			continue;
		else if(list_upper[i].get_x() >= x2)
			break;	
		else
		{	
			precision_type x = list_upper[i].get_x();
			precision_type y;
			if(i==0)
			{
				y = list_upper[i].get_y_right();
				if(y < y1+(y2-y1)/(x2-x1)*(x-x1))
					return false;
				if(upper.is_left_closed())
				{
					y = list_upper[i].get_y();
					if(y < y1+(y2-y1)/(x2-x1)*(x-x1))
						return false;
				}	
			}
			else if(i==list_upper.size()-1)
			{
				y = list_upper[i].get_y_left();
				if(y < y1+(y2-y1)/(x2-x1)*(x-x1))
					return false;
				if(upper.is_right_closed())
				{
					y = list_upper[i].get_y();
					if(y < y1+(y2-y1)/(x2-x1)*(x-x1))
						return false;
				}
			}
			else
			{
				y = list_upper[i].get_y_left();
				if(y < y1+(y2-y1)/(x2-x1)*(x-x1))
					return false;
				y = list_upper[i].get_y();
				if(y < y1+(y2-y1)/(x2-x1)*(x-x1))
					return false;
				y = list_upper[i].get_y_right();
				if(y < y1+(y2-y1)/(x2-x1)*(x-x1))
					return false;		
			}		
		}		
	}
	if(x2!=upper.get_domain().upper() || (x2==upper.get_domain().upper() && upper.is_right_closed()))
	{
		if(upper.value_at(x2) < y2)
			return false;
	}		
	if(upper.left_lim_at(x2) < y2)
		return false;	
	
	return true;
}

/**
 * Finds the Eucledian shortest path from FROM to TO, which must be breakpoints of either of the functions
 * the path lies in between the LOWER and UPPER PLF
 * The PLF's must be well behaved in the sense that left limits, value and right limits at the breakpoints must match
 * FROM and TO must be well defined as well (in the same sense as above)
 *
 * If min_link = true, the minimum link path is constructed that consists of vertices. Otherwise, the Euclidean distance
 * is minimized.
 *
 * @author Rajat Kateja, with modifications by Goran Frehse
 */ 
inline
piecewise_linear_function shortest_path(piecewise_linear_function lower, piecewise_linear_function upper, breakpoint from, breakpoint to, bool min_link = false)
{
	assert(dist(lower, upper) >= (precision_type)0);
	const std::vector<breakpoint>& list_lower = lower.get_list();
	const std::vector<breakpoint>& list_upper = upper.get_list();
	std::vector<breakpoint> list;
	list.clear();
	int lower_size = list_lower.size();
	int upper_size = list_upper.size();
	int i, j;
	int source, target;
	for(i=0;i<lower_size;i++)
	{
		breakpoint p;
		if(i==0)
		{
			p.set_x(list_lower[i].get_x());
			p.set_y_left(list_lower[i].get_y_right());
			p.set_y(list_lower[i].get_y_right());
			p.set_y_right(list_lower[i].get_y_right());
			list.push_back(p);
			if(lower.is_left_closed() && list_lower[i].get_y() != list_lower[i].get_y_right())
			{
				p.set_y_left(list_lower[i].get_y());
				p.set_y(list_lower[i].get_y());
				p.set_y_right(list_lower[i].get_y());
				list.push_back(p);
			}
		}
		else if(i==lower_size-1)
		{
			p.set_x(list_lower[i].get_x());
			p.set_y_left(list_lower[i].get_y_left());
			p.set_y(list_lower[i].get_y_left());
			p.set_y_right(list_lower[i].get_y_left());
			list.push_back(p);
			if(lower.is_right_closed() && list_lower[i].get_y() != list_lower[i].get_y_left())
			{
				p.set_y_left(list_lower[i].get_y());
				p.set_y(list_lower[i].get_y());
				p.set_y_right(list_lower[i].get_y());
				list.push_back(p);
			}
		}
		else
		{
			p.set_x(list_lower[i].get_x());
			p.set_y_left(list_lower[i].get_y());
			p.set_y(list_lower[i].get_y());
			p.set_y_right(list_lower[i].get_y());
			list.push_back(p);
			if(list_lower[i].get_y_left() != list_lower[i].get_y())
			{
				p.set_y_left(list_lower[i].get_y_left());
				p.set_y(list_lower[i].get_y_left());
				p.set_y_right(list_lower[i].get_y_left());
				list.push_back(p);
			}
			if(list_lower[i].get_y_right() != list_lower[i].get_y())
			{
				p.set_y_left(list_lower[i].get_y_right());
				p.set_y(list_lower[i].get_y_right());
				p.set_y_right(list_lower[i].get_y_right());
				list.push_back(p);
			}
		}
	}
	for(i=0;i<upper_size;i++)
	{
		breakpoint p;
		if(i==0)
		{
			p.set_x(list_upper[i].get_x());
			p.set_y_left(list_upper[i].get_y_right());
			p.set_y(list_upper[i].get_y_right());
			p.set_y_right(list_upper[i].get_y_right());
			list.push_back(p);
			if(upper.is_left_closed() && list_upper[i].get_y() != list_upper[i].get_y_right())
			{
				p.set_y_left(list_upper[i].get_y());
				p.set_y(list_upper[i].get_y());
				p.set_y_right(list_upper[i].get_y());
				list.push_back(p);
			}
		}
		else if(i==upper_size-1)
		{
			p.set_x(list_upper[i].get_x());
			p.set_y_left(list_upper[i].get_y_left());
			p.set_y(list_upper[i].get_y_left());
			p.set_y_right(list_upper[i].get_y_left());
			list.push_back(p);
			if(upper.is_right_closed() && list_upper[i].get_y() != list_upper[i].get_y_left())
			{
				p.set_y_left(list_upper[i].get_y());
				p.set_y(list_upper[i].get_y());
				p.set_y_right(list_upper[i].get_y());
				list.push_back(p);
			}
		}
		else
		{
			p.set_x(list_upper[i].get_x());
			p.set_y_left(list_upper[i].get_y());
			p.set_y(list_upper[i].get_y());
			p.set_y_right(list_upper[i].get_y());
			list.push_back(p);
			if(list_upper[i].get_y_left() != list_upper[i].get_y())
			{
				p.set_y_left(list_upper[i].get_y_left());
				p.set_y(list_upper[i].get_y_left());
				p.set_y_right(list_upper[i].get_y_left());
				list.push_back(p);
			}
			if(list_upper[i].get_y_right() != list_upper[i].get_y())
			{
				p.set_y_left(list_upper[i].get_y_right());
				p.set_y(list_upper[i].get_y_right());
				p.set_y_right(list_upper[i].get_y_right());
				list.push_back(p);
			}
		}
	
	}
	int graph_size = list.size();
	precision_type graph[graph_size][graph_size];
	precision_type distance[graph_size];
	int pred[graph_size];
	bool reached[graph_size];

	for(i=0;i<graph_size;i++)
	{
		for(j=0;j<=i;j++)
		{
			if(i==j)
			{
				graph[i][j] = graph[j][i] = 0;
			}
			else if(list[i].get_x() == list[j].get_x())
			{
				if (min_link) {
					graph[i][j] = graph[j][i] = 1;
				} else {
					graph[i][j] = graph[j][i] = std::max((list[i].get_y() - list[j].get_y()), (list[j].get_y() - list[i].get_y()));
				}
			}
			else if(is_visible(lower, upper, list[i], list[j]))
			{
				if (min_link) {
					graph[i][j] = graph[j][i] = 1;
				} else {
					// Euclidean distance
					graph[i][j] = graph[j][i] = length(list[i], list[j]);
				}
			}
			else
			{
				graph[i][j] = graph[j][i] = -1;
			}
		}
	}
	
	for(i=0;i<graph_size;i++)
	{
		distance[i] = POS_INFTY;
		pred[i] = -1;
		reached[i] = false;
		if(list[i] == from)
		{
			distance[i] = 0;
			source = i;
		}	
		if(list[i] == to)
		{
			target = i;
		}
	}

	int min = 0;
	while(!reached[target])
	{
		precision_type min_dist = POS_INFTY;
		for(i=0;i<graph_size;i++)
		{
			if(distance[i] < min_dist && !reached[i])
			{
				min = i;
				min_dist = distance[i];
			}	
		}

		reached[min] = true;

		if(distance[min] == POS_INFTY)
			break;
			

		for(i=0;i<graph_size;i++)
		{
			if(graph[min][i] != -1)
			{
				precision_type dist = min_dist + graph[min][i];
				if(dist < distance[i])
				{
					distance[i] = dist;
					pred[i] = min;
				}
			}

		}
	}

	piecewise_linear_function shortest_path;

	i = target;
	while(pred[i]!=-1)
	{

		shortest_path.insert_point(list[i]);
		i = pred[i];
	}
	shortest_path.insert_point(list[source]);
	std::vector<breakpoint>& shortest_list = shortest_path.get_list();
	std::vector<breakpoint>::iterator iter;

	for(iter=shortest_list.begin();iter!=shortest_list.end();)
	{
		if((iter+1)!=shortest_list.end() && (*iter).get_x() == (*(iter+1)).get_x())
		{
			(*iter).set_y_right((*(iter+1)).get_y());
			iter = shortest_list.erase(iter+1);
		}
		else
			iter = iter + 1;
	}
	shortest_path.set_domain(interval(from.get_x(), to.get_x()));
	shortest_path.set_left_closed(true);
	shortest_path.set_right_closed(true);
	shortest_path.set_left_bounded();
	shortest_path.set_right_bounded();
	return shortest_path;
}

/** Comparator of points in lexicographical order */
struct bb_pair_compare {
	typedef std::pair<breakpoint,bool> bb_pair;
	bool operator()(const bb_pair& bp1, const bb_pair& bp2) {
		if (bp1.first.get_x()<bp2.first.get_x())
			return true;
		if (bp1.first.get_x()>bp2.first.get_x())
			return false;
		// x is equal, so compare y
		if (bp1.first.get_y()<bp2.first.get_y())
			return true;
		if (bp1.first.get_y()>bp2.first.get_y())
			return false;
		// x and y are equal, prefer bool=true (lower)
		if (bp1.second && !bp2.second)
			return true;
		else
			return false;
	}
};

/** Compute the slope between two points */
inline
precision_type slope(const std::pair<breakpoint,bool>& p1, const std::pair<breakpoint,bool>& p2) {
	// for points on the lower bound, use max of value and both limits,
	// for points on the upper bound, use min
	precision_type p1_x = p1.first.get_x();
	precision_type p2_x = p2.first.get_x();

	precision_type p1_y;
	if (p1.second)
		p1_y = max(p1.first.get_y_left(), p1.first.get_y_right(),
				p1.first.get_y());
	else
		p1_y = min(p1.first.get_y_left(), p1.first.get_y_right(),
				p1.first.get_y());
	precision_type p2_y;
	if (p2.second)
		p2_y = max(p2.first.get_y_left(), p2.first.get_y_right(),
				p2.first.get_y());
	else
		p2_y = min(p2.first.get_y_left(), p2.first.get_y_right(),
				p2.first.get_y());
	if (p1_x == p2_x) {
		return POS_INFTY * sgn(p2_y - p1_y);
	} else {
		return (p2_y - p1_y) / (p2_x - p1_x);
	}
}


/** Return whether the points p1,p2,p3 are (maybe) collinear.
 *
 * @author Goran Frehse */
inline bool is_collinear(const breakpoint& p1,
		const breakpoint& p2,const breakpoint& p3) {
	precision_type x1=p1.get_x();
	precision_type x2=p2.get_x();
	precision_type x3=p3.get_x();
	precision_type y1=p1.get_y_right();
	precision_type y2=p2.get_y();
	precision_type y3=p3.get_y_left();

	precision_type tolerance = max_abs(x1,x2,x3,y1,y2,y3) * 1e-6;
	return std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) < tolerance;
}

/** Return the intersection of lines p1,p2 and q1,q2
 *
 * @author Goran Frehse */
inline breakpoint compute_intersection_precise(const breakpoint& p1,
		const breakpoint& p2,const breakpoint& q1,
		const breakpoint& q2, bool& intersect) {

	precision_type ay_cy, ax_cx;
	precision_type dx_cx = q2.get_x() - q1.get_x();
	precision_type dy_cy = q2.get_y() - q1.get_y();
			precision_type bx_ax = p2.get_x() - p1.get_x();
			precision_type by_ay = p2.get_y() - p1.get_y();

	precision_type de = (bx_ax) * (dy_cy) - (by_ay) * (dx_cx);

	precision_type tolerance = max_abs(p2.get_x(), p1.get_x(), p2.get_y(),
			p1.get_y(), q2.get_x(), q1.get_x(), q2.get_y(), q1.get_y()) * 1e-6;

	if (std::abs(de)<tolerance)
		intersect = false;

	ax_cx = p1.get_x() - q1.get_x(); //L1StartPoint.X - L2StartPoint.X;
	ay_cy = p1.get_y() - q1.get_y(); //L1StartPoint.Y - L2StartPoint.Y;
	precision_type r = ((ay_cy) * (dx_cx) - (ax_cx) * (dy_cy)) / de;
	precision_type s = ((ay_cy) * (bx_ax) - (ax_cx) * (by_ay)) / de;
	precision_type px = p1.get_x() + r * (bx_ax);
	precision_type py = p1.get_y() + r * (by_ay);

	intersect = true;
	return breakpoint(px,py); //indicate there is intersection
}

/** Return the intersection of lines p1,p2 and q1,q2
 *
 * @author Goran Frehse */
inline breakpoint compute_intersection(const breakpoint& p1,
		const breakpoint& p2,const breakpoint& q1,
		const breakpoint& q2) {
	precision_type x1=p1.get_x();
	precision_type x2=p2.get_x();
	precision_type x3=q1.get_x();
	precision_type x4=q2.get_x();
	precision_type y1=p1.get_y_right();
	precision_type y2=p2.get_y_left();
	if (x2<x1) {
		y1=p1.get_y_left();
		y2=p2.get_y_right();
	}
	precision_type y3=q1.get_y_right();
	precision_type y4=q2.get_y_left();
	if (x4<x3) {
		y3=q1.get_y_left();
		y4=q2.get_y_right();
	}

	precision_type x12, x34, y12, y34, a, b, c;

	x12 = x1 - x2;
	x34 = x3 - x4;
	y12 = y1 - y2;
	y34 = y3 - y4;

	c = (x12 * y34) - (y12 * x34);
	a = x1 * y2 - y1 * x2;
	b = x3 * y4 - y3 * x4;

	precision_type x = (a * x34 - b * x12) / c;
	precision_type y = (a * y34 - b * y12) / c;

	const precision_type eps(1e-10);
	// snap to limits if close
	if (x - x1 < eps && x1 - x < eps) {
		x = x1;
		y = y1;
	} else if (x - x2 < eps && x2 - x < eps) {
		x = x2;
		y = y2;
	} else if (x - x3 < eps && x3 - x < eps) {
		x = x3;
		y = y3;
	} else if (x - x4 < eps && x4 - x < eps) {
		x = x4;
		y = y4;
	}

	return breakpoint(x,y);
}

/** Return the intersection of lines p1,p2 and q1,q2
 *
 * @author Goran Frehse */
inline breakpoint compute_intersection(const std::pair<breakpoint, bool>& p1,
		const std::pair<breakpoint, bool>& p2,const std::pair<breakpoint, bool>& q1,
		const std::pair<breakpoint, bool>& q2) {

	return compute_intersection(p1.first,p2.first,q1.first,q2.first);
}


/**
 * Find a greedy piecewise concave function from the left between lower and upper
 * @author Goran Frehse
 */
inline
piecewise_linear_function pw_convex_approx_greedy_left(const piecewise_linear_function& lower,
		const piecewise_linear_function& upper) {
	// for interval functions

	//piecewise_linear_function res = piecewise_linear_function::create_closed();

	// get a sorted list of breakpoints, each with the info whether it's a lower or upper bound point
	typedef std::pair<breakpoint, bool> bb_pair;

	size_t N = lower.get_list().size();
	size_t M = upper.get_list().size();

	assert(N>0);
	assert(M>0);
	// @todo figure out what to do if one of the lists is empty (do we care?)

	// get vectors with bool

	std::vector<bb_pair> lower_points;
	lower_points.reserve(2*lower.get_list().size());
	for (piecewise_linear_function::list_of_points_type::const_iterator it =
			lower.get_list().begin(); it != lower.get_list().end(); ++it) {
		if (it->get_y() != it->get_y_left())
			lower_points.push_back(
					std::make_pair(breakpoint(it->get_x(), it->get_y_left()),
							true));
		lower_points.push_back(
				std::make_pair(breakpoint(it->get_x(), it->get_y()), true));
		if (it->get_y() != it->get_y_right())
			lower_points.push_back(
					std::make_pair(breakpoint(it->get_x(), it->get_y_right()),
							true));
	}
	std::vector<bb_pair> upper_points;
	upper_points.reserve(2*upper.get_list().size());
	for (piecewise_linear_function::list_of_points_type::const_iterator it =
			upper.get_list().begin(); it != upper.get_list().end(); ++it) {
		if (it->get_y() != it->get_y_left())
			upper_points.push_back(
					std::make_pair(breakpoint(it->get_x(), it->get_y_left()),
							false));
		upper_points.push_back(
				std::make_pair(breakpoint(it->get_x(), it->get_y()), false));
		if (it->get_y() != it->get_y_right())
			upper_points.push_back(
					std::make_pair(breakpoint(it->get_x(), it->get_y_right()),
							false));
	}
	// merge two vectors supposing they're already sorted
	std::vector<bb_pair> points(lower_points.size()+upper_points.size());
	std::merge(lower_points.begin(), lower_points.end(), upper_points.begin(),
			upper_points.end(), points.begin(),bb_pair_compare());
//	std::vector<bb_pair> points(lower_points.begin(), lower_points.end());
//	assert(N == lower_points.size());
//	assert(M == upper_points.size());

//	std::vector<bb_pair> points(lower_points.size()+upper_points.size());
//	std::copy(lower_points.begin(),lower_points.end(),points.begin());
//	std::copy(upper_points.begin(),upper_points.end(),points.begin()+lower_points.size());
//	std::sort(points.begin(),points.end(),bb_pair_compare());


	std::vector<breakpoint> f_list;
	f_list.reserve(points.size());
	if (points.size() > 0) {
		breakpoint lC = points[0].first;
		f_list.push_back(lC);
		size_t pC(0), lv(0), uv(0), lP(0), lC_index(0);
		while (pC + 1 < points.size()) {
			if (f_list.size()>2*points.size()) {
				throw std::runtime_error("greed algo found too many points");
			}
			interval visible = interval::whole();
			while (!empty(visible) && pC + 1 < points.size()) {
				++pC; // go to the next point
//				std::cout << "lC:" << lC << " pC:" << points[pC].first << std::endl;
				precision_type sC = slope(std::make_pair(lC,true), points[pC]);
				if (in(sC, visible)) {
					if (points[pC].second) {
						// pC in L
						lv = pC;
						lP = pC;
						visible = interval(sC, visible.upper()); // restrict visibility from below
						//std::cout << "lower " << points[pC].first << " visible, restricting from below"<< std::endl;
					} else {
						// pC in U
						uv = pC;
						visible = interval(visible.lower(), sC); // restrict visibility from above
						//std::cout << "upper " << points[pC].first << " visible, restricting from above"<< std::endl;
					}
				} else if (sC < visible.lower()) {
					if (points[pC].second) {
						// pC in L
						lP = pC;
						//std::cout << "lower " << points[pC].first << " invisible below, becomes last lower"<< std::endl;
					} else {
						// pC in U
						lC_index = lv;
						lC = points[lv].first;
						if (lC.get_x() > f_list.rbegin()->get_x()
								|| lC.get_y() != f_list.rbegin()->get_y())
							f_list.push_back(lC);
						visible = interval(visible.lower(), sC); // restrict visibility from above
						assert(empty(visible));
//						std::cout << "upper " << points[pC].first << " invisible below, stop, insert " << lC << std::endl;
					}
				} else if (sC > visible.upper()) {
					if (points[pC].second) {
//						std::cout << "lower " << points[pC].first << " invisible above from " << lC;
						// pC in L
						lC_index = lP;
						//lC = compute_intersection(points[lv], points[uv], points[lP], points[pC]);
						// GF: start ray from last chosen lower point
						lC = compute_intersection(std::make_pair(lC,true), points[uv], points[lP], points[pC]);
						if (lC.get_x() > f_list.rbegin()->get_x()
								|| lC.get_y() != f_list.rbegin()->get_y())
							f_list.push_back(lC);
						visible = interval(sC, visible.upper()); // restrict visibility from below
//						std::cout << ", found intersection " << lC << ", restricting from below"<< std::endl;
					} else {
						//std::cout << "upper invisible above" << points[pC].first << " ignored"<< std::endl;
					}
				} else {
					//std::cout << "point " << points[pC].first << " ignored"<< std::endl;
				}
			}
			if (pC + 1 == points.size()) {
				// last point
				if (lC.get_x() > f_list.rbegin()->get_x()
						|| lC.get_y() != f_list.rbegin()->get_y())
					f_list.push_back(lC);

//				// Colas' fix for the transparent L hills problem
//				lC_index = lv;
//				lC = points[lv].first;
//				res.insert_point(lC);
//				pC = lv;
				//std::cout << "reached end, keeping " << lC << std::endl;
			} else {
				// pC = lC
				pC = lC_index;
				// GF: set pC to the next point to the left of lC;
				// ++pC will advance it to the point right of lC
				while (pC + 1 < points.size() && (!points[(pC+1)].second || points[(pC+1)].first.get_x()
						< lC.get_x())) {
					++pC;
				}
				//std::cout << "advancing to " << points[pC].first << std::endl;
			}
		}
		// GF: if no point on last x coordinate selected, add lower bound
		breakpoint last_lower = lower_points.back().first;
		//std::cout << last_lower << " vs " << *res.get_list().rbegin() << std::endl;
		if (!f_list.empty() && f_list.back().get_x()<last_lower.get_x())
			f_list.push_back(last_lower);
	}
	// copy the points to a vector
//	std::vector<breakpoint> res(f_list.begin(),f_list.end());
	// sort because it seems they're not all in order
	// std::sort(res.begin(),res.end());

//	return piecewise_linear_function::create_closed_from_temp(res);
	return piecewise_linear_function::create_closed_from_temp(f_list);
}

/** Check if three breakpoints are in counterclockwise order
 *
 * It is assumed that the breakpoints are continuous.
 *
 * Returns a negative value if the points are definitely ccw,
 * zero if they are maybe ccw and a positive value if they are
 * definitely not ccw.
 *
 * @author Goran Frehse */
inline
int is_counterclockwise_order(const breakpoint& A,
		const breakpoint& B, const breakpoint& C, precision_type rel_tolerance = 1e-8, precision_type abs_tolerance = 1e-12) {
//	precision_type m = A.get_x() * B.get_y() - A.get_y() * B.get_x() + A.get_y() * C.get_x() - A.get_x() * C.get_y() + B.get_x() * C.get_y()
//			- C.get_x() * B.get_y();
//	return m >= 1e-8; //1e-6;

//	bool result = ((x[0] - y[0]) * (z[1] - y[1]) - (x[1] - y[1])
//			* (z[0] - y[0]) < scalar_type(0));

	// Numerically more robust test:
	precision_type m = (A.get_x() - B.get_x()) * (C.get_y() - B.get_y()) - (A.get_y() - B.get_y()) * (C.get_x() - B.get_x());
	precision_type tolerance = max_abs(A.get_x(), A.get_y(),B.get_x(), B.get_y(),C.get_x(), C.get_y()) * rel_tolerance;
	tolerance = max(tolerance,abs_tolerance);
	tolerance = abs_tolerance;

	if (m < -tolerance) {
		return -1;
	} else if (m > tolerance) {
		return 1;
	} else {
		return 0;
	}
}
;


/** Find inflection breakpoints in bounded PLF, as well as their left neighbors
 *
 * If from_left=false, then the right neighbor is returned.
 * @attention If two consecutive points are identical, the result may be different than expected;
 *
 * @remark Used some ideas from
 * @cite Jerzy W. Jaromczyk, G.W. Wasilkowski, Computing convex hull in a floating point arithmetic, Computational Geometry, Volume 4, Issue 5, October 1994, Pages 283-292
 *
 * @author Goran Frehse */
inline
std::pair<std::vector<precision_type>, std::vector<precision_type> > inflection_points(
		const piecewise_linear_function& f, bool from_left = true, precision_type rel_tolerance = 1e-12, precision_type abs_tolerance = 1e-12) {
	std::vector<precision_type> pts;
	std::vector<precision_type> neighb_pts; // point to the left of the inflection point
	if (f.get_list().size() >= 3) {
		piecewise_linear_function::list_of_points_type::const_iterator i1,i2,i3;
		i1 = f.get_list().begin();
		i2 = i1;
		++i2;
		i3 = i2;
		++i3;
		bool first_unknown = true;
		while (i3 != f.get_list().end()) {
			int ccw = is_counterclockwise_order(*i1, *i2, *i3,rel_tolerance,abs_tolerance);
			if (ccw == 0) {
				// result unknown, move only i3
				if (first_unknown && i2->get_x()-i1->get_x()<i3->get_x()-i2->get_x()) {
					i2 = i3;
				}
				++i3;
			} else {
				first_unknown = false;
				// result known
				if (ccw < 0) {
					pts.push_back(i2->get_x());
					if (from_left) {
						neighb_pts.push_back((i2 - 1)->get_x());
					} else {
						neighb_pts.push_back((i2 + 1)->get_x());
					}
				}
				// move on
				i1 = i2;
				i2 = i3;
				++i3;
			}
		}
	}
	return make_pair(pts, neighb_pts);
}

/** An error class for numerical interval matching */
class interval_matching_error: public std::runtime_error
{
    public:
	interval_matching_error(std::string const& msg):
            std::runtime_error(msg)
        {}
};

/** Construct intervals from two vectors of interval bounds, on a given domain
 * @author Goran Frehse
 */
inline
std::vector<interval> compute_intervals(
		const std::vector<precision_type>& as,
		const std::vector<precision_type>& bs, const interval& domain) {
	std::vector<interval> itvs;


	// Build pairs of as and bs
	std::vector<precision_type>::const_iterator ait=as.begin();
	std::vector<precision_type>::const_iterator bit=bs.begin();
	// If there's a b without matching a in the beginning, take the left domain
	if (bit!=bs.end() && (ait==as.end() || *ait>*bit)) {
		// don't include an interval that just consists of the beginning of the domain
		if (*bit > domain.lower()) {
			itvs.push_back(interval(domain.lower(),*bit));
		}
		++bit;
	}
	while (bit!=bs.end()) {
		bool had_problems = false;
		// @todo This could be made more numerically stable
		// @attention If there's an interval end at the end of the domain, we'll suppose the
		// corresponding interval start is also the end (and therefore we don't need to include it)
		if (ait==as.end() && *bit < domain.upper()) {
//				std::cerr << std::endl << "domain:" << domain << std::endl;
//				std::cerr << "inflection starts:" << as << std::endl;
//				std::cerr << "inflection ends:" << bs << std::endl;
//				std::cerr << "intervals so far:" << itvs << std::endl;
//				std::cerr << "a:" << *ait << " b:" << *bit << " diff: " << *bit-*ait << std::endl;
			throw interval_matching_error("no matchin interval start");
		}
		if (*ait>*bit) {
			had_problems = true;
			++bit;
			// todo: should throw a warning here, but we don't have a warning mechanism
		}
		if ((ait+1)!=as.end() && *(ait+1)<=*bit) {
			had_problems = true;
			++ait;
			// todo: should throw a warning here, but we don't have a warning mechanism
		}

		if (!had_problems) {
			// don't include an interval that just consists of the beginning of the domain
			if (*ait < domain.upper() && *bit > domain.lower()) {
				itvs.push_back(interval(*ait,*bit));
			}
			++ait;
			++bit;
		}
	}
	// If there's an a without matching b in the end, take the right domain
	if (ait!=as.end()) {
		if (*ait < domain.upper()) {
			itvs.push_back(interval(*ait,domain.upper()));
		}
		++ait;
	}
	if (ait!=as.end()) {
//			std::cerr << std::endl << "domain:" << domain << std::endl;
//			std::cerr << "inflection starts:" << as << std::endl;
//			std::cerr << "inflection ends:" << bs << std::endl;
//			std::cerr << "intervals so far:" << itvs << std::endl;
		throw interval_matching_error("no matching interval end");
	}

	return itvs;
}

/** Find inflection point intervals for bounded PLF
 *
 * @author Goran Frehse */
inline
std::pair<std::vector<interval>, std::vector<interval> > get_inflection_intervals(
		const piecewise_linear_function& lower,
		const piecewise_linear_function& upper, precision_type rel_tolerance = 1e-8, precision_type abs_tolerance = 1e-12) {
	std::vector<interval> infl_itvs,neighb_itvs;

	piecewise_linear_function greedyleft = pw_convex_approx_greedy_left(lower,
			upper);
	// now in the other direction
	piecewise_linear_function greedyright = x_mirror(
			pw_convex_approx_greedy_left(x_mirror(lower), x_mirror(upper)));

//	std::cout << std::endl << "lower:"<<lower<<" upper:"<<upper << " greedyleft:" << greedyleft << std::endl;

	static size_t count(0);
	++count;
	std::stringstream ss;
	ss << count;

//	plf_graph(greedyleft, "gif", "/tmp/test_greedy_" + ss.str());
//	plf_graph(greedyright, "gif", "/tmp/test_greedy2_" + ss.str());
//	piecewise_linear_interval_function margin(lower, upper);
//	plif_graph(
//			margin,
//			"gif",
//			"/tmp/test_spacetime_margin_greedy_" + ss.str(),
//			"-m3 -W 0.00 /tmp/test_greedy_" + ss.str()
//					+ ".txt -s -m4 /tmp/test_greedy2_" + ss.str()
//					+ ".txt -s -W 0 ");

	bool had_problems = false;
	do {
		std::pair<std::vector<precision_type>, std::vector<precision_type> > a_infl_info =
				inflection_points(greedyright, false, rel_tolerance, abs_tolerance);
		std::vector<precision_type>& as = a_infl_info.first;
		std::pair<std::vector<precision_type>, std::vector<precision_type> > b_infl_info =
				inflection_points(greedyleft, true, rel_tolerance, abs_tolerance);
		std::vector<precision_type>& bs = b_infl_info.first;

		try {
			had_problems = false;
			infl_itvs = compute_intervals(as, bs, lower.get_domain());
			// GF: the following is very buggy, so I'm taking it out
//	neighb_itvs = compute_intervals(b_infl_info.second, a_infl_info.second,
//			lower.get_domain());
			neighb_itvs = std::vector<interval>(infl_itvs.size(),
					interval::whole());
		} catch (interval_matching_error& e) {
			had_problems = true;
			rel_tolerance = rel_tolerance * 1e2;

			//std::cerr << "increasing tolerance to " << rel_tolerance << std::endl;
			/*
			std::cerr << std::endl << "lower:" << lower << " upper:" << upper
					<< " greedyleft:" << greedyleft << " greedyright:"
					<< greedyright << std::endl;

			plf_graph(greedyleft, "gif", "/tmp/test_greedy_" + ss.str());
			plf_graph(greedyright, "gif", "/tmp/test_greedy2_" + ss.str());
			piecewise_linear_interval_function margin(lower, upper);
			plif_graph(margin, "X",
					"/tmp/test_spacetime_margin_greedy_" + ss.str(),
					"-m3 -W 0.00 /tmp/test_greedy_" + ss.str()
							+ ".txt -s -m4 /tmp/test_greedy2_" + ss.str()
							+ ".txt -s -W 0 ");
			*/
			/** @attention This is really a numerical error, so we should throw. */
			//throw e;
		}
	} while (had_problems && rel_tolerance < 1.0);

//	std::cout << "inflection intervals:" << infl_itvs << std::endl;
//	std::cout << "neighbor intervals:" << neighb_itvs << std::endl;

	assert(infl_itvs.size()==neighb_itvs.size());
	for (size_t i=0;i<infl_itvs.size();++i) {
//		assert(boost::numeric::subset(infl_itvs[i],neighb_itvs[i]));
		if (!subset(infl_itvs[i],neighb_itvs[i])) {
			/** @todo This is an ad hoc fix, we shouldn't even get here */
			neighb_itvs[i] = hull(neighb_itvs[i],infl_itvs[i]);
		}
	}

	return std::make_pair(infl_itvs,neighb_itvs);
}

/**
 * Find a greedy minimum link piecewise concave function from the left between lower and upper
 * @attention doesn't work, output can go below lower.
 * @author Goran Frehse
 */
inline
piecewise_linear_function min_link_greedy_left(const piecewise_linear_function& lower,
		const piecewise_linear_function& upper) {
	// for interval functions

	//piecewise_linear_function res = piecewise_linear_function::create_closed();

	// get a sorted list of breakpoints, each with the info whether it's a lower or upper bound point
	typedef std::pair<breakpoint, bool> bb_pair;

	size_t N = lower.get_list().size();
	size_t M = upper.get_list().size();

	assert(N>0);
	assert(M>0);
	// @todo figure out what to do if one of the lists is empty (do we care?)

	// get vectors with bool

	std::vector<bb_pair> lower_points;
//	lower_points.reserve(2*lower.get_list().size());
	for (piecewise_linear_function::list_of_points_type::const_iterator it =
			lower.get_list().begin(); it != lower.get_list().end(); ++it) {
		if (it->get_y() != it->get_y_left())
			lower_points.push_back(
					std::make_pair(breakpoint(it->get_x(), it->get_y_left()),
							true));
		lower_points.push_back(
				std::make_pair(breakpoint(it->get_x(), it->get_y()), true));
		if (it->get_y() != it->get_y_right())
			lower_points.push_back(
					std::make_pair(breakpoint(it->get_x(), it->get_y_right()),
							true));
	}
	std::vector<bb_pair> upper_points;
//	upper_points.reserve(2*upper.get_list().size());
	for (piecewise_linear_function::list_of_points_type::const_iterator it =
			upper.get_list().begin(); it != upper.get_list().end(); ++it) {
		if (it->get_y() != it->get_y_left())
			upper_points.push_back(
					std::make_pair(breakpoint(it->get_x(), it->get_y_left()),
							false));
		upper_points.push_back(
				std::make_pair(breakpoint(it->get_x(), it->get_y()), false));
		if (it->get_y() != it->get_y_right())
			upper_points.push_back(
					std::make_pair(breakpoint(it->get_x(), it->get_y_right()),
							false));
	}
	// merge two vectors supposing they're already sorted
//	std::vector<bb_pair> points(N + M);
//	std::merge(lower_points.begin(), lower_points.end(), upper_points.begin(),
//			upper_points.end(), points.begin());
//	std::vector<bb_pair> points(lower_points.begin(), lower_points.end());
//	assert(N == lower_points.size());
//	assert(M == upper_points.size());

	std::vector<bb_pair> points(lower_points.size()+upper_points.size());
	std::copy(lower_points.begin(),lower_points.end(),points.begin());
	std::copy(upper_points.begin(),upper_points.end(),points.begin()+lower_points.size());
	std::sort(points.begin(),points.end(),bb_pair_compare());


	std::deque<breakpoint> f_list;
	if (points.size() > 0) {
		size_t last_visible_lower(0);
		size_t last_visible_upper(0);
		breakpoint last_chosen(points[0].first);
		f_list.push_back(last_chosen);
		size_t pC(1);
		while (pC+1 < points.size()) {
			if (f_list.size()>2*points.size()) {
				throw std::runtime_error("greedy min link algo found too many points");
			}
			interval visible = interval::whole();
			while (!empty(visible) && pC < points.size()) {
//				std::cout << "lC:" << lC << " pC:" << points[pC].first << std::endl;
				precision_type sC = slope(std::make_pair(last_chosen,true), points[pC]);
				if (in(sC, visible)) {
					if (points[pC].second) {
						// pC in L
						last_visible_lower = pC;
						visible = interval(sC, visible.upper()); // restrict visibility from below
//						std::cout << "lower " << points[pC].first << " visible, restricting from below"<< std::endl;
					} else {
						// pC in U
						last_visible_upper = pC;
						visible = interval(visible.lower(), sC); // restrict visibility from above
//						std::cout << "upper " << points[pC].first << " visible, restricting from above"<< std::endl;
					}
				} else if (sC < visible.lower()) {
					if (points[pC].second) {
						// pC in L
						//std::cout << "lower " << points[pC].first << " invisible below, becomes last lower"<< std::endl;
					} else {
//						std::cout << "upper " << points[pC].first << " invisible above from " << last_chosen;
						// pC in U
						// GF: start ray from last chosen point
						breakpoint lC = compute_intersection(std::make_pair(last_chosen, true),
								points[last_visible_lower],
								points[last_visible_upper], points[pC]);
						if (lC.get_x() > f_list.rbegin()->get_x()
								|| lC.get_y() != f_list.rbegin()->get_y()) {
							last_chosen = lC;
							f_list.push_back(lC);
						}
						visible = interval::empty();
//						visible = interval::whole();
//						visible = interval(visible.lower(),
//								slope(points[last_visible_upper], std::make_pair(last_chosen,false))); // restrict visibility from above
//						std::cout << ", found intersection " << lC << ", restarting"<< std::endl;
					}
				} else if (sC > visible.upper()) {
					if (points[pC].second) {
//						std::cout << "lower " << points[pC].first << " invisible above from " << last_chosen;
						// pC in L
						// GF: start ray from last chosen point
						breakpoint lC = compute_intersection(std::make_pair(last_chosen, true),
								points[last_visible_upper],
								points[last_visible_lower], points[pC]);
						if (lC.get_x() > f_list.rbegin()->get_x()
								|| lC.get_y() != f_list.rbegin()->get_y()) {
							last_chosen = lC;
							f_list.push_back(lC);
						}
						visible = interval::empty();
//						visible = interval::whole();
//						visible = interval(slope(points[last_visible_lower], std::make_pair(last_chosen,false)),visible.upper()); // restrict visibility from above
//						std::cout << ", found intersection " << lC << ", restarting"<< std::endl;
					} else {
						//std::cout << "upper invisible above" << points[pC].first << " ignored"<< std::endl;
					}
				} else {
					//std::cout << "point " << points[pC].first << " ignored"<< std::endl;
				}
				++pC; // go to the next point
			}
			// back up pC to last point to the right of last chosen
			while (pC < points.size() && pC>0 && points[pC-1].first.get_x() >= last_chosen.get_x()) {
				--pC;
			}
		}
		// GF: if no point on last x coordinate selected, add lower bound
		breakpoint last_lower = lower_points.back().first;
		//std::cout << last_lower << " vs " << *res.get_list().rbegin() << std::endl;
		if (!f_list.empty() && f_list.back().get_x()<last_lower.get_x())
			f_list.push_back(last_lower);
	}
	// copy the points to a vector
	std::vector<breakpoint> res(f_list.begin(),f_list.end());
	// sort because it seems they're not all in order
	// std::sort(res.begin(),res.end());

	piecewise_linear_function res_f = piecewise_linear_function::create_closed_from_temp(res);

	return res_f;
}

/**
 * Find a greedy minimum link piecewise concave function by relaxing between lower and upper
 * @attention It is assumed that both lower and upper are continuous.
 * @author Goran Frehse
 */
inline
piecewise_linear_function min_link_concave_relax(const piecewise_linear_function& lower,
		const piecewise_linear_function& upper) {

	//piecewise_linear_function res = piecewise_linear_function::create_closed();

	size_t N = lower.get_list().size();
	size_t M = upper.get_list().size();

	if (N<=4)
		return lower;
	assert(N>0);
	assert(M>0);
	// @todo figure out what to do if one of the lists is empty (do we care?)

	// get vectors with bool

	std::deque<breakpoint> f_list;

	typedef piecewise_linear_function::list_of_points_type::const_iterator const_iterator;
	const piecewise_linear_function::list_of_points_type& lpts = lower.get_list();

	const_iterator dit = lpts.begin();
	breakpoint a(*dit);
	++dit;
	breakpoint b(*dit);
	++dit;
	breakpoint c(*dit);
	++dit;
	breakpoint d(*dit);

	do {
		d = *dit;
		// first, test for collinearity
		if (false && is_collinear(b,c,d)) {
			// drop b
//			std::cout << "collinear" << std::endl;
			b = c;
			++dit;
		} else {
			bool intersects;
			breakpoint f = compute_intersection_precise(a, b, c, d, intersects);
			bool erase_c = false;
			if (!intersects) {
				// erase b and c
//				std::cout << "parallel" << std::endl;
				b = c;
				c = d;
				++dit;
				if (dit != lpts.end()) {
					b = c; // c = d
					c = *dit;
					++dit;
				}
			} else {
				if (f.get_x() >= b.get_x() && f.get_x() <= c.get_x()) {
					// get the value of the upper bound at this point
					precision_type ey = upper.value_at(f.get_x());
//					std::cout << "f.x:" << f.get_x() << " f.y:" << f.get_y()
//							<< " e.y:" << ey << std::endl;
					if (f.get_y() <= ey) {
						erase_c = true;
					}
				}
				if (erase_c) {
					b = f;
					c = d;
					++dit;
				} else {
					// advance to next
					f_list.push_back(a);
					a = b;
					b = c;
					c = d;
					++dit;
				}
			}
		}
	} while (dit != lpts.end());
	// push back a,b,c (d is beyond)
	f_list.push_back(a);
	f_list.push_back(b);
	f_list.push_back(c);

	// copy the points to a vector
	std::vector<breakpoint> res(f_list.begin(),f_list.end());
	// sort because it seems they're not all in order
	// std::sort(res.begin(),res.end());

	piecewise_linear_function res_f = piecewise_linear_function::create_closed_from_temp(res);

	return res_f;
}

}

#endif /* PIECEWISE_LINEAR_FUNCTION_PATHS_H_ */
