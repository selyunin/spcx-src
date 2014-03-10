/*
 * piecewise_linear_function_hulls.h
 *
 *  Created on: Oct 18, 2012
 *      Author: kateja
 */

#ifndef PIECEWISE_LINEAR_FUNCTION_HULLS_H_
#define PIECEWISE_LINEAR_FUNCTION_HULLS_H_

#include "piecewise_linear_function_operators.h"
#include "piecewise_linear_interval_function_utility.h" // for debug plots

namespace plif {

/**
 * Convex Hull of a PLF
 *
 * The convex hull of a PLF is its largest convex lower bound.
 */
inline
piecewise_linear_function convex_hull(const piecewise_linear_function& f)
{
	assert(!f.is_left_unbounded() && !f.is_right_unbounded());
	std::vector<breakpoint> convex_hull_list;
	convex_hull_list.reserve(f.get_list().size());

	const std::vector<breakpoint>& list = f.get_list();
	int i;
	precision_type slope_prev, slope_next, x1, y1, x2, y2;
	/**
	 *When the PLF has only two points, the convex hull is the same as the PLF
	 */
	if(list.size() == 0)
	{
		// don't do anything
	}
	else if(list.size() == 1)
	{
		x1 = list[0].get_x();
		y1 = list[0].get_y();
		convex_hull_list.push_back(breakpoint(x1, y1, y1, y1));
	}
	else if(list.size() == 2)
	{
		x1 = list[0].get_x();
		if(f.is_left_closed())
			y1 = list[0].inf_right();
		else
			y1 = list[0].get_y_right();
		x2 = list[1].get_x();
		if(f.is_right_closed())
			y2 = list[1].inf_left();
		else
			y2 = list[1].get_y_left();
		convex_hull_list.push_back(breakpoint(x1, y1));
		convex_hull_list.push_back(breakpoint(x2, y2));
	}
	/**
	 *Else Case: The PLF has more than two points, sequentially check each of the points
	 *according to the slope it can be decided whether the function need to be included or not
	 */
	else
	{
		precision_type start_x;
		precision_type start_y;
		start_x = list[0].get_x();
		if(f.is_left_closed())
			start_y = list[0].inf_right();
		else
			start_y = list[0].get_y_right();
		x1 = start_x;
		y1 = start_y;
		x2 = list[1].get_x();
		y2 = list[1].inf();
		slope_prev = (y2-y1)/(x2-x1);
		convex_hull_list.push_back(breakpoint(x1, y1));
		convex_hull_list.push_back(breakpoint(x2, y2));
		x1 = x2;
		y1 = y2;
		for(i=2;i<list.size();i++)
		{
			x2 = list[i].get_x();
			if(i==list.size()-1 && !f.is_right_closed())
				y2 = list[i].get_y_left();
			else
				y2 = list[i].inf();
			slope_next = (y2-y1)/(x2-x1);
			if(slope_next>=slope_prev)
			{
				convex_hull_list.push_back(breakpoint(x2, y2));
				slope_prev = slope_next;
				x1 = x2;
				y1 = y2;
			}
			else
			{
				convex_hull_list.pop_back();
				if(convex_hull_list.size()>=2)
				{
					x1 = convex_hull_list[convex_hull_list.size()-1].get_x();
					y1 = convex_hull_list[convex_hull_list.size()-1].get_y();
					slope_prev = (y1-convex_hull_list[convex_hull_list.size()-2].get_y())/(x1-convex_hull_list[convex_hull_list.size()-2].get_x());
					if (i>0) {
						i--;
					}
				}
				else
				{
					x1 = x2;
					y1 = y2;
					slope_prev = (y1 - start_y)/(x1 - start_x);
					convex_hull_list.push_back(breakpoint(x1, y1));
				}

			}
		}
	}
	piecewise_linear_function convex_hull=piecewise_linear_function::create_closed_from_temp(convex_hull_list);

//	f.display();
//	convex_hull.display();
//	static size_t counter = 0;
//	++counter;
//plif::piecewise_linear_interval_function margin(convex_hull,f);
//plif::plif_graph(margin,"png","/tmp/test_chull_"+to_string(counter));

	return convex_hull;
}

/**
 * Concave Hull of a PLF
 *
 * The concave hull of a PLF is its smallest concave upper bound.
 */
inline
piecewise_linear_function concave_hull(const piecewise_linear_function& f)
{
	return -convex_hull(-f);
}

}

#endif /* PIECEWISE_LINEAR_FUNCTION_HULLS_H_ */
