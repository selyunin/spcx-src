/*
 * piecewise_linear_interval_function_operators.h
 *
 *  Created on: Oct 18, 2012
 *      Author: kateja
 */

#ifndef PIECEWISE_LINEAR_INTERVAL_FUNCTION_OPERATORS_H_
#define PIECEWISE_LINEAR_INTERVAL_FUNCTION_OPERATORS_H_

#include "geometry.h"
#include "piecewise_linear_interval_function.h"


namespace plif {

/**
 *
 * Pointwise addition
 *
 * both have bounded domains
 * addition is performed on the intersection of the two functions domain
 *
 */

inline
piecewise_linear_interval_function operator+(const piecewise_linear_interval_function& f, const piecewise_linear_interval_function& g)
{
	piecewise_linear_interval_function h(f.get_lower()+g.get_lower(),f.get_upper()+g.get_upper());
	return h;
}

/**
 *
 * Negation
 *
 */

inline
piecewise_linear_interval_function operator-(const piecewise_linear_interval_function& f)
{
	piecewise_linear_interval_function h(-f.get_upper(),-f.get_lower());
	return h;
}


/**
 * Scalar Multiplication for the PLIF
 */
inline
piecewise_linear_interval_function operator*(const piecewise_linear_interval_function& f, precision_type lambda)
{
	piecewise_linear_interval_function lambda_f(f.get_lower()*lambda,f.get_upper()*lambda);
	return lambda_f;
}

/**
 * Scalar Multiplication for the PLIF
 */
inline
piecewise_linear_interval_function operator*(precision_type lambda, const piecewise_linear_interval_function& f)
{
	return f*lambda;
}


/**
 * Intersection of two PLIF's f and g
 */
inline
piecewise_linear_interval_function intersection(const piecewise_linear_interval_function& f, const piecewise_linear_interval_function& g)
{
	if (!empty(f.get_domain()) && !empty(f.get_domain())) {
		piecewise_linear_function lower = pointwise_maximum(f.get_lower(),
				g.get_lower());
		piecewise_linear_function upper = -pointwise_maximum(-f.get_upper(),
				-g.get_upper());
		piecewise_linear_interval_function intersection_f_g(lower, upper);
		return intersection_f_g;
	} else {
		piecewise_linear_interval_function intersection_f_g;
		intersection_f_g.set_domain(interval::empty());
		return intersection_f_g;
	}
}


/**
 * Splits a bounded PLIF into connected components, returns the vector<PLIF> as the vector of all the connected components
 */
/*
inline
std::vector<piecewise_linear_interval_function> split(const piecewise_linear_interval_function& f)
{
	assert(	   !f.get_lower().is_left_unbounded() 
			&& !f.get_lower().is_right_unbounded() 
			&& !f.get_upper().is_left_unbounded() 
			&& !f.get_upper().is_right_unbounded() );
	std::vector<piecewise_linear_interval_function> split_f;
	std::vector<breakpoint> list_lower, list_upper;
	std::vector<interval> domains;
	piecewise_linear_function lower, upper;
	std::vector<precision_type> list;
	int i, j;
	lower = f.get_lower();
	upper = f.get_upper();
	interval domain = intersect(lower.get_domain(), upper.get_domain());

	list_lower.clear();
	list_upper.clear();
	split_f.clear();
	domains.clear();
	lower.get_list(list_lower);
	upper.get_list(list_upper);
	for(i=0, j=0;i<list_lower.size() || j<list_upper.size();)
	{
		if(j>=list_upper.size() || (i<list_lower.size() && list_lower[i].get_x() < list_upper[j].get_x()))
		{
			if(belongs(list_lower[i].get_x(), domain))
				list.push_back(list_lower[i].get_x());
			i++;

		}
		else if(i>=list_lower.size() || (j<list_upper.size() && list_lower[i].get_x() > list_upper[j].get_x()))
		{
			if(belongs(list_upper[j].get_x(), domain))
				list.push_back(list_upper[j].get_x());
			j++;

		}
		else
		{
			if(belongs(list_upper[i].get_x(), domain))
				list.push_back(list_upper[j].get_x());
			i++;
			j++;

		}
	}

	for(i=0;i<list.size()-1;i++)
	{

		piecewise_linear_function low, up;
		precision_type prev = list[i];
		precision_type next = list[i+1];
		bool R, L;
		bool left_closed = false, right_closed = false;
		R = lower.right_lim_at(prev) <= upper.right_lim_at(prev);
		L = lower.left_lim_at(next) <= upper.left_lim_at(next);
		low.set_left_closed(false);
		low.set_right_closed(false);
		up.set_left_closed(false);
		up.set_right_closed(false);
		if(i!=0 && lower.value_at(prev)<=upper.value_at(prev))
		{
			left_closed = true;
		}
		else if( i==0 && ((list_lower[i].get_x() > list_upper[i].get_x() && lower.is_left_closed())
						 ||(list_lower[i].get_x() < list_upper[i].get_x() && upper.is_left_closed() )
						 ||(list_lower[i].get_x()==list_upper[i].get_x() && lower.is_left_closed() && upper.is_left_closed()))
					  && (lower.value_at(prev)<=upper.value_at(prev)) )
		{
			left_closed = true;
		}
		int lindex = list_lower.size()-1;
		int uindex = list_upper.size()-1;
		if(i!=list.size()-2 && lower.value_at(next)<=upper.value_at(next))
		{
			right_closed = true;
		}

		else if(i==list.size()-2 && ((list_lower[lindex].get_x() < list_upper[uindex].get_x() && lower.is_right_closed())
									||(list_lower[lindex].get_x() > list_upper[uindex].get_x()&&upper.is_right_closed())
									||(list_lower[lindex].get_x()==list_upper[uindex].get_x()&&lower.is_right_closed()&&upper.is_right_closed()))
								 &&	(lower.value_at(next)<=upper.value_at(next)))
		{
			right_closed = true;
		}
		if(R && L)
		{

			piecewise_linear_interval_function func;
			precision_type start = prev;
			precision_type end = next;

			breakpoint start_up, end_up, start_low, end_low;
			start_up.set_x(start);
			start_low.set_x(start);
			end_up.set_x(end);
			end_low.set_x(end);
			interval domain(start, end);
			if(left_closed)
			{
				start_up.set_y(upper.value_at(start));
				start_low.set_y(lower.value_at(start));
			}
			if(right_closed)
			{
				end_up.set_y(upper.value_at(end));
				end_low.set_y(lower.value_at(end));
			}
			start_up.set_y_right(upper.right_lim_at(start));
			start_low.set_y_right(lower.right_lim_at(start));
			end_up.set_y_left(upper.left_lim_at(end));
			end_low.set_y_left(lower.left_lim_at(end));


			low.copy_on_domain(lower, start, end);
			up.copy_on_domain(upper, start, end);

			low.insert_point(start_low);
			low.insert_point(end_low);
			up.insert_point(start_up);
			up.insert_point(end_up);

			low.set_left_closed(left_closed);
			up.set_left_closed(left_closed);
			low.set_right_closed(right_closed);
			up.set_right_closed(right_closed);
			func.set_lower(low);
			func.set_upper(up);
			func.set_domain(domain);
			func.set_left_closed(left_closed);
			func.set_right_closed(right_closed);
			split_f.push_back(func);
		}
		else if(!R && L)
		{
			piecewise_linear_interval_function func;
			precision_type start, end;

			start = intersection(prev, next, lower.right_lim_at(prev), lower.left_lim_at(next), upper.right_lim_at(prev), upper.left_lim_at(next), false, false);
			end = next;
			interval domain(start, end);

			if(end >= start)
			{
				left_closed = true;
				low.set_left_closed(left_closed);
				up.set_left_closed(left_closed);
				low.set_right_closed(right_closed);
				up.set_right_closed(right_closed);




				breakpoint start_up, end_up, start_low, end_low;
				start_up.set_x(start);
				start_low.set_x(start);
				end_up.set_x(end);
				end_low.set_x(end);
				if(left_closed)
				{
					start_up.set_y(upper.value_at(start));
					start_low.set_y(lower.value_at(start));
				}
				if(right_closed)
				{
					end_up.set_y(upper.value_at(end));
					end_low.set_y(lower.value_at(end));
				}
				start_up.set_y_right(upper.right_lim_at(start));
				start_low.set_y_right(lower.right_lim_at(start));
				end_up.set_y_left(upper.left_lim_at(end));
				end_low.set_y_left(lower.left_lim_at(end));


				low.copy_on_domain(lower, start, end);
				up.copy_on_domain(upper, start, end);

				low.insert_point(start_low);
				low.insert_point(end_low);
				up.insert_point(start_up);
				up.insert_point(end_up);



				func.set_lower(low);
				func.set_upper(up);
				func.set_domain(domain);
				func.set_left_closed(left_closed);
				func.set_right_closed(right_closed);
				split_f.push_back(func);
		}
		}
		else if (R && !L)
		{

			piecewise_linear_interval_function func;
			precision_type start, end;
			start = prev;
			end = intersection(prev, next, lower.right_lim_at(prev), lower.left_lim_at(next), upper.right_lim_at(prev), upper.left_lim_at(next), false, false);
			interval domain(start, end);
			if(end>=start)
			{
				right_closed = true;
				low.set_left_closed(left_closed);
				up.set_left_closed(left_closed);
				low.set_right_closed(right_closed);
				up.set_right_closed(right_closed);

				breakpoint start_up, end_up, start_low, end_low;
				start_up.set_x(start);
				start_low.set_x(start);
				end_up.set_x(end);
				end_low.set_x(end);
				if(left_closed)
				{
					start_up.set_y(upper.value_at(start));
					start_low.set_y(lower.value_at(start));
				}
				if(right_closed)
				{
					end_up.set_y(upper.value_at(end));
					end_low.set_y(lower.value_at(end));
				}
				start_up.set_y_right(upper.right_lim_at(start));
				start_low.set_y_right(lower.right_lim_at(start));
				end_up.set_y_left(upper.left_lim_at(end));
				end_low.set_y_left(lower.left_lim_at(end));


				low.copy_on_domain(lower, start, end);
				up.copy_on_domain(upper, start, end);

				low.insert_point(start_low);
				low.insert_point(end_low);
				up.insert_point(start_up);
				up.insert_point(end_up);


				func.set_lower(low);
				func.set_upper(up);
				func.set_domain(domain);
				func.set_left_closed(left_closed);
				func.set_right_closed(right_closed);
				split_f.push_back(func);
			}
		}
		else if(!R && !L)
		{}

	}
	std::vector<piecewise_linear_interval_function>::iterator pointer = split_f.begin();
	for(int i=0;i+1<split_f.size();i++)
	{
		pointer++;
		if(split_f[i].get_domain().upper() == split_f[i+1].get_domain().lower() && split_f[i].is_right_closed() && split_f[i+1].is_left_closed())
		{
			int j;
			interval new_domain(split_f[i].get_domain().lower(), split_f[i+1].get_domain().upper());
			bool left_closed = split_f[i].is_left_closed();
			bool right_closed = split_f[i+1].is_right_closed();
			std::vector<breakpoint> copy_points;
			piecewise_linear_function low_prev = split_f[i].get_lower();
			piecewise_linear_function up_prev = split_f[i].get_upper();
			piecewise_linear_function low_next = split_f[i+1].get_lower();
			piecewise_linear_function up_next = split_f[i+1].get_upper();
			low_next.get_list(copy_points);
			for(j=0;j<copy_points.size();j++)
			{
				
				if(j==0)
					low_prev.set_y_right_at(copy_points[j].get_x(), copy_points[j].get_y_right());
				low_prev.insert_point(copy_points[j]);
			}
			up_next.get_list(copy_points);
			for(j=0;j<copy_points.size();j++)
			{
				if(j==0)
					up_prev.set_y_right_at(copy_points[j].get_x(), copy_points[j].get_y_right());
				up_prev.insert_point(copy_points[j]);
			}

			low_prev.set_right_closed(low_next.is_right_closed());
			low_prev.set_domain(new_domain);
			up_prev.set_right_closed(up_next.is_right_closed());
			up_prev.set_domain(new_domain);
			split_f[i].set_lower(low_prev);
			split_f[i].set_upper(up_prev);

			split_f[i].set_domain(new_domain);
			split_f[i].set_left_closed(left_closed);
			split_f[i].set_right_closed(right_closed);


			pointer = split_f.erase(pointer);
			pointer--;
			i--;
		}
	}
	return split_f;
}
*/

/** Divides h into subdomains only overlapping at a given set of cut points.
  *
  * It is assumed that all cut points lie within the domain of h.
  * The domain of the first piece is the start of the domain of h to the 
  * first cut point , the second from the first cut point to
  * the second cut point, etc. The domain of the last piece is the last 
  * cut point to the end of the domain of h.
  */
inline
std::vector<piecewise_linear_interval_function> dissect(const piecewise_linear_interval_function& h, const std::vector<precision_type> cut_points)
{
	assert(	   !h.get_lower().is_left_unbounded() 
			&& !h.get_lower().is_right_unbounded() 
			&& !h.get_upper().is_left_unbounded() 
			&& !h.get_upper().is_right_unbounded() );
	std::vector<piecewise_linear_interval_function> dissection_list;
	// we know its size, so let's reserve the space
	dissection_list.reserve(cut_points.size()+1);

	// catch case of no cutpoints
	if (cut_points.size()==0) {
		dissection_list.push_back(h);
	} else {
		const piecewise_linear_function& low_h = h.get_lower();
		const piecewise_linear_function& up_h = h.get_upper();

		std::vector<piecewise_linear_function> lower_cuts = dissect(low_h,cut_points);
		std::vector<piecewise_linear_function> upper_cuts = dissect(up_h,cut_points);

		assert(lower_cuts.size()==upper_cuts.size());

		/** All the sub-domains in a loop */
		std::vector<piecewise_linear_function>::const_iterator low_it=lower_cuts.begin();
		std::vector<piecewise_linear_function>::const_iterator upp_it=upper_cuts.begin();
		while (low_it!=lower_cuts.end()) {
			piecewise_linear_interval_function h_piece(*low_it,*upp_it);
			dissection_list.push_back(h_piece);
			++low_it;
			++upp_it;
		}
	}

	return dissection_list;
}

/** Divides a PLIF f into subdomains given by a set of intervals
  *
  * @author Goran Frehse
  */
inline
std::vector<piecewise_linear_interval_function> dissect(
		const piecewise_linear_interval_function& h,
		const std::vector<interval> subdomains) {
	assert( !h.get_lower().is_left_unbounded()
			&& !h.get_lower().is_right_unbounded()
			&& !h.get_upper().is_left_unbounded()
			&& !h.get_upper().is_right_unbounded() );
	std::vector<piecewise_linear_interval_function> dissection_list;
	// we know its size, so let's reserve the space
	dissection_list.reserve(subdomains.size());

	// catch case of no subdomains
	if (subdomains.size() != 0) {
		const piecewise_linear_function& low_h = h.get_lower();
		const piecewise_linear_function& up_h = h.get_upper();

		std::vector<piecewise_linear_function> lower_cuts = dissect(low_h,
				subdomains);
		std::vector<piecewise_linear_function> upper_cuts = dissect(up_h,
				subdomains);

		assert(lower_cuts.size()==upper_cuts.size());

		/** All the sub-domains in a loop */
		std::vector<piecewise_linear_function>::const_iterator low_it =
				lower_cuts.begin();
		std::vector<piecewise_linear_function>::const_iterator upp_it =
				upper_cuts.begin();
		while (low_it != lower_cuts.end()) {
			piecewise_linear_interval_function h_piece(*low_it, *upp_it);
			dissection_list.push_back(h_piece);
			++low_it;
			++upp_it;
		}
	}

	return dissection_list;
}

/** Mirror f along the y axis
 *
 * @author Goran Frehse */
inline
piecewise_linear_interval_function x_mirror(const piecewise_linear_interval_function& f) {
	return piecewise_linear_interval_function(x_mirror(f.get_lower()),x_mirror(f.get_upper()));
}

/** Splits h into a vector of nonempty, connected pieces
 *
 * @author Goran Frehse
 */

inline std::vector<piecewise_linear_interval_function> split(
		const piecewise_linear_interval_function& h) {

	piecewise_linear_function e = h.get_upper() - h.get_lower();

//	std::cout << "e:" << e << std::endl;

	std::vector<interval> nonempty_intervals = above_threshold_subdomains(e,
			precision_type(0));

//	std::cout << "intervals:" << nonempty_intervals << std::endl;

	std::vector<piecewise_linear_interval_function> res;
	res = dissect(h, nonempty_intervals);

//	std::cout << "dissected:" << res << std::endl;

	return res;
}
;

/** Returns a map from time points in s to the earliest time point in r where s<=r
 *
 * Both r and s are assumed to be continuous, and r is assumed to be unimodal.
 *
 * @author Goran Frehse */
inline piecewise_linear_function compute_first_inclusion_points(
		const piecewise_linear_function& s, const piecewise_linear_function& r,
		const precision_type& beyond_value) {
	typedef piecewise_linear_function::list_of_points_type list_type;

	list_type lpts, upts;

	const list_type& spts = s.get_list();
	const list_type& rpts = r.get_list();

	//std::cout << "Computing inclusion points from domain " << s.get_domain() << "(" << spts.size() << ") to " << r.get_domain() << "(" << rpts.size() << ")" << std::endl;
	//std::cout << spts << std::endl;
	//precision_type r_max = rit->get_y(); // largest value of r_y = last_ry

	// loop over points in s
	list_type::const_iterator rit = rpts.begin();

	precision_type last_rx = rit->get_x();
	precision_type last_ry = rit->get_y(); // to detect negative slope
	list_type::const_iterator max_rit = rit; // argmax r
	precision_type max_rx = rit->get_x(); // argmax r
	precision_type max_ry = rit->get_y(); // max r

	precision_type last_sx = spts.begin()->get_x(); // to detect negative slope
	precision_type last_sy = spts.begin()->get_y(); // to detect negative slope

	bool is_beyond = false;
	tribool intersects, dummy;

	for (list_type::const_iterator sit = spts.begin(); sit != spts.end();
			++sit) {
		const precision_type& sx = sit->get_x();
		const precision_type& sy = sit->get_y();
		// PART I: LOWER BOUND
		// if sy is increasing, search in increasing rx direction
		if (sy >= last_sy) {
//			std::cout << "searching forward for sy=" << sy << " at "
//					<< rit->get_x() << std::endl;
			// find the first point in r that is larger or equal
			// iterate until either we're at the end or ry >= sy or ry is decreasing
			while (rit != rpts.end() && definitely_is_LT(rit->get_y(),sy)
					&& !definitely_is_LT(rit->get_y(),max_ry)) {
				// if we're already beyond last_sy, add the interpolation points
				if (sit != spts.begin() && definitely_is_GT(rit->get_y(),last_sy)) {
					precision_type sx_interp = crosses_threshold(intersects,
							*(sit - 1), *sit, rit->get_y());
//					std::cout << "interpolating intermediate r between "
//							<< *(sit - 1) << " and " << *sit << " at y="
//							<< rit->get_y() << ", got " << std::boolalpha
//							<< dummy << " at " << sx_interp << std::endl;
					if (intersects)
					if (lpts.empty() || definitely_is_GT(sx_interp,lpts.back().get_x())) {
						lpts.push_back(breakpoint(sx_interp, rit->get_x()));
					}
				}

				last_rx = rit->get_x();
				last_ry = rit->get_y();
				if (max_ry < rit->get_y()) {
					max_ry = rit->get_y();
					max_rx = rit->get_x();
					max_rit = rit;
				}
				++rit;
			}
			if (rit == rpts.end())
				--rit;
			// if we found a dominating ry, insert result
			if (rit != rpts.end() && !definitely_is_LT(rit->get_y(),sy)) {
				// get a precise rx value through interpolation
				precision_type rx = rit->get_x();
				if (rit != rpts.begin()) {
					precision_type rinters = crosses_threshold(intersects, *(rit - 1), *rit, sy);
					if (intersects)
						rx =  rinters;
				}
				if (lpts.empty() || definitely_is_GT(sx,lpts.back().get_x()))
					lpts.push_back(breakpoint(sx, rx));
			} else {
				// if the previous value of sy was ok, find the interpolation point
				if (!definitely_is_GT(last_sy,max_ry)) {
					// couldn't find a matching r
					// find sx_interp where s reaches the max of r
					precision_type sx_interp = last_sx;
					if (sit != spts.begin()
							&& !definitely_is_GT((sit - 1)->get_y_right(),max_ry)) {
						precision_type sx_intersect = crosses_threshold(intersects, *(sit - 1), *sit,
								max_ry);
						if (intersects)
							sx_interp = sx_intersect;
					}
					// do something to mark the empty set...
					// Here: let's set the lower bound to past the highest value
					// Only insert another point if it's to the right of the last one
					if (lpts.empty() || definitely_is_GT(sx_interp,lpts.back().get_x())) {
						lpts.push_back(
								breakpoint(sx_interp, max_rx, max_rx,
										beyond_value));
					} else {
						// otherwise just redefine the right limit
						lpts.back().set_y_right(beyond_value);
					}
//					std::cout << "found exit point for ry=" << max_ry
//							<< " at sx=" << sx_interp << ", rx=" << max_rx
//							<< std::endl;
				} else {
					// insert at least one point to mark the beginning of the domain
					if (lpts.empty()) {
						lpts.push_back(breakpoint(sx, beyond_value));
					}
				}
				is_beyond = true;
			}
		} else {
//			std::cout << "searching backward for sy=" << sy << " at rx="
//					<< rit->get_x() << std::endl;
			// if sy is decreasing, search in increasing rx direction

			if (is_beyond) {
				// don't look backward, just wait for sy to go back down
				if (!definitely_is_LT(max_ry,sy)) {
					// it went back down, so get where it crossed the max
					precision_type sx_interp = sx;
					precision_type sx_intersect = crosses_threshold(intersects,
							*(sit - 1), *sit, max_ry);
					if (intersects) {
						sx_intersect = sx_interp;
					}
//					std::cout << "interpolating between " << *(sit - 1)
//							<< " and " << *sit << " at y=" << max_ry << ", got "
//							<< std::boolalpha << dummy << " at " << sx_interp
//							<< std::endl;
					if (lpts.empty() || definitely_is_GT(sx_interp,lpts.back().get_x())) {
						lpts.push_back(
								breakpoint(sx_interp, beyond_value, max_rx,
										max_rx));
					}
					rit = max_rit;
					last_rx = rit->get_x();
					last_ry = rit->get_y();
					is_beyond = false;
//					std::cout << "found re-entry point for ry=" << max_ry
//							<< " at sx=" << sx_interp << ", rx=" << max_rx
//							<< std::endl;
				}
			}
			if (!is_beyond) {
				// find the first point in r that is smaller or equal
				// iterate until either we're at the beginning or ry <= sy
				while (rit != rpts.begin() && definitely_is_GT(rit->get_y(),sy)) {
					// if we're already below last_sy, add the interpolation points
					if (sit != spts.begin() && definitely_is_LT(rit->get_y(),last_sy)) {
						precision_type sx_interp = last_sx;
						precision_type sx_intersect = crosses_threshold(intersects,
								*(sit - 1), *sit, rit->get_y());
						if (intersects)
							sx_interp = sx_intersect;
//						std::cout << "interpolating intermediate r between "
//								<< *(sit - 1) << " and " << *sit << " at y="
//								<< rit->get_y() << ", got " << std::boolalpha
//								<< dummy << " at " << sx_interp << std::endl;
						if (lpts.empty() || definitely_is_GT(sx_interp,lpts.back().get_x())) {
							lpts.push_back(breakpoint(sx_interp, rit->get_x()));
						}
					}

					last_rx = rit->get_x();
					last_ry = rit->get_y();
					--rit;
				}
				// if we found a minoring ry, insert result
				precision_type rx = rit->get_x();
				if (!definitely_is_GT(rit->get_y(),sy)) {
					// if point to the left was above, interpolate
					if ((rit + 1) != rpts.end() && !definitely_is_LT((rit + 1)->get_y(),sy)) {
						// get a precise rx value through interpolation
						rx = crosses_threshold(dummy, *rit, *(rit + 1), sy);
//						std::cout << "found interpolation point " << dummy
//								<< " for sy=" << sy << " at sx=" << sx
//								<< ", rx=" << rx << std::endl;
					}
				}
//				std::cout << "found a point inside point for sy=" << sy
//						<< " at sx=" << sx << ", rx=" << rx << std::endl;
				if (lpts.empty() || definitely_is_GT(sx,lpts.back().get_x())) {
					lpts.push_back(breakpoint(sx, rx));
				}
			}
		}
		last_sx = sx;
		last_sy = sy;
	}

//	if (lpts.front().get_x()>s.get_domain().lower()) {
//		std::cout << s << r << lpts << std::endl;
//	}
	if (definitely_is_LT(lpts.back().get_x(),s.get_domain().upper())) {
		lpts.push_back(breakpoint(s.get_domain().upper(), beyond_value));
	}

	piecewise_linear_function tl =
			piecewise_linear_function::create_closed_from_temp(lpts);
	return tl;
}
;

/** Returns a interval-valued map from time points in s to time points in r where s<=r
 *
 * Both s and r are assumed to be unimodal (otherwise there could be several intervals for
 * each time point) and continuous.
 *
 * @author Goran Frehse */
inline piecewise_linear_interval_function compute_inclusion_map(
		const piecewise_linear_function& s,
		const piecewise_linear_function& r) {
	precision_type tr_low = r.get_domain().lower();
	precision_type tr_upp = r.get_domain().upper();
	precision_type beyond_upp = tr_upp + (tr_upp - tr_low)
			+ precision_type(1.0); // to mark beyond upper bound
	precision_type beyond_low = tr_low - (tr_upp - tr_low)
			- precision_type(1.0); // to mark beyond upper bound

	piecewise_linear_function tl = compute_first_inclusion_points(s, r,
			beyond_upp);
	// to get the upper bound (latest points),
	// use r(-x), then negate the resulting function
	piecewise_linear_function tu = -compute_first_inclusion_points(s,
			x_mirror(r), -beyond_low);
	return piecewise_linear_interval_function(tl, tu);
}
;

/** Simplify the lower and upper bounds */
inline piecewise_linear_interval_function simplify(
		const piecewise_linear_interval_function& s) {
	piecewise_linear_function f = simplify(s.get_lower());
	piecewise_linear_function g = simplify(s.get_upper());
	return piecewise_linear_interval_function(f,g);
}
;

}

#endif /* PIECEWISE_LINEAR_INTERVAL_FUNCTION_OPERATORS_H_ */
