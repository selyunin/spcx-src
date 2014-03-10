/*
 * piecewise_linear_function_operators.h
 *
 *  Created on: Oct 18, 2012
 *      Author: kateja
 */

#ifndef PIECEWISE_LINEAR_FUNCTION_OPERATORS_H_
#define PIECEWISE_LINEAR_FUNCTION_OPERATORS_H_

#include "piecewise_linear_function.h"
#include "geometry.h"
#include "numeric_comp.h"
#include "misc_operators.h"

#include <iostream>
#include <iterator>
#include <algorithm>

namespace plif {

/**
 * Scalar Multiplication.
 *
 * Keeps domain same, multiplies left and right slope by lambda, and the y value and left and right limit values
 * at each breakpoint by lambda
 */
inline
piecewise_linear_function operator*(const piecewise_linear_function& f, const precision_type& lambda)
{
	piecewise_linear_function lambda_f = f;
	lambda_f *= lambda;
	return lambda_f;
}

/**
 * Scalar Multiplication.
 *
 * Keeps domain same, multiplies left and right slope by lambda, and the y value and left and right limit values
 * at each breakpoint by lambda
 */
inline
piecewise_linear_function operator*(const precision_type& lambda, const piecewise_linear_function& f)
{
	return f*lambda;
}

/**
 * Scalar Addition.
 *
 * Keeps domain and left and right slops the same, adds lambda to the y value and left and right limit values
 * at each breakpoint.
 */
inline
piecewise_linear_function operator+(const piecewise_linear_function& f, const precision_type& lambda)
{
	piecewise_linear_function lambda_f = f;
	lambda_f += lambda;
	return lambda_f;
}

/**
* Overloading the negative operator for PLF's.
*
* Returns the negative of the PLF, i.e. -f
*/
inline
piecewise_linear_function operator-(const piecewise_linear_function& f)
{
	piecewise_linear_function g = f*precision_type(-1);
	g.set_computation_parameters(negate(f.get_computation_parameters()));
	return g;
}

/**
 * Infimum of a PLF function.
 *
 * Checks for the infimum at the end points and returns that as the infimum of the function
 * in case of unbounded domain with megative slopes, returns negative infinity
 *
 */
inline precision_type infimum(const piecewise_linear_function& func) {
	if (func.get_list().size() <= 0) {
		throw std::runtime_error("calling infimum on plf without breakpoints");
	}

	if ((func.is_left_unbounded() && func.get_left_slope() < 0)
			|| (func.is_right_unbounded() && func.get_right_slope() < 0)) {
		return NEG_INFTY;
	}

	precision_type ret_val;
	const std::vector<breakpoint>& list = func.get_list();
	piecewise_linear_function::list_of_points_type::const_iterator it =
			list.begin();

	// now we know that the infimum is finite
	// check special case in which we don't look at neighborhoods
	if (func.get_list().size() == 1) {
		if (func.is_left_closed()) {
			if (func.is_right_closed()) {
				return it->get_y();
			} else {
				return it->inf_right();
			}
		} else {
			if (func.is_right_closed()) {
				return it->inf_left();
			} else {
				return it->inf();
			}
		}
	}

	// now we know that there's at least two breakpoints

	// take the right inf of the first point if it's closed
	if (func.is_left_closed()) {
		ret_val = it->inf_right();
	} else {
		ret_val = it->inf();
	}
	++it;
	// take the left and right inf of the points in the interior
	if (it != list.end()) {
		piecewise_linear_function::list_of_points_type::const_iterator last_it =
				--list.end();
		while (it != last_it) {
			ret_val = std::min(ret_val, it->inf());
			++it;
		}
		// take the left inf of the last point if it's closed
		if (func.is_right_closed()) {
			ret_val = std::min(ret_val, it->inf_left());
		} else {
			ret_val = std::min(ret_val, it->inf());
		}
	}

	return ret_val;
}

/**
 * Supremum of a PLF.
 *
 * @author Goran Frehse
 */
inline
precision_type supremum(const piecewise_linear_function& func) {
	if (func.get_list().size() <= 0) {
		throw std::runtime_error("calling supremum on plf without breakpoints");
	}

	if ((func.is_left_unbounded() && func.get_left_slope() < 0)
			|| (func.is_right_unbounded() && func.get_right_slope() < 0)) {
		return POS_INFTY;
	}

	precision_type ret_val;
	const std::vector<breakpoint>& list = func.get_list();
	piecewise_linear_function::list_of_points_type::const_iterator it =
			list.begin();

	// now we know that the infimum is finite
	// check special case in which we don't look at neighborhoods
	if (func.get_list().size() == 1) {
		if (func.is_left_closed()) {
			if (func.is_right_closed()) {
				return it->get_y();
			} else {
				return it->sup_right();
			}
		} else {
			if (func.is_right_closed()) {
				return it->sup_left();
			} else {
				return it->sup();
			}
		}
	}

	// now we know that there's at least two breakpoints

	// take the right sup of the first point if it's closed
	if (func.is_left_closed()) {
		ret_val = it->sup_right();
	} else {
		ret_val = it->sup();
	}
	++it;
	// take the left and right sup of the points in the interior
	if (it != list.end()) {
		piecewise_linear_function::list_of_points_type::const_iterator last_it =
				--list.end();
		while (it != last_it) {
			ret_val = std::max(ret_val, it->sup());
			++it;
		}
		// take the left sup of the last point if it's closed
		if (func.is_right_closed()) {
			ret_val = std::max(ret_val, it->sup_left());
		} else {
			ret_val = std::max(ret_val, it->sup());
		}
	}

	return ret_val;
}

/** Combines the y values of p and q by applying the function bin_op
 *
 * The elements y, y_left, y_right of p are modified using
 * p.y = bin_op(p.y,q.y), etc.
 *
 * @author Goran Frehse */
typedef precision_type (*scalar_bin_op) (const precision_type&, const precision_type&);
inline
void apply_bin_op(breakpoint& p,const breakpoint& q, scalar_bin_op bin_op) {
	p.get_y_left() = bin_op(p.get_y_left(),q.get_y_left());
	p.get_y() = bin_op(p.get_y(),q.get_y());
	p.get_y_right() = bin_op(p.get_y_right(),q.get_y_right());
}

/**
 *
 * Combines the breakpoints of F and G by applying bin_op to their y-values.
 * Missing x-values are interpolated.
 * 
 * both have bounded domains
 * addition is performed on the intersection of the two functions domain
 *
 * @note The breakpoints consist of the union of the breakpoints of f and g.
 *
 * @author Complete reimplementation by Goran Frehse
 */
inline
piecewise_linear_function breakpoint_bin_op(const piecewise_linear_function& f, const piecewise_linear_function& g, scalar_bin_op bin_op )
{
	const computation_parameters& f_params = f.get_computation_parameters();
	const computation_parameters& g_params = g.get_computation_parameters();

	if (f.is_left_unbounded() || f.is_right_unbounded() || g.is_left_unbounded() || g.is_right_unbounded()  ) {
		throw std::runtime_error("breakpoint_bin_op not implemented for unbounded functions");
	}
	/*
	if (f.get_list().size()<2 || g.get_list().size()<2 ) {
		// acceptable if the domain is a singleton
		const precision_type eps = 1e-11;
		if (f.get_domain().upper()-f.get_domain().lower()<eps && g.get_domain().upper()-g.get_domain().lower()<eps) {
			std::vector<breakpoint> pts;
			breakpoint p = f.get_list().front();
			apply_bin_op(p,g.get_list().front(),bin_op);
			pts.push_back(p);
			return piecewise_linear_function::create_closed_from_temp(pts);
		} else {
			f.display(std::cerr); g.display(std::cerr);
			throw std::runtime_error("breakpoint_bin_op doesn't know how to interpolate on lists with less than 2 points");
		}
	}
	*/
//	piecewise_linear_function sum;
	interval domain;
	domain = intersect(f.get_domain(), g.get_domain());
//	sum.set_domain(domain);

	const std::vector<breakpoint>& f_list=f.get_list();
	const std::vector<breakpoint>& g_list=g.get_list();

	std::vector<breakpoint>::const_iterator it=f_list.begin();
	std::vector<breakpoint>::const_iterator jt=g_list.begin();

//  std::cout << std::endl << "adding " << f_list << " and " << g_list << std::endl;

	std::vector<breakpoint> pts;
	// reserve space for at least as many points as in f or g (the result could have less)
	pts.reserve(f_list.size() + g_list.size());
	breakpoint p;

	// find lowest element in domain
	breakpoint bsta(domain.lower(), precision_type(0));
	it = std::lower_bound(it, f_list.end(), bsta);
	jt = std::lower_bound(jt, g_list.end(), bsta);

//	std::cout << "starting from " << *it << " and " << *jt << std::endl;

	while ((it->get_x() <= domain.upper() || jt->get_x() <= domain.upper())
			&& it != f_list.end() && jt != g_list.end()) {
		// the current point is the left one of both
		if ((it->get_x()<jt->get_x())) {
			// interpolate g on f's position
			p = *it;
			breakpoint q = interpolate(*(jt - 1), *jt, it->get_x(),g_params);
			apply_bin_op(p,q,bin_op);
//			std::cout << "interpolated " << *(jt-1) << " and " << *(jt)
//					<< " at " << it->get_x() << " to yield " << p
//					<< std::endl;
			pts.push_back(p);
			++it;
		} else if ((jt->get_x()<it->get_x())) {
			// interpolate f on g's position
			p = interpolate(*(it - 1), *it, jt->get_x(),f_params);
			apply_bin_op(p,*jt,bin_op);
//			std::cout << "interpolated " << *(it-1) << " and " << *(it)
//					<< " at " << jt->get_x() << " to yield " << p
//					<< std::endl;
			pts.push_back(p);
			++jt;
		} else {
			p = *it;
			apply_bin_op(p,*jt,bin_op);
			pts.push_back(p);
//			std::cout << "combined " << *it << " and " << *jt << " to yield " << p
//					<< std::endl;
			++it;
			++jt;
		}
	}

//	std::cout << "last point:" << *pts.rbegin() << std::endl;

	/** Make sure we didn't loose any points */
	if (pts.size()<std::min(f.get_list().size(),g.get_list().size())) {
		std::cerr << "missing points: " <<
				pts.size() << " < " << f.get_list().size() << " and " << g.get_list().size() << std::endl;
		throw std::runtime_error(__FUNCTION__);
	}
//	assert(pts.size()>=std::min(f.get_list().size(),g.get_list().size()));

	//computation_parameters res_params = combine(f_params,g_params);
	return piecewise_linear_function::create_closed_from_temp(pts,f_params);
}

/** Function wrappers for built-in functions
 *
 * Needed since there are no function pointers to built-in types.
 * @author Goran Frehse*/
inline
precision_type sum(const precision_type& a, const precision_type& b) {
	return a+b;
}
inline
precision_type diff(const precision_type& a, const precision_type& b) {
	return a-b;
}
inline
precision_type max(const precision_type& a, const precision_type& b) {
	return std::max(a,b);
}
inline
precision_type min(const precision_type& a, const precision_type& b) {
	return std::min(a,b);
}

/** Pointwise sum of two closed and bounded PLF
 *
 * @author Goran Frehse
 */
inline
piecewise_linear_function operator+(const piecewise_linear_function& f, const piecewise_linear_function& g)
{
	return breakpoint_bin_op(f,g,sum);
}

/** Pointwise max of two closed and bounded PLF
 *
 * @author Goran Frehse
 */
inline
piecewise_linear_function pointwise_maximum(const piecewise_linear_function& f, const piecewise_linear_function& g)
{
	return breakpoint_bin_op(f,g,max);
}

/** Pointwise max of two closed and bounded PLF
 *
 * @author Goran Frehse
 */
inline
piecewise_linear_function pointwise_minimum(const piecewise_linear_function& f, const piecewise_linear_function& g)
{
	return breakpoint_bin_op(f,g,min);
}

/**
 *
 * Finds the pointwise subtracted PLF of the two PLF's F and G
 *
 * both have bounded domains
 * substraction is performed on the intersection of the two functions domain
 *
 */
inline
piecewise_linear_function operator-(const piecewise_linear_function& f, const piecewise_linear_function& g)
{
	// return breakpoint_bin_op(f,g,diff);
	/** Note: the rounding needs to be taken into account, so use negation */
	return f+(-g);
}

/** Distance function: returns the minimum (over range x) of the PLF given as F-G */
inline
precision_type dist(const piecewise_linear_function& f, const piecewise_linear_function& g)
{
	return infimum(g-f);
}

/** Convert a list of cut points to a list of intervals on the domain intv
 *
 * @author Goran Frehse */
inline std::vector<interval> cut_points_to_intervals(
		const std::vector<precision_type> cut_points, const interval& intv) {
	std::vector<interval> intv_vec;
	intv_vec.reserve(cut_points.size() + 2); // can't be more than that

	precision_type last_x = intv.lower(); // first interval starts at domain
	for (std::vector<precision_type>::const_iterator it = cut_points.begin();
			it != cut_points.end() && *it < intv.upper(); ++it) {
		if (*it > last_x) {
			intv_vec.push_back(interval(last_x, *it));
			last_x = *it;
		}
	}
	// last interval ends at domain
	intv_vec.push_back(interval(last_x, intv.upper()));

	return intv_vec;
}

/** Divides a PLF f into subdomains only overlapping at a given set of cut points.
  *
  * It is assumed that all cut points lie within the domain of f.
  * The domain of the first piece is the start of the domain of f to the
  * first cut point , the second from the first cut point to
  * the second cut point, etc. The domain of the last piece is the last
  * cut point to the end of the domain of f.
  *
  * @ The breakpoints consist of the breakpoints of f union the cut points.
  *
  * @author Rajat Kateja, Goran Frehse
  */
inline
std::vector<piecewise_linear_function> dissect(const piecewise_linear_function& f, const std::vector<precision_type> cut_points)
{
	assert( !f.is_left_unbounded() && !f.is_right_unbounded());

	const computation_parameters& f_params = f.get_computation_parameters();

	typedef piecewise_linear_function::list_of_points_type list_of_points_type;
	std::vector<piecewise_linear_function> dissection_list;
	dissection_list.reserve(cut_points.size() + 1);

//std::cout << "dissect on "; f.display();

	// catch case of no cutpoints
	if (cut_points.size() == 0) {
		dissection_list.push_back(f);
	} else {
		precision_type start, end;
		piecewise_linear_function cut_function;
		int i;

		const list_of_points_type& list = f.get_list();
		list_of_points_type::const_iterator bit,eit;
		bit = list.begin();
		eit = bit;

		list_of_points_type list_of_points;
		list_of_points.reserve(list.size());
		breakpoint bp,ep;
		/** All the sub-domains in a loop */
		// compute whether there will be an additional point
		// at the beginning and the end
		int i_start = 0;
		int i_end = cut_points.size();
		// for now: always return cut_points.size()+1 pieces (required by spacetime_plif)
		//if (*cut_points.begin()>f.get_domain().lower())
			i_start = -1;
		//if (*(cut_points.end()-1)<f.get_domain().upper())
			i_end = cut_points.size();
		for (i = i_start; i < i_end; i++) {
			if (i >= 0) {
				start = cut_points[i];
			} else {
				start = f.get_domain().lower();
			}
			if (i + 1 < cut_points.size()) {
				end = cut_points[i + 1];
			} else {
				end = f.get_domain().upper();
			}

			// find lowest element in domain
			breakpoint bsta(start, precision_type(0));
			breakpoint bend(end, precision_type(0));
			bit = std::lower_bound(
					eit, list.end(), bsta);
			eit = std::upper_bound(bit,
					list.end(), bend);
			// if bit doesn't point to the exact value, interpolate from the previous
			// we need to know the number of points before we can allocate the vector
			// so we compute this beforehand
			size_t lsize = eit-bit; // size of the list
			bool use_bp = false;
			if (bit>list.begin() && start < bit->get_x()) {
				++lsize;
				bp = interpolate(*(bit - 1), *bit, start, f_params);
				use_bp = true;
			}
			// if eit-1 doesn't point to the exact value, interpolate from the previous
			bool use_ep = false;
			if (eit<list.end() && (eit-1)->get_x() < end) {
				++lsize;
				ep = interpolate(*(eit - 1), *eit, end, f_params);
				use_ep = true;
			}
			// allocate the vector and fill with values
			list_of_points = list_of_points_type(lsize);
			list_of_points_type::iterator lit = list_of_points.begin();
			if (use_bp) {
				*lit = bp;
				++lit;
			}
			std::copy(bit, eit, lit);
			if (use_ep) {
				*list_of_points.rbegin() = ep;
			}

			cut_function = piecewise_linear_function::create_closed_from_temp(list_of_points,f_params);
			dissection_list.push_back(cut_function);

//std::cout << "cutting from " << start << " to " << end << std::boolalpha << " bp:"<<use_bp << " ep:" << use_ep << std::endl;
//cut_function.display();

			if (eit != list.begin()) {
				--eit; // this is the position of the next beginning; it could start at the same point
			}
		}
	}

	return dissection_list;
}

/** Divides a PLF f into subdomains given by a set of intervals with disjoint interiors
  *
  * @note Resulting breakpoints are a subset of the breakpoints of f union the start and
  * end points of the subdomains.
  *
  * @author Goran Frehse
  */
inline
std::vector<piecewise_linear_function> dissect(
		const piecewise_linear_function& f,
		const std::vector<interval> subdomains) {
	assert( !f.is_left_unbounded()
			&& !f.is_right_unbounded() );

	std::vector<piecewise_linear_function> plf_result;

	if (subdomains.size() > 0) {
		precision_type start, end;
		piecewise_linear_function cut_function;
		int i;

		/** All the sub-domains in a loop */
		for (i = 0; i < subdomains.size(); i++) {
			start = subdomains[i].lower();
			end = subdomains[i].upper();
			if (start < f.get_domain().lower() || end > f.get_domain().upper()) {
				throw std::runtime_error("cannot dissect plf outside of domain");
			}
			cut_function.set_domain(interval(start, end));
			cut_function.set_left_closed(true);
			cut_function.set_right_closed(true);
			cut_function.set_left_bounded();
			cut_function.set_right_bounded();
			cut_function.set_computation_parameters(f.get_computation_parameters());
			cut_function.copy_on_domain(f, start, end);

			plf_result.push_back(cut_function);
		}
	}

	return plf_result;
}

/** Restrict a PLF f into a subdomain
  *
  * @author Goran Frehse
  */
inline
piecewise_linear_function restrict(const piecewise_linear_function& f,
		const interval& subdomain) {
	if (subset(f.get_domain(), subdomain)) {
		return f;
	} else {
		std::vector<interval> subvec;
		subvec.push_back(subdomain);
		std::vector<piecewise_linear_function> avec = dissect(f, subvec);
		assert(avec.size()==1);
		return *avec.begin();
	}
}

/** Mirror f along the y axis
 *
 * @author Goran Frehse */
inline
piecewise_linear_function x_mirror(const piecewise_linear_function& f) {
	piecewise_linear_function g;
	g.set_computation_parameters(f.get_computation_parameters());
	g.set_domain(-f.get_domain());
	if (f.is_left_unbounded())
		g.set_right_unbounded(-f.get_left_slope());
	else
		g.set_right_bounded();
	if (f.is_right_unbounded())
		g.set_left_unbounded(-f.get_right_slope());
	else
		g.set_left_bounded();

	g.set_right_closed(f.is_left_closed());
	g.set_left_closed(f.is_right_closed());

	piecewise_linear_function::list_of_points_type g_points(f.get_list().size());
	// careful: we also need to reverse the order
	std::transform(f.get_list().rbegin(),f.get_list().rend(),g_points.begin(),breakpoint_x_mirror);
	g.set_list_from_temp(g_points);
	return g;
}

/** Standard stream output operator
 *
 * @author Goran Frehse */
inline
std::ostream& operator<<(std::ostream& os, const piecewise_linear_function& f) {
	f.display(os);
	return os;
}

/** Returns how many consecutive points in the function are numerically indistinguishable in time
 */
inline
unsigned int indistinguishable_count(const piecewise_linear_function& f) {
	const piecewise_linear_function::list_of_points_type& lop = f.get_list();

	unsigned int warning_superimp_count = 0;

	for (piecewise_linear_function::list_of_points_type::const_iterator it =
			lop.begin(); (it != lop.end()) && ((it + 1) != lop.end()); ++it) {

		const precision_type& x1 = it->get_x();
		const precision_type& x2 = (it + 1)->get_x();
		// if for some reason the x values are the same don't do anything
		if (!definitely_is_LT(x1, x2)) {
			std::cout << std::endl << "offended by: " << x1 << " !< " << x2 << std::endl;
			++warning_superimp_count;
		}
	}

	return warning_superimp_count;
}

/** Removes indistinguishable points from f
 *
 * Rounding is done conservatively according to the rounding mode.
 *
 * @note Does not remove the first or the last breakpoint of f,
 * so bounded plf keep their domain.
 *
 * @todo Handle unbounded plfs.
 */
inline
piecewise_linear_function simplify(const piecewise_linear_function& f) {
	assert(!f.is_left_unbounded() && !f.is_right_unbounded());

	typedef piecewise_linear_function::list_of_points_type list;
	const list& lop = f.get_list();

	unsigned int warning_superimp_count = 0;
	list pts;
	pts.reserve(lop.size());

	// take the first point
	if (lop.size() >= 1) {
		pts.push_back(*lop.begin());
		if (lop.size() > 1) {
			list::const_iterator it = ++lop.begin();
			list::iterator pit = pts.begin();
			while (it != lop.end()) {
				const precision_type& xlast = pit->get_x();
				const precision_type& x = it->get_x();
				// if for some reason the x values are the same don't do anything
				if (!definitely_is_LT(xlast, x)) {
					// merge the two into one
					++warning_superimp_count;
					// new x is left of xlast, max of xlast.y,xlast.yr,x.yl,x.y
					pit->set_y(
							max(pit->get_y(), it->get_y(), it->get_y_left()));
					pit->set_y_right(
							max(pit->get_y_right(), it->get_y_right()));
					// if it's the last point, use the last (rightmost) x
					if ((it+1) == lop.end()) {
						pit->set_x(max(pit->get_x(),it->get_x()));
					}
				} else {
					pts.push_back(*it);
					++pit;
				}
				++it;
			}
		}
	}

	piecewise_linear_function res = piecewise_linear_function::create_closed_from_temp(pts,
			f.get_computation_parameters());
	if (warning_superimp_count > 0) {
//		std::cout << "removed " << warning_superimp_count
//				<< " points in simplification" << std::endl;
//		std::cout << "from "; f.display(std::cout); std::cout << "to "; res.display(std::cout);
	}

	return res;
}
;

/** Returns a list of intervals in the domain on which the PLF f satisfies f(x)>=c for a given threshold c.
 *
 * The intervals are disjoint and x is in exactly one of them if f(x)>=c.
 */

inline std::vector<interval> above_threshold_subdomains(
		const piecewise_linear_function& f, const precision_type& c) {
	/*
	 assert(!f.is_left_unbounded() && !f.is_right_unbounded());
	 interval domain = f.get_domain();
	 piecewise_linear_function upper = f;
	 piecewise_linear_function lower;
	 lower.set_domain(domain);
	 lower.set_left_bounded();
	 lower.set_right_bounded();
	 lower.set_left_closed(f.is_left_closed());
	 lower.set_right_closed(f.is_right_closed());
	 lower.insert_continuous_point(domain.lower(), c);
	 lower.insert_continuous_point(domain.upper(), c);
	 piecewise_linear_interval_function g(lower, upper);
	 std::vector<piecewise_linear_interval_function> list_of_plif = split(g);
	 std::vector<interval> subdomains;
	 subdomains.clear();
	 int i;
	 for(i=0;i<list_of_plif.size();i++)
	 {
	 interval subdomain;
	 subdomain = list_of_plif[i].get_domain();
	 subdomains.push_back(subdomain);
	 }
	 return subdomains;
	 */

	// Walk through the breakpoints
	std::vector<interval> subdomains;
	typedef piecewise_linear_function::list_of_points_type list_type;
	const list_type& list = f.get_list();
	if (list.empty())
		return subdomains;
	bool intv_open = false;
	precision_type a = -1.0;
	precision_type b = -1.0;
	list_type::const_iterator it = list.begin();

	for (; it != list.end(); ++it) {
//std::cout << "at " << it->get_x() << ":" << it->get_y() << ">=" << c << std::boolalpha << ":" << intv_open << std::endl;
		// check if we need to close at current c (and possibly reopen)
		if (intv_open) {
			if (definitely_is_LT(it->get_y(),c)) {
				b = it->get_x();
				subdomains.push_back(interval(a, b));
				intv_open = false;
				if (!definitely_is_LT(it->get_y_right(),c)) {
					a = it->get_x();
					intv_open = true;
				}
			}
		} else {
			if (!definitely_is_LT(it->get_y(),c)) {
				a = it->get_x();
				intv_open = true;
				if (definitely_is_LT(it->get_y_right(),c)) {
					b = it->get_x();
					subdomains.push_back(interval(a, b));
					intv_open = false;
				}
			}
		}

		if (intv_open) {
			// try to find the closing point
			if (it + 1 != list.end() && definitely_is_LT((it + 1)->get_y_left(),c)) {
				// @note: y2-y1 can't be zero, otherwise the interval wouldn't be open
				// find the intersection point and assign it to b
				tribool dummy;
				b = crosses_threshold(dummy, *it, *(it + 1), c);
				subdomains.push_back(interval(a, b));
				intv_open = false;
			}
		} else {
			if (it + 1 != list.end() && !definitely_is_LT((it + 1)->get_y_left(),c)) {
				// @note: y2-y1 can't be zero, otherwise the interval wouldn't be open
				// find the intersection point and assign it to a
				tribool dummy;
				a = crosses_threshold(dummy, *it, *(it + 1), c);
				intv_open = true;
			}
		}
	}
	if (intv_open) {
		// we're at the end, let's close
		b = list.back().get_x();
		subdomains.push_back(interval(a, b));
	}

	return subdomains;
}

}

#endif /* PIECEWISE_LINEAR_FUNCTION_OPERATORS_H_ */
