/*
 * interval.h
 *
 *  Created on: Oct 19, 2012
 *      Author: kateja
 */

#ifndef INTERVAL_H_
#define INTERVAL_H_

#include <boost/numeric/interval.hpp>

#include "numeric_comp.h"

namespace plif {

/**
 * These two are symbols to be used for positive and negative infinity, using the definitions in C++ standard library
 */
#define NEG_INFTY (-1*std::numeric_limits<precision_type>::infinity())
#define POS_INFTY (std::numeric_limits<precision_type>::infinity())

typedef double precision_type;

/** Define an interval type that does proper rounding, disallows NaN, but accepts empty intervals.
 *
 * Note that by default the boost interval cannot be empty, and an exception is thrown if it is.
 */
typedef boost::numeric::interval_lib::rounded_math<precision_type> rounding_policy;
typedef boost::numeric::interval_lib::checking_no_nan<precision_type> checking_policy;
typedef boost::numeric::interval_lib::policies<rounding_policy,checking_policy> interval_policies;
//typedef boost::numeric::interval<precision_type,interval_policies> interval;

/** A class for representing time intervals */
class interval{
public:
	typedef precision_type scalar_type;
	interval() {
		my_lower = NEG_INFTY;
		my_upper = POS_INFTY;
	}
	explicit interval(scalar_type x) :
		my_lower(x), my_upper(x) {
	}
	interval(scalar_type l, scalar_type u) :
		my_lower(l), my_upper(u) {
		if (u<l && is_MEQ(l,u)) {
			my_lower=u;
			my_upper=l;
		}
	}
	scalar_type lower() const{
		return my_lower;
	}
	scalar_type upper() const{
		return my_upper;
	}
	void set_lower(scalar_type l){
		my_lower = l;
	}
	void set_upper(scalar_type u){
		my_upper = u;
	}
	void set(scalar_type l, scalar_type u) {
		*this = interval(l,u);
	}
	void set_empty(){
		my_lower = scalar_type(1);
		my_upper = scalar_type(0);
	}
	scalar_type width() const{
		return my_upper - my_lower;
	}
	bool is_empty() const{
		return definitely_is_LT(my_upper,my_lower);
	}
	bool is_finite() const {
		return my_lower!=NEG_INFTY && my_upper!=POS_INFTY;
	}
	/** Returns the negated interval (values mirrored around zero) */
	interval operator-() const {
		if (is_empty()) {
			return *this;
		} else {
			return interval(-my_upper, -my_lower);
		}
	}
	/** Returns the interval shifted up by a scalar */
	interval operator+(scalar_type x) const {
		if (!is_empty()) {
			return interval(my_lower + x, my_upper + x);
		} else {
			return *this;
		}
	}
	/** Returns the interval shifted down by a scalar */
	interval operator-(scalar_type x) const {
		if (!is_empty()) {
			return interval(my_lower - x, my_upper - x);
		} else {
			return *this;
		}
	}
	/** Returns the pointwise sum of *this with another interval */
	interval operator+(const interval& intv) const {
		if (!is_empty() && !intv.is_empty()) {
			return interval(my_lower+intv.my_lower,my_upper+intv.my_upper);
		} else {
			return empty();
		}
	}
	/** Assigns to *this the pointwise sum with an interval */
	interval& operator+=(const interval& intv) {
		if (!is_empty() && !intv.is_empty()) {
			*this = interval(my_lower+intv.my_lower,my_upper+intv.my_upper);
		} else {
			*this = empty();
		}
		return *this;
	}
	static interval empty() {
		interval intv;
		intv.set_empty();
		return intv;
	}
	static interval whole() {
		interval intv;
		return intv;
	}
private:
	bool is_finite(scalar_type x) {
		if (x==POS_INFTY || x==NEG_INFTY)
			return false;
		else
			return true;
	}
	scalar_type my_lower;
	scalar_type my_upper;
};

/** Returns true if the interval is empty.
 *
 * @note Included for interface compatibility with boost::numeric::interval. */
inline
bool empty(interval itv) {
	return itv.is_empty();
}

/** Returns the size of the interval.
 *
 * @note Included for interface compatibility with boost::numeric::interval. */
inline
interval::scalar_type width(interval itv) {
	return itv.width();
}

/** Returns true if the first interval is maybe contained in the second one.
 *
 * @note Included for interface compatibility with boost::numeric::interval. */
inline
bool subset(interval itv1, interval itv2) {
	if (definitely_is_LT(itv1.upper(),itv2.lower()) || definitely_is_LT(itv2.upper(),itv1.lower())) {
		return false;
	} else {
		return true;
	}
}

/** Returns the convex hull of two intervals.
 *
 * @note Included for interface compatibility with boost::numeric::interval. */
inline
interval hull(interval A, interval B) {
	if (A.is_empty())
		return B;
	else if (B.is_empty())
		return A;
	else {
		precision_type a = std::min(A.lower(), B.lower());
		precision_type b = std::max(A.upper(), B.upper());
		return interval(a, b);
	}
}


/**
 *Returns the interval as intersection of two intervals. Simply the maximum of the lower limits and minimum of the upper limits.
 */
inline
interval intersect(interval A, interval B)
{
//	return boost::numeric::intersect(A,B);
	/** Don't use boost intersection to allow for numerical errors */
	if (A.is_empty() || B.is_empty()) {
		return interval::empty();
	} else {
		precision_type a = std::max(A.lower(), B.lower());
		precision_type b = std::min(A.upper(), B.upper());
		if (definitely_is_LT(b, a)) {
			return interval::empty();
		} else if (a <= b) {
			return interval(a, b);
		} else {
			return interval(b, a);
		}
	}
}
/**
 * Returns whether x belongs to the interval A
 *
 * @note Included for interface compatibility with boost::numeric::interval. */
inline
bool in(precision_type x, interval A)
{
	if (A.is_empty())
		return false;
	else {
		return !(definitely_is_LT(x,A.lower()) || definitely_is_LT(A.upper(),x));
	}
}

/** 
 *Comparing operator for intervals
 */

inline
bool strict_weak_interval_order (interval A, interval B)
{
	return (A.upper() < B.lower());
} 

/** Functor for ordering indices
 * @author Goran Frehse
 */
class interval_vector_order {
private:
	const std::vector<interval>& my_itvs;
public:
	interval_vector_order(const std::vector<interval>& itvs) : my_itvs(itvs) {};
	bool operator()(size_t i,size_t j) {
		return strict_weak_interval_order(my_itvs[i],my_itvs[j]);
	};
};

/** 
 * 
 * Function for optimal partitioning over several functions
 *
 * Takes as argument a vector of intervals
 * Returns the vector of overlap intervals as well as a vector of indices for each
 * overlap interval to the original intervals.
 *
 */
 
inline
std::pair<std::vector<interval>,std::vector<std::vector<size_t> > > partition(const std::vector<interval>& I)
{
	/** Construction of comparability graph */
	interval_vector_order comp(I);

	size_t graph_size = I.size();
	typedef std::vector<bool> bool_vector;
	typedef std::vector<bool_vector> bool_matrix;
	bool_matrix graph(graph_size,bool_vector(graph_size));			/** graph[i][j] = true implies an edge from interval i to j, implying they dont overlap */
	size_t edge_count = 0;
	for(size_t i=0;i<graph_size;i++)
	{
		for(size_t j=0;j<graph_size;j++)
		{
			bool is_edge = comp(i,j);
			graph[i][j] = is_edge;
			if (is_edge) {
				++edge_count;
			}
		}
	}
//	std::cout << "Graph with " << graph_size << " vertices and " << edge_count << " edges." <<std::endl;
	
	// Create a sorted index vector for the greedy coloring
	std::vector<size_t> sorted_indices(graph_size);
	for (size_t i=0;i<graph_size;++i) {
		sorted_indices[i]=i;
	}
	std::sort(sorted_indices.begin(), sorted_indices.end(), comp);

	std::vector<bool> colour_code(2);
	colour_code[0] = true; // zero means uncolored
	colour_code[1] = false;
	int colour[graph_size];				/** Stores the colour of the vertex (i.e. the interval) */
	
	for(size_t i=0;i<graph_size;i++)
	{
		colour[i] = 0;
	}
	int max_c = 0;
	// proceed by order of the sorted indices
	for (size_t si = 0; si < graph_size; si++) {
		size_t i = sorted_indices[si];
		for(size_t j=0;j<graph_size;j++)
		{
			if(graph[i][j] || graph[j][i])
			{
				if(colour[j] != 0)
				{
					colour_code[colour[j]] = true;
				}
			}
		}
		int c = 1;
		while(colour_code[c])
		{
			c++;
			if (c>=colour_code.size())
				colour_code.push_back(false);
		}
		for(size_t j=0;j<graph_size;j++)
		{
			if(graph[i][j] || graph[j][i])
			{
				if(colour[j] != 0)
				{
					colour_code[colour[j]] = false;
				}		
			}
		}
		colour[i] = c;
		max_c = std::max(max_c,c);
	}
	
	
	std::vector<interval> partition_of_intervals(max_c);
//	bool only_one_colour = true;
//	for(i=0;i<graph_size && only_one_colour;i++)
//	{
//		if(colour[i] != colour[0])
//		{
//			only_one_colour = false;
//		}
//	}
	if(max_c==1)
	{
		interval J(NEG_INFTY, POS_INFTY);
		for(size_t i=0;i<graph_size;i++)
		{
			J = intersect(J, I[i]);
		}
		partition_of_intervals[max_c-1]=J;
	}
	else
	{
		bool marked[graph_size];
		for(size_t i=0;i<graph_size;i++)
		{
			marked[i] = false;	
		}
		for(size_t i=0;i<graph_size;i++)
		{
			if(!marked[i])
			{
				interval J(NEG_INFTY, POS_INFTY);
				for(size_t j=0;j<graph_size;j++)
				{
					if(colour[i] == colour[j])
					{
						J = intersect(J, I[j]);
						marked[j] = true;
					}
				}
				partition_of_intervals[colour[i]-1]=J;
			}
		}
	}

	// retain the groups
	// std::vector<std::vector<size_t> > group_indices(size_t(max_c));
	std::vector<std::vector<size_t> > group_indices(max_c);
	for(size_t i=0;i<graph_size;i++)
	{
		assert(colour[i]>0); // no uncolored nodes
		size_t c = colour[i]-1;
		std::vector<size_t>& group=group_indices[c];
		group.push_back(i);
	}

	return std::make_pair(partition_of_intervals,group_indices);
} 

/**
 * Return the set of common subintervals of a vectors of intervals and an interval.
 *
 * The set is empty if and only if an empty vector is returned.
 * @author Goran Frehse
 */

inline std::vector<interval> intersect_intervals(const std::vector<interval>& itv_vec1,
		const interval& itv2) {
	std::vector<interval> res_vec;
	for (std::vector<interval>::const_iterator it1 = itv_vec1.begin();
			it1 != itv_vec1.end(); ++it1) {
		interval res_itv = intersect(*it1, itv2);
		if (!empty(res_itv)) {
			res_vec.push_back(res_itv);
		}
	}
	return res_vec;
}

/**
 * Return the set of common subintervals of two vectors of intervals
 *
 * The set is empty if and only if an empty vector is returned.
 * @author Goran Frehse
 */

inline std::vector<interval> intersect_intervals(const std::vector<interval>& itv_vec1,
		const std::vector<interval>& itv_vec2) {
	std::vector<interval> res_vec;
	for (std::vector<interval>::const_iterator it1=itv_vec1.begin();it1!=itv_vec1.end();++it1) {
		for (std::vector<interval>::const_iterator it2=itv_vec2.begin();it2!=itv_vec2.end();++it2) {
			interval res_itv=intersect(*it1,*it2);
			if (!empty(res_itv)) {
				res_vec.push_back(res_itv);
			}
		}
	}
	return res_vec;
}

/** Returns the convex hull of a vector of intervals
 *
 * @attention Assumes that the vector is ordered in the following
 * sense: The smallest value is in the first element
 * of the vector and the largest value is in the last element.
 *
 * @author Goran Frehse
 */
inline interval convex_hull(const std::vector<interval>& intv_vec) {
	if (intv_vec.size()<1) {
		return interval::empty();
	} else if (intv_vec.size()==1) {
		return *intv_vec.begin();
	} else {
		precision_type t_beg = intv_vec.begin()->lower();
		precision_type t_end = intv_vec.rbegin()->upper();
		return interval(t_beg,t_end);
	}
}


/**
 * Returns the tightened set of coarse subintervals
 *
 * The intended application is that the intervals in vec_coarse has been computed using coarse
 * precision, and vec_tight using finer precision. If vec_tight has more intervals than vec_coarse,
 * it may be preferable to use vec_coarse, but tightened to exclude parts that don't have a
 * correspondance in vec_tight.
 *
 * It is assumed that every interval in vec_tight is fully contained in one of the intervals
 * in vec_coarse.
 *
 * The set is empty if and only if an empty vector is returned.
 * @author Goran Frehse
 */
inline std::vector<interval> tighten_subintervals(
		const std::vector<interval>& vec_coarse,
		const std::vector<interval>& vec_tight) {
	typedef std::vector<interval> interval_vec;

	if (vec_tight.size() > 1) {
		interval_vec res_vec;
		for (interval_vec::const_iterator it = vec_coarse.begin();
				it != vec_coarse.end(); ++it) {
			interval_vec inters_vec = intersect_intervals(vec_tight, *it);
			if (!inters_vec.empty()) {
				interval tightened_intv = convex_hull(inters_vec);
				res_vec.push_back(tightened_intv);
			}
		}
		return res_vec;
	} else {
		return vec_tight;
	}
}

/** Stream interval
 * @author Goran Frehse
 */
inline std::ostream& operator<<(std::ostream& os,
		const interval& itv) {
	os << "[" << itv.lower() << "," << itv.upper() << "]";
	return os;
}

/** Stream interval vector
 * @author Goran Frehse
 */
inline std::ostream& operator<<(std::ostream& os,
		const std::vector<interval>& itvs) {
	for (std::vector<interval>::const_iterator it = itvs.begin(); it
			!= itvs.end(); ++it) {
		os << *it;
		if (it != itvs.end())
			os << ",";
	}
	return os;
}

/** Prints an STL vector in the format:
 * [0:element0,1:element1,...,n-1:element(n-1)]
 *
 *  @author Goran Frehse
 */
template<typename T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	os << "[";
	for (unsigned int i = 0; i < v.size(); ++i) {
		if (i > 0)
			os << ",";
		os << i << ":" << v[i];
	}
	os << "]";
	return os;
}

}


#endif /* INTERVAL_H_ */
