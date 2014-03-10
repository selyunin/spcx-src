/*
 * piecewise_linear_interval_function.h
 *
 *  Created on: Oct 18, 2012
 *      Author: kateja
 */

#ifndef PIECEWISE_LINEAR_INTERVAL_FUNCTION_H_
#define PIECEWISE_LINEAR_INTERVAL_FUNCTION_H_

#include "piecewise_linear_function_operators.h"

namespace plif {

/**
 * Piecewise Linear Interval Function class
 *
 * Has two data fields, lower and upper to identify the lower and upper bound PLF's of this PLIF
 */
class piecewise_linear_interval_function
{
private:
	piecewise_linear_function lower, upper;		/** lower is lower bounding PLF, upper is upper bounding PLF */
	interval domain;
	bool left_closed;
	bool right_closed;
	bool left_unbounded;
	bool right_unbounded;
public:
	/**
	 *Default Constructor
	 */
	piecewise_linear_interval_function()
	{}
	/**
	 *Constructor
	 */
	piecewise_linear_interval_function(piecewise_linear_function low, piecewise_linear_function up)
	{
		set_lower(low);
		set_upper(up);
		//domain = low.get_domain();
		//if (low.get_domain().lower() != up.get_domain().lower() || low.get_domain().upper() != up.get_domain().upper() ) {
		if (!is_MEQ(low.get_domain().lower(), up.get_domain().lower())
				|| !is_MEQ(low.get_domain().upper(), up.get_domain().upper())) {
			std::cerr << "lower: " << lower << ", upper: " << upper
					<< ", diff lower: "
					<< low.get_domain().lower() - up.get_domain().lower()
					<< ", diff upper: "
					<< low.get_domain().upper() - up.get_domain().upper()
					<< std::endl;
			throw std::runtime_error(
					"initializing PLIF with different domains");
		}
		domain = interval(intersect(low.get_domain(), up.get_domain()));
	}
	/** Get the number of breakpoints (max over lower and upper)
	 */
	size_t size() {
		return std::max(lower.size(),upper.size());
	}
	/**
	 * Set lower bound PLF
	 */
	void set_lower(piecewise_linear_function f) {
		lower = f;
		lower.set_rounding_down();
	}
	/**
	 * Get lower bound PLF
	 */
	const piecewise_linear_function& get_lower() const {
		return lower;
	}
	/**
	 * Set upper bound PLF
	 */
	void set_upper(piecewise_linear_function f) {
		upper = f;
		upper.set_rounding_up();
	}
	/**
	 * Get upper bound PLF
	 */
	const piecewise_linear_function& get_upper() const {
		return upper;
	}
	void set_domain(const interval& A)
		{domain = A;}
	const interval& get_domain() const
		{return domain;}
	void set_left_closed(bool b)
		{left_closed = b;}
	const bool& is_left_closed() const
		{return left_closed;}
	void set_right_closed(bool b)
		{right_closed = b;}
	const bool& is_right_closed() const
		{return right_closed;}
	void set_left_unbounded(bool b)
		{left_unbounded = b;}
	const bool& is_left_unbounded() const
		{return left_unbounded;}
	void set_right_unbounded(bool b)
		{right_unbounded = b;}
	const bool& is_right_unbounded() const
		{return right_unbounded;}
	/**
	 *Function to display the PLIF
	 */
	 void display(std::ostream& os = std::cout) const
	 {
	 	os << "Domain: [" << domain.lower() << ", " << domain.upper() << "]" << std::endl;
	 	os << "Left Closed: " << left_closed << std::endl;
	 	os << "Right Closed: " << right_closed << std::endl;
	 	os << "Lower Function" << std::endl;
	 	lower.display();
	 	os << "Upper Function" << std::endl;
	 	upper.display();
	 }
};

/** Stream to output */
inline
std::ostream& operator<<(std::ostream& os, const piecewise_linear_interval_function& h) {
	h.display(os);
	return os;
}

/** Stream to output */
inline
std::ostream& operator<<(std::ostream& os, const std::vector<piecewise_linear_interval_function>& h_vec) {
	size_t i=0;
	for (std::vector<piecewise_linear_interval_function>::const_iterator it = h_vec.begin(); it!= h_vec.end(); ++it) {
		os << i<<": " << *it << std::endl;
		++i;
	}
	return os;
}

}

#endif /* PIECEWISE_LINEAR_INTERVAL_FUNCTION_H_ */
