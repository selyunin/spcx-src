/*
 * parse_exception.h
 *
 *  Created on: Jun 9, 2011
 *      Author: notroot
 */

#ifndef PARSE_EXCEPTION_H_
#define PARSE_EXCEPTION_H_

#include <sstream>
#include "utility/basic_exception.h"
//#include <boost/spirit/include/phoenix_function.hpp>

namespace parser {

class parse_exception : public basic_exception {
public:
	/** Construct an exception with error message */
	explicit parse_exception(const std::string& msg);

	/** Construct an exception with error message and attach another
	 * exception as cause. */
	parse_exception(const std::string& msg, const parse_exception& cause);

	/** Construct an exception with error message and attach another
	 * exception as cause.
	 *
	 * Tries to cast cause to basic_exception first. */
	parse_exception(const std::string& msg, const std::exception& cause);
};

struct parse_exception_handler_impl {
	template <class, class, class, class, class>
	   struct result { typedef void type; };

	template <class M, class T1, class T2, class T3, class T4>
	void operator()(const M& msg, const T1& _1, const T2& _2, const T3& _3, const T4& _4 ) const {
		std::stringstream ss;
		ss << "Error parsing " << msg << ": \n" << std::string(_1, _2) << std::endl;
		ss << std::string(_3-_1,' ') << "^" << std::string(_3, _2) << std::endl;
		ss << "Internal grammar state:" << _4;
		throw parse_exception(ss.str());
	}
};

}

#endif /* PARSE_EXCEPTION_H_ */
