/*
 * math_exception.h
 *
 *  Created on: Oct 14, 2010
 *      Author: frehse
 */

#ifndef MATH_EXCEPTION_H_
#define MATH_EXCEPTION_H_

#include "utility/basic_exception.h"

namespace math {

/** A basic class for handling math exceptions. */
class math_exception: public basic_exception {
public:
	math_exception(const std::string& msg) :
			basic_exception(msg) {
	}
	;
	/** Construct an exception with error message and attach another
	 * exception as cause. */
	math_exception(const std::string& msg, const basic_exception& cause) :
			basic_exception(msg, cause) {
	}
	;

	/** Construct an exception with error message and attach another
	 * exception as cause.
	 *
	 * Tries to cast cause to basic_exception first. */
	math_exception(const std::string& msg, const std::exception& cause) :
			basic_exception(msg, basic_exception(cause)) {
	}
	;
};

}

#endif /* MATH_EXCEPTION_H_ */
