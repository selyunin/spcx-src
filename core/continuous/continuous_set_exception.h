/*
 * continuous_set_exception.h
 *
 *  Created on: Oct 15, 2010
 *      Author: frehse
 */

#ifndef CONTINUOUS_SET_EXCEPTION_H_
#define CONTINUOUS_SET_EXCEPTION_H_

#include "utility/basic_exception.h"

namespace continuous {

/** An exception due to a manipulation of a continuous set. */
class continuous_set_exception: public basic_exception {
public:
	continuous_set_exception(const std::string& msg) :
		basic_exception(msg) {
	}
	;
	/** Construct an exception with error message and attach another
	 * exception as cause. */
	continuous_set_exception(const std::string& msg, const basic_exception& cause) :
		basic_exception(msg, cause) {
	}
	;
};

/** An exception due to a manipulation of an unbounded continuous set. */
class unbounded_set_exception: public continuous_set_exception {
public:
	unbounded_set_exception(const std::string& msg) :
		continuous_set_exception(msg) {
	}
	;
	/** Construct an exception with error message and attach another
	 * exception as cause. */
	unbounded_set_exception(const std::string& msg, const continuous_set_exception& cause) :
		continuous_set_exception(msg, cause) {
	}
	;
};

/** An exception due to a manipulation of an empty continuous set. */
class empty_set_exception: public continuous_set_exception {
public:
	empty_set_exception(const std::string& msg) :
		continuous_set_exception(msg) {
	}
	;
	/** Construct an exception with error message and attach another
	 * exception as cause. */
	empty_set_exception(const std::string& msg, const continuous_set_exception& cause) :
		continuous_set_exception(msg, cause) {
	}
	;
};

}

#endif /* CONTINUOUS_SET_EXCEPTION_H_ */
