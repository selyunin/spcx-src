/*
 * invalid_number_exception.h
 *
 *  Created on: Oct 30, 2012
 *      Author: notroot
 */

#ifndef INVALID_NUMBER_EXCEPTION_H_
#define INVALID_NUMBER_EXCEPTION_H_

#include "utility/basic_exception.h"

#include <boost/math/special_functions/fpclassify.hpp>

namespace math {

/** An exception due to a nan or inf */
class invalid_number_exception: public basic_exception {
public:
	invalid_number_exception(const std::string& msg) :
		basic_exception(msg) {
	}
	;
	/** Construct an exception with error message and attach another
	 * exception as cause. */
	invalid_number_exception(const std::string& msg, const basic_exception& cause) :
		basic_exception(msg, cause) {
	}
	;
};

/** A function to check whether a number is finite (neither inf nor nan) */
template<class scalar_type>
bool is_finite(const scalar_type& a) {
	return boost::math::isfinite(a);
}
;

}

#endif /* INVALID_NUMBER_EXCEPTION_H_ */
