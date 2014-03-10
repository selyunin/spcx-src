/*
 * shared_ptr_output.h
 *
 *  Created on: Sep 1, 2009
 *      Author: frehse
 */

#ifndef SHARED_PTR_OUTPUT_H_
#define SHARED_PTR_OUTPUT_H_

#include <boost/shared_ptr.hpp>

/** Output the object pointed to.
 */
template<typename T> std::ostream& operator<<(std::ostream& os,
		const boost::shared_ptr<T>& p) {
	assert(p);
	os << (*p);
	return os;
}

#endif /* SHARED_PTR_OUTPUT_H_ */
