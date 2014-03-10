/*
 * printable.cpp
 *
 *  Created on: Sep 1, 2009
 *      Author: frehse
 */

#include "utility/printable.h"

#include <cassert>

std::ostream& operator<<(std::ostream& os, const printable& p) {
	p.print(os);
	return os;
}

std::ostream& operator<<(std::ostream& os, const printable*& p) {
	assert(p);
	p->print(os);
	return os;
}

//std::ostream& operator<<(std::ostream& os, const boost::shared_ptr<printable> p) {
//	assert(p);
//	p->print(os);
//	return os;
//}
//
//std::ostream& operator<<(std::ostream& os, const boost::shared_ptr<const printable> p) {
//	assert(p);
//	p->print(os);
//	return os;
//}
