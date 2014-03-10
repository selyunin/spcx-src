#ifndef PRINTABLE_H_
#define PRINTABLE_H_

#include <iostream>
//#include "boost/shared_ptr.hpp"

/** Provides an overloadable version of operator<< for easy
 * output of derived classes. The class printable provides the
 * abstract interface print(ostream).
 * */
class printable {
public:
	virtual ~printable() {
	}
	;

	/**
	 * Output as a stream of characters.
	 */
	virtual void print(std::ostream& os) const = 0;

};

/**
 * Output as a stream of characters. Calls the print() method of the class.
 */
std::ostream& operator<<(std::ostream& os, const printable& p);
/**
 * Output as a stream of characters. Calls the print() method of the class.
 */
std::ostream& operator<<(std::ostream& os, const printable*& p);

//std::ostream& operator<<(std::ostream& os, const boost::shared_ptr<printable> p);
//std::ostream& operator<<(std::ostream& os, const boost::shared_ptr<const printable> p);

#endif /*PRINTABLE_H_*/
