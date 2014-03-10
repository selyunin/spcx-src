/*
 * basic_exception.h
 *
 *  Created on: Oct 14, 2010
 *      Author: frehse
 */

#ifndef BASIC_EXCEPTION_H_
#define BASIC_EXCEPTION_H_

//#include <iostream>
#include <stdexcept>
#include <string>

/** A basic class for exception handling. */
class basic_exception: public std::runtime_error {
public:
	/** Construct an exception with error message */
	explicit basic_exception(const std::string& msg);

	/** Construct an exception with error message and attach another
	 * exception as cause. */
	basic_exception(const std::string& msg, const basic_exception& cause);

	/** Construct an exception with error message and attach another
	 * exception as cause.
	 *
	 * Tries to cast cause to basic_exception first. */
	basic_exception(const std::string& msg, const std::exception& cause);

	/** Convert a std::exception to basic_exception.
	 *
	 * Tries to cast e to basic_exception first. */
	explicit basic_exception(const std::exception& e);

	/** Copy constructor with deep copy. */
	basic_exception(const basic_exception& e);

	/** Virtual destructor */
	virtual ~basic_exception() throw ();

	/** Assignment with deep copy. */
	basic_exception& operator=(const basic_exception& e);

	/** Returns the message of *this.
	 *
	 * Overloads std::runtime_error. */
	virtual const char* what() const throw() ;

private:
	void prepare_msg();
	basic_exception* my_cause;
	std::string my_msg;
};

#endif /* BASIC_EXCEPTION_H_ */
