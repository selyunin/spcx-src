/*
 * basic_exception.cpp
 *
 *  Created on: Oct 17, 2010
 *      Author: frehse
 */

#include "basic_exception.h"

basic_exception::basic_exception(const std::string& msg) :
	std::runtime_error(msg), my_cause(0) {
	prepare_msg();
}

basic_exception::basic_exception(const std::string& msg,
		const basic_exception& cause) :
	std::runtime_error(msg), my_cause(new basic_exception(cause)) {
	prepare_msg();
}

basic_exception::basic_exception(const std::string& msg,
		const std::exception& cause) :
	std::runtime_error(msg) {
	if (const basic_exception* b = dynamic_cast<const basic_exception*>(&cause)) {
		my_cause = new basic_exception(*b);
	} else {
		my_cause = new basic_exception(cause.what());
	}
	prepare_msg();
}

basic_exception::basic_exception(const std::exception& e) :
			std::runtime_error(e.what()) {
	if (const basic_exception* b = dynamic_cast<const basic_exception*>(&e)) {
		*this = *b;
	} else {
		*this = basic_exception(e.what());
	}
}

basic_exception::basic_exception(const basic_exception& e) :
	std::runtime_error((std::runtime_error) e), my_cause(0) {
	if (e.my_cause) {
		my_cause = new basic_exception(*e.my_cause);
	}
	prepare_msg();
}

basic_exception::~basic_exception() throw () {
	delete my_cause;
}

basic_exception& basic_exception::operator=(const basic_exception& e) {
	std::runtime_error::operator=(e);
	if (e.my_cause) {
		my_cause = new basic_exception(*e.my_cause);
	}
	prepare_msg();
	return *this;
}

const char* basic_exception::what() const throw () {
	return my_msg.c_str();
}

void basic_exception::prepare_msg() {
	my_msg = std::runtime_error::what();
	if (my_cause)
		my_msg = my_msg + "\ncaused by: " + std::string(my_cause->what());
}

