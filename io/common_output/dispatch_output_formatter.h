/*
 * dispatch_output_formatter.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef DISPATCH_OUTPUT_FORMATTER_H_
#define DISPATCH_OUTPUT_FORMATTER_H_

#include "output_formatter.h"
#include "utility/dispatching/single_dispatch.h"
#include "core/continuous/continuous_set.h"

namespace io {

template<template<typename > class implementor>
class dispatch_output_formatter_wrapper {
public:
	dispatch_output_formatter_wrapper(output_formatter& caller) :
		my_caller(caller) {
	}
	;
	template<typename T> void implement(const T* c) {
		implementor<T>::output(my_caller, *c);
	}
	;
private:
	output_formatter& my_caller;
};

/** This class simplifies (a little) the interface to a single dispatch
 * for downcasting continuous sets.
 *
 * A deriving class must provide a template class implementor, with
 * a member function void output(std::ostream&,T).
 */
template<template<typename > class implementor,
		template<typename > class caster = dispatching::no_caster>
class dispatch_output_formatter: public output_formatter {
public:
	typedef output_formatter base_class;

	dispatch_output_formatter(std::ostream& os) :
		base_class(os) {
	}
	;
	virtual ~dispatch_output_formatter() {
	}
	;
	virtual void output(const continuous::continuous_set& c) {
		typedef dispatch_output_formatter_wrapper<implementor>
				dispatch_implementor;
		dispatch_implementor f(*this);
		dispatching::single_dispatch<dispatch_implementor,
				continuous::continuous_set_typelist,
				continuous::continuous_set, caster>(f, &c);
	}
	;
	virtual void output(const hybrid_automata::symbolic_state& sstate) {
		base_class::output(sstate);
	}
	;
	virtual void output(const hybrid_automata::symbolic_state_collection& sstates) {
		base_class::output(sstates);
	}
	;
};

}

#endif /* DISPATCH_OUTPUT_FORMATTER_H_ */
