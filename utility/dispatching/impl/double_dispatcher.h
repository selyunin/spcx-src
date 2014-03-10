#ifndef DOUBLE_DISPATCHER_H_
#define DOUBLE_DISPATCHER_H_

/** Recursive double dispatching, visitor-based implementation of multimethods.
 * It converts two base pointers to the corresponding derived pointers using
 * a visitor (dispatcher).
 * See Shopyrin: Multimethods in C++ Using Recursive Deferred Dispatching
 * */

#include "utility/dispatching/dispatcher.h"
#include "utility/value_owner.h"

namespace dispatching {

template<typename return_type, class implementor, class derived_class1,
		typename types_to_visit, typename types_to_instantiate> class double_dispatcher :
	public double_dispatcher<return_type, implementor, derived_class1, types_to_visit, typename types_to_instantiate::tail> {
public:
	double_dispatcher(const derived_class1* derived1,
			value_owner<return_type>* rec) :
		double_dispatcher< return_type, implementor, derived_class1,
				types_to_visit, typename types_to_instantiate::tail>(derived1,
				rec) {
	}
	;
	virtual ~double_dispatcher() {
	}
	;
	virtual void dispatch(const typename types_to_instantiate::head* derived2) {
		this->my_recipient->set(implementor::implement(this->my_derived1,
				derived2));
	}
	;
};

template<typename return_type, class implementor, class derived_class1,
		typename types_to_visit> class double_dispatcher< return_type,
		implementor, derived_class1, types_to_visit, null_typelist> :
	public dispatcher<types_to_visit> {
public:
	double_dispatcher(const derived_class1* derived1,
			value_owner<return_type>* rec) :
		my_derived1(derived1), my_recipient(rec) {
	}
	;
	virtual ~double_dispatcher() {
	}
	;
	virtual void dispatch(const void* c) {
		throw std::runtime_error("void double_dispatcher call");
	}
	;

protected:
	const derived_class1* my_derived1;
	value_owner<return_type>* my_recipient;
};

}

#endif /*DOUBLE_DISPATCHER_H_*/
