#ifndef CROSS_DISPATCHER_H_
#define CROSS_DISPATCHER_H_

#include "utility/dispatching/impl/double_dispatcher.h"

namespace dispatching {

template<typename return_type, class implementor, class base_class,
		typename types_to_visit, typename types_to_instantiate> class cross_dispatcher :
	public cross_dispatcher<return_type,
		implementor, base_class, types_to_visit, typename types_to_instantiate::tail> {
public:
	cross_dispatcher(const base_class* base2) :
		cross_dispatcher<return_type, implementor, base_class, types_to_visit,
				typename types_to_instantiate::tail>(base2) {
	}
	;
	virtual ~cross_dispatcher() {
	}
	;
	virtual void dispatch(const typename types_to_instantiate::head* derived1) {
		double_dispatcher<return_type, implementor, typename types_to_instantiate::head, types_to_visit,
		types_to_visit> d(derived1, this);
		this->my_base2->accept(d);
	}
	;
};

template<typename return_type, class implementor, class base_class,
		typename types_to_visit> class cross_dispatcher<return_type,
		implementor, base_class, types_to_visit, null_typelist> :
	public dispatcher<types_to_visit>, public value_owner<return_type> {
public:
	cross_dispatcher(const base_class* base2) :
		my_base2(base2) {
	}
	;
	virtual ~cross_dispatcher() {
	}
	;
	virtual void dispatch(const void* derived1) {
		throw std::runtime_error("void cross_dispatcher call. Forgot a class in the typelist?");
	}
	;
protected:
	const base_class* my_base2;
};

}
#endif /*CROSS_DISPATCHER_H_*/
