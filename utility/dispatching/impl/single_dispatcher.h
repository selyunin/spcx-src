/*
 * single_dispatcher.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef SINGLE_DISPATCHER_H_
#define SINGLE_DISPATCHER_H_

/** Recursive single dispatching.
 * It converts a base pointer to the corresponding derived pointer using
 * a visitor (dispatcher).
 * */

#include "utility/dispatching/dispatcher.h"
#include "utility/value_owner.h"
#include "utility/dispatching/caster.h"

namespace dispatching {

/** The single dispatcher calls an implementor class that must provide
 * an member function void implement(derived_class*) for all derived classes
 * (can be a template member function).
 *
 * The implementor class serves to provide a default implementation,
 * and can be fully specialized.
 * For partial specialization or upcasting redirect the implement call
 * from the implementor class to a template class.
 * The template class can be partially specialized etc.
 *
 * Optionally, a (up)cast can be performed by passing a caster, which is a
 * template class defining a type result. The derived class is then cast to
 * caster<derived>::result. This can be used to upcast to classes like
 * polyhedron if different derived classes of polyhedron are available
 * and all of them can be treated the same way.
 */
template<class implementor, typename types_to_visit,
		template<typename > class caster = no_caster,
		typename types_to_instantiate = types_to_visit> class single_dispatcher: public single_dispatcher<
		implementor, types_to_visit, caster,
		typename types_to_instantiate::tail> {
public:
	single_dispatcher(implementor& impl) :
		single_dispatcher<implementor, types_to_visit, caster,
				typename types_to_instantiate::tail> (impl) {
	}
	;
	virtual ~single_dispatcher() {
	}
	;
	virtual void dispatch(const typename types_to_instantiate::head* derived) {
		typedef typename caster<typename types_to_instantiate::head>::result
				casted_type;
		this->my_impl.implement(static_cast<const casted_type*> (derived));
	}
	;
};

template<class implementor, typename types_to_visit,
		template<typename > class caster> class single_dispatcher<implementor,
		types_to_visit, caster, null_typelist> : public dispatcher<
		types_to_visit> {
public:
	single_dispatcher(implementor& impl) :
		my_impl(impl) {
	}
	;
	virtual ~single_dispatcher() {
	}
	;
	virtual void dispatch(const void* c) {
		throw std::runtime_error("void single_dispatcher call");
	}
	;
protected:
	implementor& my_impl;
};

}

#endif /* SINGLE_DISPATCHER_H_ */
