/*
 * single_dispatch.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef SINGLE_DISPATCH_H_
#define SINGLE_DISPATCH_H_

#include "utility/dispatching/impl/single_dispatcher.h"
#include "utility/dispatching/impl/template_class_implementor_wrapper.h"
#include "utility/dispatching/caster.h"

namespace dispatching {

/** Single dispatch to implementor class.
 *
 * The implementor must provide a function
 * void implement(d) for any derived type d.
 * A useful way to guarantee this is to let implement be a template
 * member function.
 * However, this cannot be partially specialized.
 * For partial specializations, use single_dispatch_tc.
 */
template<class implementor, typename types_to_visit, typename base_class,
		template<typename > class caster> void single_dispatch(
		implementor& impl, const base_class* base) {
	single_dispatcher<implementor, types_to_visit, caster> d(impl);
	base->accept(d);
}
;
//
///** Single dispatch to templated implementor class.
// *
// * The class implementor<T> must be templated, and will be instantiated
// * for all combinations of derived classes.
// * The class must provide
// * void implement(const T* t) for derived types T of base_class.
// * Note that the class can be partially specialized.
// */
//template<template<typename > class template_class_implementor,
//		typename types_to_visit, class base_class,
//		template<typename > class caster = no_caster> void single_dispatch_tc(
//		implementor& impl, const base_class* base) {
//	typedef template_class_implementor_wrapper<template_class_implementor>
//			implementor;
//	single_dispatch<implementor, base_class, types_to_visit> (impl, base);
//}
//;

}

#endif /* SINGLE_DISPATCH_H_ */
