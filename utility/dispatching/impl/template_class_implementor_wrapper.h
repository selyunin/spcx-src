/*
 * template_class_implementor_wrapper.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef TEMPLATE_CLASS_IMPLEMENTOR_WRAPPER_H_
#define TEMPLATE_CLASS_IMPLEMENTOR_WRAPPER_H_

#include "utility/dispatching/caster.h"

namespace dispatching {

/** Converts a call to a template member function into a call to a
 * template class.
 */
template<template<typename > class template_class_implementor>
class template_class_implementor_wrapper {
public:
	template<typename T> static bool implement(const T* t) {
		return template_class_implementor<T>::implement(t);
	}
	;
};

/** Converts a call to a template member function into a call to a
 * template class.
 *
 * It optionally takes two casters to cast the derived types.
 */
template<typename return_type,
		template<typename , typename > class template_class_implementor,
		template<typename , typename > class caster = no_caster_double>
class template_class_double_implementor_wrapper {
public:
	template<typename T1, typename T2> static return_type implement(
			const T1* t1, const T2* t2) {
		return template_class_implementor<typename caster<T1, T2>::result1,
				typename caster<T1, T2>::result2>::implement(t1, t2);
	}
	;
};

}

#endif /* TEMPLATE_CLASS_IMPLEMENTOR_WRAPPER_H_ */
