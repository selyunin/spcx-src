/*
 * caster.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef CASTER_H_
#define CASTER_H_

#include <boost/type_traits/is_convertible.hpp>
#include <boost/mpl/if.hpp>

namespace dispatching {

/** Identity caster
 *
 * This caster does not modify the type. */
template<typename T>
class no_caster {
public:
	typedef T result;
};

template<typename T1, typename T2>
class no_caster_double {
public:
	typedef T1 result1;
	typedef T2 result2;
};

/** Caster that casts to target_class if possible */
template<typename T, class target_class>
class convertible_caster {
public:
	static const bool is_target =
			boost::is_convertible<T*, target_class*>::value;
	typedef typename boost::mpl::if_c<is_target, target_class, T>::type result;
};

/** Concatenation of two casters, head and tail.
 *
 * The tail is applied first. */
template<class head, class tail>
class caster_list {
public:
	static const bool is_target = tail::is_target || head::is_target;
	typedef typename boost::mpl::if_c<tail::is_target, typename tail::result, typename head::result>::type result;
};

}

#endif /* CASTER_H_ */
