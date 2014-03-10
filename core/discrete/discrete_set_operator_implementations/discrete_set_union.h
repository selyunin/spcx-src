/*
 * discrete_set_union.h
 *
 *  Created on: Mar 22, 2011
 *      Author: mgoyal
 */

#ifndef DISCRETE_SET_UNION_H_
#define DISCRETE_SET_UNION_H_

#include "boost/shared_ptr.hpp"
#include <iostream>
#include <stdexcept>

/** Forward declarations */
namespace discrete {
class discrete_set;
typedef boost::shared_ptr<discrete_set> discrete_set_ptr;
typedef boost::shared_ptr<const discrete_set> discrete_set_const_ptr;
class singleton_set;
class discrete_set_stl_set;
}

namespace discrete {

/** Declaration of the union operator and its wrapper. */
template<typename T1, typename T2> class union_operator {
public:
	static discrete_set_ptr implement(const T1* t1, const T2* t2) {
		//throw std::runtime_error("union_operator : missing implementation");
		return discrete_set_ptr();
	}
	;
};

template<> class union_operator<discrete_set_stl_set,
		discrete_set_stl_set> {
public:
	static discrete_set_ptr implement(const discrete_set_stl_set* p1,
			const discrete_set_stl_set* p2);
};

template<> class union_operator<discrete_set_stl_set, singleton_set> {
public:
	static discrete_set_ptr implement(const discrete_set_stl_set* p1,
			const singleton_set* p2);
};

template<> class union_operator<singleton_set, discrete_set_stl_set> {
public:
	static discrete_set_ptr implement(const singleton_set* p1,
			const discrete_set_stl_set* p2);
};

template<> class union_operator<singleton_set, singleton_set> {
public:
	static discrete_set_ptr implement(const singleton_set* p1,
			const singleton_set* p2);
};

}

#endif /* DISCRETE_SET_INTERSECTION_H_ */
