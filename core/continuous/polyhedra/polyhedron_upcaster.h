/*
 * polyhedron_upcaster.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef POLYHEDRON_UPCASTER_H_
#define POLYHEDRON_UPCASTER_H_

#include <boost/type_traits/is_convertible.hpp>
#include <boost/mpl/if.hpp>

#include "global/global_types.h"

/** Forward declarations */
namespace continuous {
template<typename T> class polyhedron;
}

namespace continuous {

/** Try to upcast to a polyhedron.
 *
 * The result is available as the type polyhedron_upcaster<T>::result.
 *
 * Since boost::is_convertible only
 * applies to full types, we must try all scalar_types with which
 * polyhedron might be instantiated.
 */
template<typename T>
class polyhedron_upcaster {
public:
	typedef continuous::polyhedron<global_types::float_type> poly_float;
	typedef continuous::polyhedron<global_types::rational_type> poly_rational;
	static const bool is_poly_float =
			boost::is_convertible<T*, poly_float*>::value;
	static const bool is_poly_rational = boost::is_convertible<T*,
			poly_rational*>::value;
	static const bool is_target = is_poly_float || is_poly_rational;
	static const bool is_poly = is_poly_float || is_poly_rational;
	typedef typename boost::mpl::if_c<is_poly_float, poly_float,
			typename boost::mpl::if_c<is_poly_rational, poly_rational, T>::type>::type
			result;
};

}

#endif /* POLYHEDRON_UPCASTER_H_ */
