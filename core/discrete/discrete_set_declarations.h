/*
 * discrete_set_declarations.h
 *
 *  Created on: Sep 8, 2009
 *      Author: frehse
 */

#ifndef DISCRETE_SET_DECLARATIONS_H_
#define DISCRETE_SET_DECLARATIONS_H_

#include "utility/typelist.h"

/** Forward declaration of all classes that derive from continuous_set and
 * that wish to accept a dispatcher. */

namespace discrete {
class singleton_set;
class discrete_set_stl_set;
}

namespace discrete {
typedef typelist<discrete_set_stl_set,
        typelist<singleton_set,
null_typelist>
> discrete_set_typelist;
}

#endif /* DISCRETE_SET_DECLARATIONS_H_ */
