/*
 * support_function_provider_upcaster.h
 *
 *  Created on: Dec 22, 2009
 *      Author: frehse
 */

#ifndef SUPPORT_FUNCTION_PROVIDER_UPCASTER_H_
#define SUPPORT_FUNCTION_PROVIDER_UPCASTER_H_

#include "utility/dispatching/caster.h"

/** Forward declarations */
namespace continuous {
class support_function_provider;
}

namespace continuous {

/** Try to upcast to a support_function_provider.
 *
 * The result is available as the type support_function_provider_upcaster<T>::result.
 */
template<typename T>
class support_function_provider_upcaster : public dispatching::convertible_caster<T,support_function_provider> {
};

}

#endif /* SUPPORT_FUNCTION_PROVIDER_UPCASTER_H_ */
