#ifndef DISPATCHER_H_
#define DISPATCHER_H_

/** Generic visitor based on typelists. Can be used for double dispatching.
 */

#include "utility/typelist.h"
#include <stdexcept>

namespace dispatching {

/** Declaration of (non-const) dispatcher. */
template<typename types_to_visit> class dispatcher :
	public dispatcher<typename types_to_visit::tail> {
public:
	virtual ~dispatcher() {
	}
	;
	/** The using declaration is necessary for resolving the dispatch calls (why?). */
	using dispatcher<typename types_to_visit::tail>::dispatch;
	virtual void dispatch(const typename types_to_visit::head* c) = 0;
};

/** Dispatcher that terminates the typelist.*/
template<> class dispatcher<null_typelist> {
public:
	virtual ~dispatcher() {
	}
	;
	/** The dispatch method is necessary, otherwise the using declaration
	 * in the above template fails. */
	virtual void dispatch(const void* c) {
		throw std::runtime_error("void dispatcher call");
	}
	;
};

}

#endif /*DISPATCHER_H_*/
