/*
 * symbolic_state_acceptor.h
 *
 *  Created on: Sep 25, 2009
 *      Author: frehse
 */

#ifndef SYMBOLIC_STATE_ACCEPTOR_H_
#define SYMBOLIC_STATE_ACCEPTOR_H_

#include "boost/shared_ptr.hpp"

/** Forward declarations. */
namespace hybrid_automata {
class symbolic_state;
typedef boost::shared_ptr<symbolic_state> symbolic_state_ptr;
class adapt_discrete_set_visitor;
class adapt_continuous_set_visitor;
}

namespace hybrid_automata {

/** An interface class that accepts a symbolic state.
 *
 * It serves as a common interface to multiple objects
 * that process symbolic states, such as symbolic_state_collection,
 * PWL, etc. */

class symbolic_state_acceptor {
public:
	virtual ~symbolic_state_acceptor() {
	}
	;

	/** Process the symbolic state sstate. */
	virtual void accept(const symbolic_state_ptr& sstate) = 0;
};

}

#endif /* SYMBOLIC_STATE_ACCEPTOR_H_ */
