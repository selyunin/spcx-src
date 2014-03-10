/*
 * continuous_dynamics_base.h
 *
 *  Created on: Sep 10, 2009
 *      Author: frehse
 */

#ifndef CONTINUOUS_DYNAMICS_BASE_H_
#define CONTINUOUS_DYNAMICS_BASE_H_

#include "boost/shared_ptr.hpp"
#include <boost/enable_shared_from_this.hpp>
#include "utility/printable.h"
#include "math/vdom/variable.h"
#include "math/vdom/vdom_vector.h"
#include "utility/dispatching/double_dispatch.h"

//#include "../../../utility/index_to_variable_id_map.h"
#include "utility/tree_node.h"
//#include "../../../math/matrix.h"
//#include "../../../math/lin_expression.h"
#include "core/continuous/continuous_set.h"

#include "core/continuous/continuous_dynamics/continuous_dynamics_declarations.h"

/** Forward declaration of classes used in header file. */
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;
}

namespace continuous {

typedef tree::node dynamics_predicate;

/** Abstract interface class for representing set-based continuous dynamics, i.e.,
 * set based ODEs of the form \f$ \dot x \in F(x) \f$.
 * \todo Not quite sure yet what the difference to continuous_set_transforms is.
 * Mainly because we want to use these in locations, and not all the rest of them.
 * Well, they're purely syntactic.
 * Maybe they should be derived from continuous_set_transforms?
 *
 * The get_relation function is different from continuous_set_transforms!
 */
class continuous_dynamics : public boost::enable_shared_from_this<
continuous_dynamics>, public virtual printable {
public:
	typedef boost::shared_ptr<continuous_dynamics> ptr;
	typedef boost::shared_ptr<const continuous_dynamics> const_ptr;

	typedef dispatching::dispatcher<continuous_dynamics_typelist> const_visitor;

	/** Return a shared_ptr to *this. */
	ptr get_ptr() {
		return boost::enable_shared_from_this<continuous_dynamics>::shared_from_this();
	}
	;

	/** Return a shared_ptr to const *this. */
	const_ptr get_const_ptr() const {
		return boost::enable_shared_from_this<continuous_dynamics>::shared_from_this();
	}
	;

	virtual ~continuous_dynamics() {
	}
	;

	/** Returns the dynamics in predicate form.*/
	virtual dynamics_predicate::ptr get_predicate() const = 0;

	/** Returns the ids of all variables in the domain or codomain.
	 *  A variable x and its derivative \dot x are considered as being the same variable. */
	virtual variable_id_set get_variable_ids() const = 0;

	/** Returns the ids of variables whose derivative is unconstrained.
	 *
	 *  The set can be conservative in the sense that unconstrained variables
	 *  may erroneously be considered constrained. I.e., if a variable
	 *  is returned by this function, one can be absolutely sure that it
	 *  is unconstrained.
	 *
	 *  Recall that the flow F \subset R^\dot X \times R^X.
	 *  Let v, \dot v, \dot w be valuations over X, respectivbely \dot X.
	 *  A variable x is considered unconstrained if
	 *  for all (\dot v,v) in F, \dot w in R^\dot X with
	 *  \dot v(y)=\dot w(y) if y \neq x holds
	 *  (\dot w,v) in F.
	 *  */
	virtual variable_id_set get_unconstrained_variable_ids() const = 0;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const = 0;

	/** Accept a const_visitor. A const_visitor must provide the function
	 * <code> void dispatch(const T* c) </code>
	 * for all derived classes of continuous_dynamics that are listed in
	 * continuous_dynamics_typelist (see continuous_dynamics_declarations),
	 * plus a "dummy" function for the type void (which is never called but
	 * makes things a little easier at the back end).
	 */
	virtual void accept(const_visitor& d) const = 0;
};

}

#endif /* CONTINUOUS_DYNAMICS_BASE_H_ */
