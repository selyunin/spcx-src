#ifndef CONTINUOUS_SET_TRANSFORM_H_
#define CONTINUOUS_SET_TRANSFORM_H_

#include <list>
#include <boost/any.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "utility/printable.h"
#include "math/vdom/variable.h"
#include "utility/shared_ptr_user.h"
#include "math/vdom/index_to_variable_id_map.h"
#include "utility/tree_node.h"
#include "core/predicates/valuation_function_tree_utility.h"
#include "math/matrix.h"
//#include "../continuous_set.h"
//#include "../continuous_set_operators.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transform_declarations.h"
#include "utility/dispatching/double_dispatch.h"

namespace continuous {

class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;

typedef continuous_set_ptr relation_ptr;
typedef continuous_set_const_ptr relation_const_ptr;

/** Abstract interface class for representing transformations on continuous sets,
 * i.e., of set-valued functions \f$F:R^n \rightarrow 2^{R^n}\f$. If
 * \f$x=(x_1,...,x_n)\f$ are the (values of the ) variables before the transformation and
 * \f$x'=(x_1',...,x_n')\f$ are the variables after, then \f$x' \in F(x)\f$.
 *
 * Optimizations are possible by implementation in the derived continuous set. */
class continuous_set_transform: public boost::enable_shared_from_this<
		continuous_set_transform>, public virtual printable {
public:
	typedef boost::shared_ptr<continuous_set_transform> ptr;
	typedef boost::shared_ptr<const continuous_set_transform> const_ptr;

	typedef dispatching::dispatcher<continuous_set_transform_typelist> const_visitor;

	/** Return a shared_ptr to *this. */
	ptr get_ptr() {
		return boost::enable_shared_from_this<continuous_set_transform>::shared_from_this();
	}
	;

	/** Return a shared_ptr to const *this. */
	const_ptr get_const_ptr() const {
		return boost::enable_shared_from_this<continuous_set_transform>::shared_from_this();
	}
	;

	virtual ~continuous_set_transform() {
	}
	;

	/** Get the relation of the transformation, i.e., a continuous set describing
	 * the set of pairs \f$(x,x')\f$ for which \f$x'\f$ is one of the possible values
	 * that \f$x\f$ can be transformed into.
	 * cset is used to instantiate a universe or empty relation if necessary.
	 * */
	virtual relation_const_ptr get_relation(continuous_set_const_ptr cset) const = 0;

	/** Get the variables used by the transformation and the variables modified by the
	 * transformation. Assume that \f$x=(x_1,...,x_n)\f$. A variable x_i is used unless
	 * the transform yields the same set of
	 * values regarless of the value of x_i, i.e., x_i is not used if
	 * \f$x'\in F(x_1,...,x_i,...,x_n)  \Rightarrow \forall x_i' : x'\in F(x_1,...,x_i',...,x_n)\f$.
	 * A variable x_i is modified if it can change in the transform, i.e., it is not modified if
	 * \f$\forall x,x' : x' \in F(x)  \Rightarrow x_i'=x_i\f$.
	 */
	virtual void get_used_and_modif_variables(variable_id_set& used_vars,
			variable_id_set& modif_vars) const = 0;

	/**
	 * Output as a stream of characters.
	 */
	virtual void print(std::ostream& os) const = 0;

	/** Accept a const_visitor. A const_visitor must provide the function
	 * <code> void dispatch(const T* c) </code>
	 * for all derived classes of continuous_set_transform that are listed in
	 * continuous_set_transform_typelist (see continuous_set_transform_declarations),
	 * plus a "dummy" function for the type void (which is never called but
	 * makes things a little easier at the back end).
	 */
	virtual void accept(const_visitor& d) const = 0;
};

}

#endif /*CONTINUOUS_SET_TRANSFORM_H_*/
