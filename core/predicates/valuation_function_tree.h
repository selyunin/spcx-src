#ifndef VALUATION_FUNCTION_TREE_H_
#define VALUATION_FUNCTION_TREE_H_

//#include <boost/math/special_functions/sqrt1pm1.hpp>
#include "utility/tree_node.h"
#include "math/type_conversion.h"
#include "core/predicates/valuation_function.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/predicates/valuation_function_tree_utility.h"
#include "global/global_types.h"
#include "math/scalar_types/rational.h"

namespace valuation_functions {

/** Evaluates an arithmetic expression given by a valuation_function_tree.
 * The expressions are evaluated as objects of type scalar_type.
 * The type of constants, const_type, was fixed when the tree was constructed, so it is independent
 * of scalar_type. const_type is cast to scalar_type whenever a constant is evaluated.
 */
template <typename scalar_type, typename const_type> class arithmetic_evaluator {
public:
	virtual ~arithmetic_evaluator() {
	}
	;

	virtual scalar_type const_node_eval(const_node<const_type>* p,
			const typename variable_valuation<scalar_type>::const_ptr& v) {
		return convert_element<scalar_type>(p->my_val);
	}
	;

	virtual scalar_type variable_node_eval(variable_node* p,
			const typename variable_valuation<scalar_type>::const_ptr& v) {
		return v->get_value(p->my_id);
	}
	;

	virtual scalar_type arithmetic_node_eval(arithmetic_node* p,
			const typename variable_valuation<scalar_type>::const_ptr& v) {
		scalar_type x1=arithmetic_eval(p->child1, v);
		if (p->my_op==NEG) {
//			std::cout << "I'm negating " << x1 <<", the result is " << (-x1) << std::endl;
			return -x1;
		}
		/*else if(p->my_op==SQRT) {
			return (scalar_type)(-x1);
		}*/
		scalar_type x2=arithmetic_eval(p->child2, v);
		if (p->my_op==ADD) {
			return x1+x2;
		} else if (p->my_op==SUB) {
			return x1-x2;
		} else if (p->my_op==MUL) {
			return x1*x2;
		} else if (p->my_op==DIV) {
			return x1/x2;
		} else if (p->my_op==POW) {
			return pow(x1,x2);
		}
		else {
			throw std::runtime_error("unknown arithmetic operation");
		}
	}
	;


	virtual scalar_type arithmetic_eval(const tree::node::ptr& p,
			const typename variable_valuation<scalar_type>::const_ptr& v) {
		if (arithmetic_node* q = dynamic_cast<arithmetic_node*>(p.get())) {
			return arithmetic_node_eval(q, v);
		} else if (const_node<const_type>* q
				= dynamic_cast<const_node<const_type>*>(p.get())) {
			return const_node_eval(q, v);
		} else if (variable_node* q = dynamic_cast<variable_node*>(p.get())) {
			return variable_node_eval(q, v);
		} else
			throw std::runtime_error("unmatched node type in valuation_function_tree/arithmetic_eval");
	}
	;

};

/** Constructs a tree that evaluates to a boolean type (eval_type), taking as valuation
 * a scalar. There are no boolean variables (unless scalar_type is boolean).
 */

template <typename bool_type, typename scalar_type, typename const_type> class predicate_evaluator :
	public arithmetic_evaluator<scalar_type,const_type> {
public:
	virtual ~predicate_evaluator() {
	}
	;

	virtual bool_type const_bool_node_eval(const_node<bool_type>* p,
			const typename variable_valuation<scalar_type>::const_ptr& v) {
		return (bool_type)(p->my_val);
	}
	;

	virtual bool_type boolean_node_eval(boolean_node* p,
			const typename variable_valuation<scalar_type>::const_ptr& v) {
		bool_type x1=boolean_eval(p->child1, v);
		if (p->my_op==NOT) {
			return !x1;
		}
		bool_type x2=boolean_eval(p->child2, v);
		if (p->my_op==AND) {
			return x1 && x2;
		} else if (p->my_op==OR) {
			return x1 || x2;
		} else {
			throw std::runtime_error("unknown boolean operation");
		}
	}
	;

	virtual bool_type comparison_node_eval(comparison_node* p,
			const typename variable_valuation<scalar_type>::const_ptr& v) {
		scalar_type x1=arithmetic_eval(p->child1, v);
		scalar_type x2=arithmetic_eval(p->child2, v);
		if (p->my_op==LT) {
			return x1 < x2;
		} else if (p->my_op==LE) {
			return x1 <= x2;
		} else if (p->my_op==GT) {
			return x1 > x2;
		} else if (p->my_op==GE) {
			return x1 >= x2;
		} else if (p->my_op==EQ) {
			return x1 == x2;
		} else {
			throw std::runtime_error("unknown comparison operation");
		}
	}
	;


	virtual bool_type boolean_eval(const tree::node::ptr& p,
			const typename variable_valuation<scalar_type>::const_ptr& v) {
		if (boolean_node* q = dynamic_cast<boolean_node*>(p.get())) {
			return boolean_node_eval(q, v);
		} else if (comparison_node* q = dynamic_cast<comparison_node*>(p.get())) {
			return comparison_node_eval(q, v);
		} else if (const_node<bool_type>* q
				= dynamic_cast<const_node<bool_type>*>(p.get())) {
			return const_bool_node_eval(q, v);
		}
		/* we can' deal with bool variables at this time
		 else if (variable_node* q = dynamic_cast<variable_node*>(p.get())) {
		 return variable_node_eval(q, v);
		 } */
		else
			throw std::runtime_error("unmatched node type");
	}
	;

};

template <typename bool_type, typename scalar_type, typename const_type> class assignment_evaluator :
	public predicate_evaluator<bool_type, scalar_type,const_type> {
public:
	virtual ~assignment_evaluator() {
	}
	;

	virtual bool_type boolean_node_eval(boolean_node* p,
			const typename variable_valuation<scalar_type>::const_ptr& v) {
		bool_type x1=assignment_eval(p->child1, v);
		if (p->my_op==NOT) {
			return !x1;
		}
		bool_type x2=assignment_eval(p->child2, v);
		if (p->my_op==AND) {
			return x1 && x2;
		} else if (p->my_op==OR) {
			return x1 || x2;
		} else {
			throw std::runtime_error("unknown boolean operation");
		}
	}
	;

	virtual bool_type assignment_eval(const tree::node::ptr& p,
			const typename variable_valuation<scalar_type>::const_ptr& v) {
		if (boolean_node* q = dynamic_cast<boolean_node*>(p.get())) {
			return boolean_node_eval(q, v);
		} else if (comparison_node* q = dynamic_cast<comparison_node*>(p.get())) {
			return comparison_node_eval(q, v);
		} else if (const_node<bool_type>* q
				= dynamic_cast<const_node<bool_type>*>(p.get())) {
			return const_bool_node_eval(q, v);
		} else if (p == tree::node::null_node()) {
			return (bool_type)(true);
		}
		/* we can' deal with bool variables at this time
		 else if (variable_node* q = dynamic_cast<variable_node*>(p.get())) {
		 return variable_node_eval(q, v);
		 } */
		else
			return bool_type(arithmetic_eval(p, v));
			//throw std::runtime_error("unmatched node type");
	}
	;
};

/** A valuation_function wrapper for a predicate tree.
 * A predicate is a \p valuation_function from some scalar_type to some bool_type.
 * The constants of the tree are of type const_type.
 */
template<typename bool_type, typename scalar_type, typename const_type> class predicate_tree :
	public valuation_function<bool_type,scalar_type> {
public:
	/** We want to use \p ptr from shared_ptr_user to refer to the base class,
	 * so here we define new pointers for the tree. */
	typedef boost::shared_ptr<predicate_tree<bool_type,scalar_type,const_type> >
			tree_ptr;
	typedef boost::shared_ptr<const predicate_tree<bool_type,scalar_type,const_type> >
			tree_const_ptr;
	typedef predicate_evaluator<bool_type,scalar_type,const_type> evaluator_type;
	typedef boost::shared_ptr<evaluator_type> evaluator_ptr;

	predicate_tree(evaluator_ptr eval, const tree::node::ptr& root) :
		my_root(root), my_eval(eval) {
	}
	;
	predicate_tree(const tree::node::ptr& root) :
		my_root(root) {
		evaluator_ptr p_eval= evaluator_ptr(new evaluator_type);
		my_eval=p_eval;
	}
	;
	virtual ~predicate_tree() {
	}
	;
	virtual const tree::node::ptr& get_root() {
		return my_root;
	}
	;
	virtual void set_root(const tree::node::ptr& root) {
		my_root=root;
	}
	;
	bool_type eval(const typename variable_valuation<scalar_type>::const_ptr& v) const {
		return my_eval->boolean_eval(my_root, v);
	}
	;
//	bool_type print(const typename variable_valuation<scalar_type>::const_ptr& v,
//			std::ostream& os) const {
//		return my_print->boolean_printer(my_root, v, os);
//	}
//	;
	variable_id_set get_variable_ids() const {
		return valuation_functions::get_variable_ids(my_root);
	}
	;
private:
	tree::node::ptr my_root;
	evaluator_ptr my_eval;
};


/** A valuation_function wrapper for a predicate and assignment tree.
 */
template<typename bool_type, typename scalar_type, typename const_type> class assignment_tree :
	public valuation_function<bool_type,scalar_type> {
public:
	/** We want to use \p ptr from shared_ptr_user to refer to the base class,
	 * so here we define new pointers for the tree. */
	typedef boost::shared_ptr<assignment_tree<bool_type,scalar_type,const_type> >
			tree_ptr;
	typedef boost::shared_ptr<const assignment_tree<bool_type,scalar_type,const_type> >
			tree_const_ptr;
	typedef assignment_evaluator<bool_type,scalar_type,const_type> evaluator_type;
	typedef boost::shared_ptr<evaluator_type> evaluator_ptr;

	assignment_tree(evaluator_ptr eval, const tree::node::ptr& root) :
		my_root(root), my_eval(eval) {
	}
	;
		assignment_tree(const tree::node::ptr& root) :
		my_root(root) {
		evaluator_ptr p_eval= evaluator_ptr(new evaluator_type);
		my_eval=p_eval;
	}
	;
	virtual ~assignment_tree() {
	}
	;
	virtual const tree::node::ptr& get_root() {
		return my_root;
	}
	;
	virtual void set_root(const tree::node::ptr& root) {
		my_root=root;
	}
	;
	bool_type eval(const typename variable_valuation<scalar_type>::const_ptr& v) const {
		return my_eval->assignment_eval(my_root, v);
	}
	;
//	bool_type print(const typename variable_valuation<scalar_type>::const_ptr& v,
//			std::ostream& os) const {
//		return my_print->assignment_printer(my_root, v, os);
//	}
//	;
	variable_id_set get_variable_ids() const {
		return get_variable_ids(my_root);
	}
	;
private:
	tree::node::ptr my_root;
	evaluator_ptr my_eval;
//	print_ptr my_print;
};

/** A valuation_function wrapper for an arithmetic tree.
 * All evaluation is carried out on objects of type scalar_type.
 * The constants of the tree are of type const_type.
 */
template<typename scalar_type, typename const_type> class arithmetic_tree :
	public valuation_function<scalar_type,scalar_type> {
public:
	/** We want to use \p ptr from shared_ptr_user to refer to the base class,
	 * so here we define new pointers for the tree. */
	typedef boost::shared_ptr<arithmetic_tree<scalar_type,const_type> > tree_ptr;
	typedef boost::shared_ptr<const arithmetic_tree<scalar_type,const_type> >
			tree_const_ptr;
	typedef arithmetic_evaluator<scalar_type,const_type> evaluator_type;
	typedef boost::shared_ptr<evaluator_type> evaluator_ptr;

	arithmetic_tree(evaluator_ptr eval, const tree::node::ptr& root) :
		my_root(root), my_eval(eval) {
	}
	;
	arithmetic_tree(const tree::node::ptr& root) :
		my_root(root) {
		evaluator_ptr p_eval= evaluator_ptr(new evaluator_type);
		my_eval=p_eval;
	}
	;
	virtual ~arithmetic_tree() {
	}
	;
	virtual const tree::node::ptr& get_root() {
		return my_root;
	}
	;
	virtual void set_root(const tree::node::ptr& root) {
		my_root=root;
	}
	;
	scalar_type eval(const typename variable_valuation<scalar_type>::const_ptr& v) const {
		return my_eval->arithmetic_eval(my_root, v);
	}
	;
//	scalar_type print(const typename variable_valuation<scalar_type>::const_ptr& v,
//			std::ostream& os) const {
//		return my_print->arithmetic_printer(my_root, v, os);
//	}
//	;
	variable_id_set get_variable_ids() const {
		return valuation_functions::get_variable_ids(my_root);
	}
	;

private:
	tree::node::ptr my_root;
	evaluator_ptr my_eval;
};

}

#endif /*VALUATION_FUNCTION_TREE_H_*/
