/*
 * valuation_function_tree_ops.h
 *
 *  Created on: Mar 12, 2013
 *      Author: notroot
 */

#ifndef VALUATION_FUNCTION_TREE_OPS_H_
#define VALUATION_FUNCTION_TREE_OPS_H_

#include "valuation_function_tree_nodes.h"

namespace valuation_functions {

/** Create a node representing a constant value  */
template<typename T>
tree::node::ptr create_const_node(const T& value) {
	tree::node::ptr	result = tree::node::ptr(new const_node<T>(value));
	return result;
}
;

/** Create a node representing the product of a constant with an expression */
template<typename T>
tree::node::ptr scalar_product_node(const T& value,tree::node::ptr p) {
	using namespace tree;
	node::ptr	a = create_const_node(value);
	node::ptr result = node::ptr(new arithmetic_node(MUL, a, p));
	return result;
}
;

/** Create a node representing the addition of two expressions */
inline
tree::node::ptr node_addition(tree::node::ptr p,tree::node::ptr q) {
	using namespace tree;
	node::ptr result = node::ptr(new arithmetic_node(ADD, p, q));
	return result;
}
;

/** Create a node representing the subtraction of two expressions */
inline
tree::node::ptr node_subtraction(tree::node::ptr p,tree::node::ptr q) {
	using namespace tree;
	node::ptr result = node::ptr(new arithmetic_node(SUB, p, q));
	return result;
}
;

/** Create a node representing the scalar product of a constant vector with a vector of variables
 *
 * The variables are represented by the domain.
 * Given a vector a1,...,an, a scalar b and a domain x1,...,xn, the result is
 * the expression a1*x1+...+an*xn+b.
 */
template<typename T>
tree::node::ptr scalar_product_node(const math::vector<T>& a, const T& b, const positional_vdomain& dom) {
	assert(dom.size()==a.size());
	using namespace tree;
	node::ptr b_node = create_const_node(b);
	node::ptr p = b_node;
	if (dom.size() > 0) {
		for (size_t j = 0; j < dom.size(); ++j) {
			size_t i = dom.size()-1-j;
			if (!math::numeric::is_MEQ(a[i], T(0))) {
				node::ptr xi = node::ptr(
						new variable_node(dom.get_variable(i)));
				node::ptr ai_xi = scalar_product_node(a[i], xi);
				p = node_addition(ai_xi, p);
			}
		}
	}
	return p;
}


/** Returns true if p is a const node */
template<typename scalar_type>
bool is_const_node(tree::node::ptr p) {
	const_node<scalar_type>* y = dynamic_cast<const_node<scalar_type>*>(p.get());
	return y;
}

/** Returns true if p is a const node with value zero */
template<typename scalar_type>
bool is_const_node_zero(tree::node::ptr p) {
	if (const_node<scalar_type>* y =
			dynamic_cast<const_node<scalar_type>*>(p.get())) {
		if (y->my_val == scalar_type(0))
			return true;
	}
	return false;
}

/** Returns true if p is a const node with value one */
template<typename scalar_type>
bool is_const_node_one(tree::node::ptr p) {
	if (const_node<scalar_type>* y =
			dynamic_cast<const_node<scalar_type>*>(p.get())) {
		if (y->my_val == scalar_type(1))
			return true;
	}
	return false;
}

/** Simplify the node non-recursively */
template<typename scalar_type>
tree::node::ptr simplify_one(tree::node::ptr p) {
//	return p;

	using namespace tree;
	if (arithmetic_node* a = dynamic_cast<arithmetic_node*>(p.get())) {
		// if both operands are constants, evaluate
		if (is_const_node<scalar_type>(a->child1)
				&& is_const_node<scalar_type>(a->child2)) {
			math::vdom_vector<scalar_type> vec;
			scalar_type val = arithmetic_eval_node(p, vec);
			return create_const_node<scalar_type>(val);
		}
		// if commutative right child is constant but left child is not, switch
		if (is_const_node<scalar_type>(a->child2)
				&& !is_const_node<scalar_type>(a->child1) && is_commutative(a->my_op)) {
			p = node::ptr(new arithmetic_node(a->my_op, a->child2, a->child1));
		}
		// bring up constants in commutative operators
		// ---------------------------------------
		// pull constants up
		// (c MUL x) MUL y = c MUL (x MUL y)
		if (arithmetic_node* b =
				dynamic_cast<arithmetic_node*>(a->child1.get())) {
			if (is_commutative(a->my_op) && (a->my_op==b->my_op)) {
				if (is_const_node<scalar_type>(b->child1)) {
					// a(b,y),b(c,x) -> a(c,b),b(x,y)
					node::ptr b_new = node::ptr(
							new arithmetic_node(a->my_op, b->child2,
									a->child2));
					b_new = simplify_one<scalar_type>(b_new);
					node::ptr result = node::ptr(
							new arithmetic_node(a->my_op, b->child1,
									b_new));
					return simplify_one<scalar_type>(result);
				}
			}
		}
		// x MUL (c MUL y) = c MUL (x MUL y)
		if (arithmetic_node* b = dynamic_cast<arithmetic_node*>(a->child2.get())) {
			if (is_commutative(a->my_op) && (a->my_op==b->my_op)) {
				if (is_const_node<scalar_type>(b->child1)) {
					if (is_const_node<scalar_type>(a->child1)) {
						// a(d,b),b(c,y) -> a(q,y),q(d,c)
						node::ptr q_new = node::ptr(
								new arithmetic_node(a->my_op, a->child1, b->child1));
						q_new = simplify_one<scalar_type>(q_new);
						node::ptr result = node::ptr(
								new arithmetic_node(a->my_op, q_new, b->child2));
						return simplify_one<scalar_type>(result);
					} else {
						// a(x,b),b(c,y) -> a(c,b),b(x,y)
						node::ptr b_new = node::ptr(
								new arithmetic_node(a->my_op, a->child1, b->child2));
						b_new = simplify_one<scalar_type>(b_new);
						node::ptr result = node::ptr(
								new arithmetic_node(a->my_op, b->child1, b_new));
						return simplify_one<scalar_type>(result);
					}
				}
			}
		}
		// ---------------------------------------
		if (a->my_op == NEG) {
			// evaluate neg const
			if (is_const_node<scalar_type>(a->child1)) {
				if (const_node<scalar_type>* x = dynamic_cast<const_node<
						scalar_type>*>(a->child1.get())) {
					return create_const_node<scalar_type>(-(x->my_val));
				}
			}
			// evaluate neg const in MUL
			if (arithmetic_node* b = dynamic_cast<arithmetic_node*>(a->child1.get())) {
				if (b->my_op == MUL && is_const_node<scalar_type>(b->child1)) {
					if (const_node<scalar_type>* y = dynamic_cast<const_node<
							scalar_type>*>(b->child1.get())) {
						node::ptr scalar = create_const_node<scalar_type>(-y->my_val);
						node::ptr result = node::ptr(new arithmetic_node(MUL, scalar, b->child2));
						return simplify_one<scalar_type>(result);
					}
				}
			}
		} else if (a->my_op == ADD) {
			// adding zero leaves other
			if (is_const_node_zero<scalar_type>(a->child1)) {
				return a->child2;
			} else if (is_const_node_zero<scalar_type>(a->child2)) {
				return a->child1;
			}
			// plus neg = SUB
			if (arithmetic_node* b = dynamic_cast<arithmetic_node*>(a->child2.get())) {
				if (b->my_op == NEG) {
					node::ptr result = node::ptr(new arithmetic_node(SUB, a->child1, b->child1));
					return simplify_one<scalar_type>(result);
				}
			}
		} else if (a->my_op == SUB) {
			// substracting zero leaves other, plus negation if necessary
			if (is_const_node_zero<scalar_type>(a->child2)) {
				return a->child1;
			} else if (is_const_node_zero<scalar_type>(a->child1)) {
				node::ptr result = node::ptr(new arithmetic_node(NEG, a->child2, node::ptr()));
				return simplify_one<scalar_type>(result);
			}
			// minus neg = ADD
			if (arithmetic_node* b = dynamic_cast<arithmetic_node*>(a->child2.get())) {
				if (b->my_op == NEG) {
					node::ptr result = node::ptr(new arithmetic_node(ADD, a->child1, b->child1));
					return simplify_one<scalar_type>(result);
				}
			}
		} else if (a->my_op == MUL) {
			// multiplication with one leaves other
			if (is_const_node_one<scalar_type>(a->child1)) {
				return a->child2;
			} else if (is_const_node_one<scalar_type>(a->child2)) {
				return a->child1;
			}
			// multiplication with zero gives zero
			if (is_const_node_zero<scalar_type>(a->child2)) {
				return create_const_node<scalar_type>(scalar_type(0));
			} else if (is_const_node_zero<scalar_type>(a->child1)) {
				return create_const_node<scalar_type>(scalar_type(0));
			}
		} else if (a->my_op == DIV) {
			// divison with one leaves other
			if (is_const_node_one<scalar_type>(a->child2)) {
				return a->child1;
			}
			// division of zero by anything but zero gives zero
			if (is_const_node_zero<scalar_type>(a->child1)) {
				if (const_node<scalar_type>* y = dynamic_cast<const_node<
						scalar_type>*>(a->child2.get())) {
					if (y->my_val != scalar_type(0)) {
						return create_const_node<scalar_type>(scalar_type(0));
					}
				}
			}
		} else if (a->my_op == POW) {
			// zero to the positive power of anything is zero
			if (is_const_node_zero<scalar_type>(a->child1)) {
				if (const_node<scalar_type>* y = dynamic_cast<const_node<
						scalar_type>*>(a->child2.get())) {
					if (y->my_val >= scalar_type(0)) {
						return create_const_node<scalar_type>(scalar_type(0));
					}
				}
			} else if (is_const_node_zero<scalar_type>(a->child2)) {
				// anything (except zero) to the power of zero is one
				return create_const_node<scalar_type>(scalar_type(1));
			} else if (is_const_node_one<scalar_type>(a->child2)) {
				return a->child1;
			}

		}
	}

	// nothing simplified, return original
	return p;
}

/** Simplify the node recursively
 *
 * @todo Currently a hack */
template<typename scalar_type>
tree::node::ptr simplify(tree::node::ptr p) {
	tree::node::ptr res = p;
	do {
		p = res;
		res = simplify_one<scalar_type>(res);
	} while (res != p && p);
	return res;
}

/** Differentiate a predicate according to a variable var */
template<typename scalar_type>
tree::node::ptr compute_diff(tree::node::ptr p, const variable& var) {
	using namespace tree;
//	using namespace predicate_parser;

// if null pointer then return null pointer
	if (p == tree::node::ptr()) {
		return p;
	}

	scalar_type one(1);
	node::ptr leftchild = node::ptr();
	node::ptr rightchild = node::ptr();
	node::ptr scalar_node = node::ptr();
	node::ptr result = node::ptr();

	if (arithmetic_node* a = dynamic_cast<arithmetic_node*>(p.get())) {
		if (a->my_op == NEG || a->my_op == ADD || a->my_op == SUB) {
			leftchild = compute_diff<scalar_type>(a->child1, var);
			rightchild = compute_diff<scalar_type>(a->child2, var);

			result = node::ptr(
					new arithmetic_node(a->my_op, leftchild, rightchild));
		} else if (a->my_op == MUL) {
			node::ptr x = simplify<scalar_type>(a->child1);
			node::ptr y = simplify<scalar_type>(a->child2);
			node::ptr dx = compute_diff<scalar_type>(a->child1, var);
			node::ptr dy = compute_diff<scalar_type>(a->child2, var);

			leftchild = node::ptr(new arithmetic_node(MUL, dx, y));
			rightchild = node::ptr(new arithmetic_node(MUL, x, dy));
			// simplify operands
			leftchild = simplify_one<scalar_type>(leftchild);
			rightchild = simplify_one<scalar_type>(rightchild);
			result = node::ptr(new arithmetic_node(ADD, leftchild, rightchild));
		} else if (a->my_op == DIV) {
			node::ptr x = simplify<scalar_type>(a->child1);
			node::ptr y = simplify<scalar_type>(a->child2);
			node::ptr dx = compute_diff<scalar_type>(a->child1, var);
			node::ptr dy = compute_diff<scalar_type>(a->child2, var);

			leftchild = node::ptr(new arithmetic_node(MUL, dx, y));
			rightchild = node::ptr(new arithmetic_node(MUL, x, dy));
			// simplify operands
			leftchild = simplify_one<scalar_type>(leftchild);
			rightchild = simplify_one<scalar_type>(rightchild);
			node::ptr num = node::ptr(
					new arithmetic_node(SUB, leftchild, rightchild));
			num = simplify_one<scalar_type>(num);
			scalar_node = create_const_node<scalar_type>(scalar_type(2));
			node::ptr den = node::ptr(new arithmetic_node(POW, y, scalar_node));
			den = simplify_one<scalar_type>(den);
			result = node::ptr(new arithmetic_node(DIV, num, den));
		} else if (a->my_op == POW) {
			// @note: only handle constant exponents or constant base
			// make sure they are simplified, otherwise we can't go on
			node::ptr x = simplify<scalar_type>(a->child1);
			node::ptr y = simplify<scalar_type>(a->child2);
			if (is_const_node<scalar_type>(y)) {
				// d/dx f(x)^b = b*f'*f(x)^(b-1)
				const_node<scalar_type>* cy =
						static_cast<const_node<scalar_type>*>(y.get());
				scalar_node = create_const_node<scalar_type>(
						cy->my_val - scalar_type(1));
				node::ptr dx = compute_diff<scalar_type>(a->child1, var);
				result = node::ptr(new arithmetic_node(MUL, y, dx));
				result = simplify_one<scalar_type>(result);
				node::ptr new_pow = node::ptr(
						new arithmetic_node(POW, x, scalar_node));
				new_pow = simplify_one<scalar_type>(new_pow);
				result = node::ptr(new arithmetic_node(MUL, result, new_pow));
			} else {
				// d/dx f(x)^g(x) = f^g(f'*g/f + g'*ln(f))
				std::cerr << y << std::endl;
				throw std::runtime_error(
						"can't handle variable exponents in differentiation");
			}
		} else
			throw std::runtime_error("unknown arithmetic node type");

		// simplify the result
		result = simplify_one<scalar_type>(result);
	} else if (variable_node* b = dynamic_cast<variable_node*>(p.get())) {
		// if p is a variable
		if (b->my_id == var.get_id())
			result = create_const_node<scalar_type>(scalar_type(1));
		else
			result = create_const_node<scalar_type>(scalar_type(0));
	} else if (dynamic_cast<const_node<scalar_type> *>(p.get())) {
		// if p is constant
		result = create_const_node<scalar_type>(scalar_type(0));
	}
	else
		throw std::runtime_error("unknown node type");

	//std::cout << "d/d"<<var<<"("<<p<<") = " << result << std::endl;
	return result;
}

}

#endif /* VALUATION_FUNCTION_TREE_OPS_H_ */
