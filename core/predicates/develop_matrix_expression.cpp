/*
 * develop_matrix_expression.cpp
 *
 *  Created on: Oct 12, 2010
 *      Author: gvincent
 */

#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/predicates/node_print_visitor.h"

#include "math/matrix.h"

namespace valuation_functions {

math::matrix<tree::node::ptr> develop_matrix(const tree::node::ptr& p) {

	if (arithmetic_node* q = dynamic_cast<arithmetic_node*>(p.get())) {
		math::matrix<tree::node::ptr> child1(develop_matrix(q->child1));
		math::matrix<tree::node::ptr> child2(develop_matrix(q->child2));
		/**
		 * matrix negation
		 */
		if (q->my_op == NEG) {

			math::matrix<tree::node::ptr> new_mat(child1.size1(),
					child1.size2(), tree::node::null_node());
			for (int i = 0; i < new_mat.size1(); ++i) {
				for (int j = 0; j < new_mat.size2(); ++j) {
					new_mat(i, j) = tree::node::ptr(new arithmetic_node(NEG,
							child1(i, j)));
				}
			}
			return new_mat;
		}

		/**
		 * matrix addition or subtraction
		 */
		else if (q->my_op == ADD || q->my_op == SUB) {
			if (child1.size1() != child2.size1() || child1.size2()
					!= child2.size2()) {
				std::stringstream ss;
				ss << q->child1 << " vs " << q->child2;
//				valuation_functions::print_node(ss,q->child1);
//				ss << " vs ";
//				valuation_functions::print_node(ss,q->child2);
				throw std::runtime_error(
						"Incorrect matrix dimensions in matrix addition or subtraction: ("
								+ to_string(child1.size1()) + "," + to_string(
								child1.size2()) + ") vs (" + to_string(
								child2.size1()) + "," + to_string(
								child2.size2()) + "):\n Expressions:"+ss.str());
			}
			math::matrix<tree::node::ptr> new_mat(child1.size1(),
					child1.size2(), tree::node::null_node());
			for (int i = 0; i < new_mat.size1(); ++i) {
				for (int j = 0; j < new_mat.size2(); ++j) {
					new_mat(i, j) = tree::node::ptr(new arithmetic_node(
							q->my_op, child1(i, j), child2(i, j)));
				}
			}
			return new_mat;
		}

		else if (q->my_op == MUL) {
			/**
			 * Multiplication by a scalar.
			 */
			if (child1.size1() == 1 && child1.size2() == 1) {
				math::matrix<tree::node::ptr> new_mat(child2.size1(),
						child2.size2(), tree::node::null_node());
				for (int i = 0; i < new_mat.size1(); ++i) {
					for (int j = 0; j < new_mat.size2(); ++j) {
						new_mat(i, j) = tree::node::ptr(new arithmetic_node(
								MUL, child1(0, 0), child2(i, j)));
					}
				}
				return new_mat;
			}
			if (child2.size1() == 1 && child2.size2() == 1) {
				math::matrix<tree::node::ptr> new_mat(child1.size1(),
						child1.size2(), tree::node::null_node());
				for (int i = 0; i < new_mat.size1(); ++i) {
					for (int j = 0; j < new_mat.size2(); ++j) {
						new_mat(i, j) = tree::node::ptr(new arithmetic_node(
								MUL, child1(i, j), child2(0, 0)));
					}
				}
				return new_mat;
			}

			/**
			 * Multiplication by another matrix.
			 */

			if (child1.size2() != child2.size1())
				throw std::runtime_error(
						"Incorrect matrix dimensions in matrix multiplication");
			math::matrix<tree::node::ptr> new_mat(child1.size1(),
					child2.size2(), tree::node::null_node());

			for (int i = 0; i < new_mat.size1(); ++i) {
				for (int j = 0; j < new_mat.size2(); ++j) {
					new_mat(i, j) = tree::node::ptr(new arithmetic_node(MUL,
							child1(i, 0), child2(0, j)));
					for (int n = 1; n < child1.size2(); ++n) {
						new_mat(i, j) = tree::node::ptr(new arithmetic_node(
								ADD, new_mat(i, j), tree::node::ptr(
										new arithmetic_node(MUL, child1(i, n),
												child2(n, j)))));
					}
				}
			}
			return new_mat;
		}
		else if (q->my_op == DIV) {
			/**
			 * Division (only for scalar)
			 */
			if (child1.size1() == 1 && child1.size2() == 1 && child1.size2()
					== 1 && child2.size2() == 1) {
				math::matrix<tree::node::ptr> new_mat(1, 1,
						tree::node::null_node());
				new_mat(0, 0) = tree::node::ptr(new arithmetic_node(DIV,
						child1(0, 0), child2(0, 0)));
				return new_mat;
			} else
				throw std::runtime_error("Division of matrix not implemented");

		}

		else if (q->my_op == POW) {
		/**
		 * Power (only to scalar) - To be implemented properly
		 */
			if (child1.size1() == 1 && child1.size2() == 1 && child1.size2()
				== 1 && child2.size2() == 1) {
			math::matrix<tree::node::ptr> new_mat(1, 1,
					tree::node::null_node());
			new_mat(0, 0) = tree::node::ptr(new arithmetic_node(POW,
					child1(0, 0), child2(0, 0)));
					return new_mat;
			} else
					throw std::runtime_error("Power to the matrix not implemented");
		}

		else if(q->my_op == SQRT || q->my_op == LOG || q->my_op == SIN || q->my_op == EXP
					|| q->my_op == COS || q->my_op == TAN) {

			math::matrix<tree::node::ptr> new_mat(child1.size1(),
								child1.size2(), tree::node::null_node());
			for (int i = 0; i < new_mat.size1(); ++i) {
				for (int j = 0; j < new_mat.size2(); ++j) {
				new_mat(i, j) = tree::node::ptr(new arithmetic_node(q->my_op,
						child1(i, j)));
				}
			}
			return new_mat;
		}

		else
			throw std::runtime_error("Unknown arithmetic operator" + q->my_op);
	}

	/**
	 * If node is a variable.
	 */
	else if (variable_node* q = dynamic_cast<variable_node*>(p.get())) {
		int dim1 = q->get_dim1();
		int dim2 = q->get_dim2();
		/**
		 * If variable is not a vector, just return a matrix(1, 1) with this variable;
		 */
		if (dim1 <= 1 && dim2 <= 1) {
			return math::matrix<tree::node::ptr>(1, 1, p);
		}
		/**
		 * Else, return a matrix of variables nodes with components of this vector variable.
		 */
		math::matrix<tree::node::ptr> new_mat(dim1, dim2,
				tree::node::null_node());
		variable var(q->my_id);
//std::cout << "expanding variable "<<variable(q->my_id)+" to matrix";
		int n = 1;
		for (int i = 0; i < dim1; ++i) {
			for (int j = 0; j < dim2; ++j) {
				new_mat(i, j) = tree::node::ptr(new variable_node(
						var.get_element(n)));
				n++;
			}
		}
		return new_mat;
	}

	/**
	 * If node is a constant node with rational_matrix, return matrix of constant rational node.
	 * Only matrix of rational are implemented !
	 */
	else if (const_node<math::matrix<Rational> >* q = dynamic_cast<const_node<math::matrix<Rational> >*>(p.get())) {
		math::matrix<tree::node::ptr> new_mat(q->my_val.size1(),
				q->my_val.size2(), tree::node::null_node());

		for (int i = 0; i < new_mat.size1(); ++i) {
			for (int j = 0; j < new_mat.size2(); ++j) {
				new_mat(i, j) = tree::node::ptr(new const_node<Rational> (
						q->my_val(i, j)));
			}
		}
		return new_mat;
	}
	/**
	 * For other constant nodes, just encapsulate in matrix
	 */
	else {
		return math::matrix<tree::node::ptr>(1, 1, p);
	}
}

tree::node::ptr develop_assignment_matrix_tree(const tree::node::ptr& p) {
	tree::node::ptr new_node = p;

	if (boolean_node* q = dynamic_cast<boolean_node*>(p.get())) {
		tree::node::ptr child1 = develop_assignment_matrix_tree(q->child1);
		if (q->my_op == NOT) {
			new_node = tree::node::ptr(new boolean_node(NOT, child1));
		} else {
			tree::node::ptr child2 = develop_assignment_matrix_tree(q->child2);
			new_node = tree::node::ptr(new boolean_node(q->my_op, child1,
					child2));
		}
	}
	/**
	 * Develop vector comparisons
	 */
	else if (comparison_node* q = dynamic_cast<comparison_node*>(p.get())) {

					math::matrix<tree::node::ptr> child1(develop_matrix(
												q->child1));
					math::matrix<tree::node::ptr> child2(develop_matrix(
							q->child2));
					if (child1.size2() != 1 && child2.size2())
						throw std::runtime_error(
								"Matrix comparison is not implemented (only vector comparison)");

					if (child1.size1() != child2.size1())
						throw std::runtime_error(
								"Incorrect vector dimensions in vector variable assignment");

					new_node = tree::node::ptr(new comparison_node(q->my_op,
							child1(0, 0), child2(0, 0)));

					for (int n = 1; n < child2.size1(); ++n) {
						new_node = tree::node::ptr(new boolean_node(AND,
								new_node, tree::node::ptr(new comparison_node(
										q->my_op, child1(n, 0),	child2(n, 0)))));
					}
	}

	return new_node;
}

}
