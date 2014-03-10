/*
 * affine_map.h
 *
 *  Created on: Sep 24, 2009
 *      Author: frehse
 */

#ifndef AFFINE_MAP_H_
#define AFFINE_MAP_H_

#include "utility/tree_node.h"
#include "math/vdom/lin_expression.h"
#include "math/matrix.h"
//#include "math/vdom/index_to_variable_id_map_provider.h"
#include "math/vdom/positional_vdomain.h"
#include "math/vdom/vdom_matrix.h"
#include "math/vdom/vdom_matrix_utility.h"
#include "math/vdom/vdom_matrix_operators.h"
#include "math/vdom/vdom_vector.h"
#include "math/vdom/vdom_vector_utility.h"
#include "math/vdom/vdom_vector_operators.h"
#include "math/numeric/approx_comparator.h"

namespace math {

/** A class for representing an affine map of the form x:=A*x+b,
 * where A is a matrix and b is a vector of corresponding size.
 *
 * Internally, an empty A matrix represents the identity map.
 * A call to get_A() replaces it with a diagonal matrix to ensure
 * consistency. For optimal performance, test with if (!is_translation_coded()) before
 * using get_A().
 *  */
template<typename scalar_type> class affine_map: public virtual printable {
private:
	enum def_type {
		A_b_defined, translation, void_defined, universe_defined
	};
public:
	typedef boost::shared_ptr<affine_map<scalar_type> > ptr;
	typedef boost::shared_ptr<const affine_map<scalar_type> > const_ptr;
	typedef matrix<scalar_type> matrix_type;
	typedef vector<scalar_type> vector_type;
	typedef vdom_matrix<scalar_type> vdom_matrix_type;
	typedef vdom_vector<scalar_type> vdom_vector_type;
	typedef typename matrix_type::size_type size_type;
	typedef tree::node predicate_type;

	/** Define an empty map. */
	affine_map() :
		my_def(A_b_defined) {
	}
	;

	/** Define a map where y'=a^Tx+b and all other variables in iimap remain constant. */
	affine_map(const variable& y, const vdom_vector_type& a,
			const scalar_type& b) :
		my_def(A_b_defined) {
		my_A = diagonal_matrix(a.domain(), scalar_type(1));
		my_b = vdom_vector_type(a.domain(), scalar_type(0));
		assign_row(y, a, b);
	}
	;
	/** Define the map x'=x+b, b a scalar, over a given domain. */
	affine_map(positional_vdomain dom, const scalar_type& b) :
		my_def(translation) {
		my_b = vdom_vector_type(dom, b);
	}
	;
	/** Define the map x'=x+b, b a vdom_vector, over the domain dom.
	 *
	 * If the dom is different from the domain of b, b gets reordered
	 * to match dom. */
	affine_map(const positional_vdomain& dom, const vdom_vector_type& b) :
		my_def(translation) {
		my_b = b;
		if (dom != b.domain()) {
			my_b.reorder(dom);
		}
	}
	;

	/** Define the map x'=Ax, where A is a matrix.
	 *
	 * If the codomain of A is not equal to the domain of b, they
	 * are remapped to a common domain. */
	explicit affine_map(const vdom_matrix_type& A);

	/** Define the map x'=Ax+b, where A is a matrix and b is a vector.
	 *
	 * If the codomain of A is not equal to the domain of b, they
	 * are remapped to a common domain. */
	affine_map(const vdom_matrix_type& A, const vdom_vector_type& b);

	/** Define the map x'=Ax+b, where A is a square matrix and b is a vector.
	 */
	affine_map(const positional_vdomain& dom, const matrix_type& A,
			const vector_type& b) :
		my_def(A_b_defined) {
		my_A = math::vdom_matrix<scalar_type>(dom, dom, A);
		my_b = math::vdom_vector<scalar_type>(dom, b);
	}
	;

	virtual ~affine_map() {
	}
	;

	/** Obtain the domain. */
	const positional_vdomain& domain() const {
		if (is_translation_coded())
			return my_b.domain();
		else
			return my_A.domain();
	}
	;

	/** Obtain the codomain. */
	const positional_vdomain& codomain() const {
		if (is_translation_coded())
			return my_b.domain();
		else
			return my_A.codomain();
	}
	;

	/** Returns the map in predicate form */
	virtual predicate_type::ptr get_predicate() const {
		throw std::runtime_error(
				"missing implementation of get_predicate for affine_map");
	}
	;

	/** Returns true if all variables in the codomain obtain the value zero.
	 *
	 * @note By the above definition, the void map is considered zero.
	 * The universe map is not considered zero. */
	virtual bool is_zero() const {
		if (is_universe()) {
			return false;
		} else {
			return my_A.is_zero() && my_b.is_zero();
		}
	}
	;

	/** Returns whether the map is void. */
	virtual bool is_empty() const {
		return my_A.size1() == 0 && my_b.size() == 0;
	}
	;

	/** Returns whether the map is void.
	 *
	 * A map is void if it maps a non-empty domain to the empty set or if it corresponds to the predicate "false". */
	virtual bool is_void() const {
		//return my_A.size1() == 0 && my_A.size2() != 0;
		return my_def == void_defined;
	}
	;

	/** Returns whether the map is universe.
	 *
	 * A map is universe if it corresponds to the predicate "true". */
	virtual bool is_universe() const {
		//return my_A.size1() == 0 && my_A.size2() != 0;
		return my_def == universe_defined;
	}
	;

	/** Create a map that maps the domain dom to the empty set. */
	static affine_map<scalar_type> void_map(const positional_vdomain& codom = positional_vdomain(), const positional_vdomain& dom = positional_vdomain()) {
		vdom_matrix_type A(codom, dom);
		vdom_vector_type b(codom);
		affine_map<scalar_type> M(A, b);
		M.my_def = void_defined;
		return M;
	}
	;

	/** Create a map that maps the domain dom to the universe set. */
	static affine_map<scalar_type> universe_map(const positional_vdomain& codom = positional_vdomain(), const positional_vdomain& dom = positional_vdomain()) {
		vdom_matrix_type A(codom, dom);
		vdom_vector_type b(codom);
		affine_map<scalar_type> M(A, b);
		M.my_def = universe_defined;
		return M;
	}
	;

	/** Create a map that projects a vector from the domain dom onto the codomain codom.
	 *
	 * codom must be a subset of dom. */
	static affine_map<scalar_type> projection_map(const positional_vdomain& codom, const positional_vdomain& dom) {
		vdom_matrix_type A(codom, dom);
		vdom_vector_type b(codom);
		for (size_type i=0;i<codom.size();++i) {
			size_type row_index = dom.pos(codom.get_variable(i));
			A(i,row_index)=scalar_type(1);
		}
		affine_map<scalar_type> M(A, b);
		return M;
	}
	;

	/** Assign variable x to be mapped to x_vid'=a^Tx+b, without modifying the other variables. */
	void assign_row(const variable& x, const vdom_vector_type& a,
			const scalar_type& b) {
		instantiate_identity_A();
		typename vdom_matrix_type::size_type x_pos;
		x_pos = my_A.assign_row(x, a);
		my_b[x_pos] = b;
	}
	;

	/** Add a constant vector to the map */
	void operator+=(const vdom_vector_type& b) {
		if (b.domain()!=codomain()) {
			vdom_vector_type bb = b;
			bb.reorder(codomain());
			my_b+=bb;
		} else {
			my_b+=b;
		}
	}
	;

	/** Reorder the elements of the map to obtain a given codomain and domain.
	 *
	 * Throws if there is a nonzero coeffients for a variable
	 * that is not in the domain or codomain,
	 * unless the coefficent concerns the same variable in both domains.
	 * If dom or codom have new variables, the
	 * corresponding coefficients are set to zero.
	 * */
	virtual void reorder(const positional_vdomain& codom,
			const positional_vdomain& dom) {
		instantiate_identity_A();
		my_A.reorder(codom, dom);
		my_b.reorder(codom);
	}

	/**
	 * Converts the elements of the vector to the type result_type.
	 * Uses the template function convert_element, which can be specialized to
	 * implement specific conversions between types.
	 */
	template<typename result_type> affine_map<result_type> convert_to() const {
		typename affine_map<result_type>::vdom_matrix_type A =
				my_A.template convert_to<result_type> ();
		typename affine_map<result_type>::vdom_vector_type b =
				my_b.template convert_to<result_type> ();
		return affine_map<result_type> (A, b);
	}
	;

	/** Returns the ids of all variables in the domain or codomain. */
	variable_id_set get_variable_ids() const {
		variable_id_set vis = domain().get_variable_ids();
		variable_id_set vis2 = codomain().get_variable_ids();
		vis.insert(vis2.begin(), vis2.end());
		return vis;
	}

	/** Output the map in the form
	 * x1' == a11*x1+...+a1n*xn+b1 &
	 * ...
	 * xn' == an1*x1+...+ann*xn+bn
	 */
	virtual void print(std::ostream& os) const {
		if (is_translation_coded()) {
			for (unsigned int i = 0; i < my_b.size(); ++i) {
				if (i > 0) {
					os << " & ";
				}
				os << variable(
						variable::get_id_primedness_increased(
								codomain().get_variable(i).get_id()));
				os << " == ";
				// Use the printer of linear_expressions to make things easier
				vdom_vector<scalar_type> unit_vector(domain());
				unit_vector[i] = scalar_type(1);
				math::lin_expression<scalar_type> le(unit_vector, my_b[i]);
				os << le;
			}
		} else {
			for (unsigned int i = 0; i < my_A.size1(); ++i) {
				if (i > 0) {
					os << " & ";
				}
				os << variable(
						variable::get_id_primedness_increased(
								codomain().get_variable(i).get_id()));
				os << " == ";
				// Use the printer of linear_expressions to make things easier
				math::lin_expression<scalar_type> le(
						my_A.vdom_vector_from_row(i), my_b[i]);
				os << le;
			}
		}
	}
	;

	const vdom_matrix_type& get_A() const {
		instantiate_identity_A();
		return my_A;
	}
	;
	const vdom_vector_type& get_b() const {
		return my_b;
	}
	;

	/** Compute used and modified variables according to A and b. */
	void compute_used_and_modif_variables(variable_id_set& used_vars,
			variable_id_set& modif_vars) {
		if (is_translation_coded()) {
			used_vars = domain().get_variable_ids();
			modif_vars = codomain().get_variable_ids();
		} else {
			used_vars = variable_id_set();
			modif_vars = variable_id_set();
			for (unsigned int i = 0; i < my_A.size1(); ++i) {
				bool non_ident = false;
				for (unsigned int j = 0; j < my_A.size2(); ++j) {
					if (!math::numeric::approx_comparator<scalar_type>::is_maybe_zero(
							my_A(i, j))) {
						used_vars.insert(my_A.domain().get_variable(j).get_id());
						// nonzero coeff only allowed if i==j and it has to be 1
						if (i != j || !math::numeric::approx_comparator<
								scalar_type>::is_maybe_equal(my_A(i, j),
								scalar_type(1)))
							non_ident = true;
					}
				}
				if (non_ident
						|| !math::numeric::approx_comparator<scalar_type>::is_maybe_zero(
								my_b[i]))
					modif_vars.insert(my_A.codomain().get_variable(i).get_id());
			}
		}
	}
	;

	/** Compute A*x+b for a given vector x. */
	vdom_vector<scalar_type> map(const vdom_vector<scalar_type>& x) const;

	/** Returns true if the map is a pure translation
	 *
	 * The map is a translation if A is a diagonal matrix. */
	bool is_translation() const {
		//return (my_A.size1()==0 && my_A.size2()==0);
		if (my_def == translation)
			return true;
		else
			return my_A.size1()==my_A.size2() && is_diagonal_matrix<scalar_type>(my_A,scalar_type(1));
	}

	/** Returns true if the map is represented internally without A matrix */
	bool is_translation_coded() const {
		//return (my_A.size1()==0 && my_A.size2()==0);
		return (my_def == translation);
	}


private:

	void instantiate_identity_A() const {
		if (is_translation_coded()) {
			// this operation can be hidden from the user, so it needs to be const
			affine_map<scalar_type>* nonconst_this = const_cast<affine_map<
					scalar_type>*> (this);
			nonconst_this->my_A = diagonal_matrix(domain(), scalar_type(1));
			set_def_type(A_b_defined);
		}
	}
	;
	void set_def_type(def_type t) const {
		affine_map<scalar_type>* nonconst_this = const_cast<affine_map<
				scalar_type>*> (this);
		nonconst_this->my_def = A_b_defined;
	}
	;

private:
	vdom_matrix_type my_A;
	vdom_vector_type my_b;

	def_type my_def;
};

}

#include "affine_map_utility.h"

namespace math {

template<typename scalar_type> affine_map<scalar_type>::affine_map(
		const vdom_matrix_type& A) :
	my_def(A_b_defined) {
	my_A = A;
	my_b = vdom_vector_type(A.codomain());
}

template<typename scalar_type> affine_map<scalar_type>::affine_map(
		const vdom_matrix_type& A, const vdom_vector_type& b) :
	my_def(A_b_defined) {
	my_A = A;
	my_b = b;
	if (my_A.codomain() != my_b.domain()) {
		map_to_common_domain(my_A, my_b);
	}
}

template<typename scalar_type>
vdom_vector<scalar_type> affine_map<scalar_type>::map(
		const vdom_vector<scalar_type>& x) const {
	if (is_translation_coded()) {
		return x + my_b;
	} else {
		return my_A * x + my_b;
	}
}

}

#endif /* AFFINE_MAP_H_ */
