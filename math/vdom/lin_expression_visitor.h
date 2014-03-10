#ifndef LIN_EXPRESSION_VISITOR_H_
#define LIN_EXPRESSION_VISITOR_H_

namespace math {

/** Forward declaration. */
template <typename scalar_type> class lin_expression;

/** Visitor for processing lin_expression objects */
template <typename scalar_type> class lin_expression_visitor {
public:
	typedef unsigned int index_type;
	virtual ~lin_expression_visitor() {
	}
	;
	/** visit is called by the lin_expression's accept. */
	virtual void visit(const lin_expression<scalar_type>& l) {
		lin_expression_prologue(l);
		visit_coeffs(l);
		lin_expression_epilogue(l);
	}
	;
	/** The prologue is called before processing the linear expression. */
	virtual void lin_expression_prologue(const lin_expression<scalar_type>& l) {
	}
	;
	/** If it returns true the visit_coeff function will be called for every coefficient. */
	virtual void visit_coeffs(const lin_expression<scalar_type>& l) {
		for (unsigned int i=0; i<l.size(); i++) {
			visit_coeff(l, l[i], i);
		}
		visit_inh_coeff(l, l.get_inh_coeff());
	}
	;
	/** The visit_coeff is called for every coefficient with its index i and the iimap of the lin_expression. */
	virtual void visit_coeff(const lin_expression<scalar_type>& l, const scalar_type& coeff,
			const index_type& i) {
	}
	;
	/** The visit_inh_coeff is called for the inhomogeneous coefficient b. */
	virtual void visit_inh_coeff(const lin_expression<scalar_type>& l, const scalar_type& coeff) {
	}
	;
	/** The epilogue is called after processing the linear expression. */
	virtual void lin_expression_epilogue(const lin_expression<scalar_type>& l) {
	}
	;
};

}

#endif /*LIN_EXPRESSION_VISITOR_H_*/
