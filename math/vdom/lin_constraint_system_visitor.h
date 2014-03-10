/*
 * lin_constraint_system_visitor.h
 *
 *  Created on: Dec 23, 2009
 *      Author: frehse
 */

#ifndef LIN_CONSTRAINT_SYSTEM_VISITOR_H_
#define LIN_CONSTRAINT_SYSTEM_VISITOR_H_

namespace math {

/* Forward declaration. */
template<typename scalar_type> class lin_constraint_system;

/** Visitor for processing lin_constraint objects */
template<typename scalar_type> class lin_constraint_system_visitor: public math::lin_constraint_visitor<
		scalar_type> {
public:
	virtual ~lin_constraint_system_visitor() {
	}
	;
	/** visit calls the prologue, then accept for every constraint, and then the epiogue. */
	virtual void visit(const lin_constraint_system<scalar_type>& con_sys);
	/** The prologue is called before processing the linear expression. */
	virtual void lin_constraint_system_prologue() {
	}
	;
	/** The epilogue is called after processing the linear expression. */
	virtual void lin_constraint_system_epilogue() {
	}
	;
};

/* -----------------------------------------------------------
 * Implementations
 * ----------------------------------------------------------- */

template<typename scalar_type> void lin_constraint_system_visitor<scalar_type>::visit(
		const lin_constraint_system<scalar_type>& con_sys) {
	lin_constraint_system_prologue();
	for (typename lin_constraint_system<scalar_type>::const_iterator it =
			con_sys.begin(); it != con_sys.end(); ++it) {
		it->accept(*this);
	}
	lin_constraint_system_epilogue();
}
;

}

#endif /* LIN_CONSTRAINT_SYSTEM_VISITOR_H_ */
