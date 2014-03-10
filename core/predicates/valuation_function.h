#ifndef VALUATION_FUNCTION_H_
#define VALUATION_FUNCTION_H_

#include <map>
#include <string>
#include <iostream>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
//#include "../../utility/variable.h"
#include "core/predicates/valuation.h"

namespace valuation_functions {

/** \brief A general interface for functions that evaluate valuations.
 * A valuation_function takes a valuation of type valuation_type and returns
 * a value of type eval_type. */
template<typename eval_type, typename valuation_type, typename codomain_type=variable_id> class valuation_function :
	public shared_ptr_user<valuation_function<eval_type,valuation_type> > {
public:
	valuation_function() {
	}
	;
	virtual ~valuation_function() {
	}
	;
	virtual eval_type
			eval(const typename variable_valuation<valuation_type>::const_ptr& v) const = 0;

//	virtual eval_type
//			print(const typename variable_valuation<valuation_type>::const_ptr& v, std::ostream& os) const = 0;


	/** Returns the variables used in the evaluation of the function. */
	virtual variable_id_set get_variable_ids() const = 0;
};

/** An assignment attributes a valuation_function to a variable. */
template<typename eval_type, typename valuation_type> class assignment {
private:
	std::string my_variable;
	typename valuation_function<eval_type,valuation_type>::ptr my_function_ptr;

public:
	assignment(
			const std::string& x,
			const typename valuation_function<eval_type,valuation_type>::ptr& f_ptr) {
		my_variable = x;
		my_function_ptr = f_ptr;
	}
	;
	assignment() {
	}
	;
	virtual ~assignment() {
	}
	;

	std::string get_assigned_variable() {
		return my_variable;
	}
	;
	const typename valuation_function<eval_type,valuation_type>::ptr& get_function() {
		return my_function_ptr;
	}
	;

};


/* An assignment attributes a valuation_function to a variable.
template<typename eval_type, typename valuation_type> class assignment {
private:
	//string my_variable;
	tree::node::ptr my_variable;
	typename valuation_function<eval_type,valuation_type>::ptr my_function_ptr;

public:
	assignment(
			const string& x,
			const typename valuation_function<eval_type,valuation_type>::ptr& f_ptr) {
		my_variable = tree::node::ptr(new variable_node(x));
		my_function_ptr = f_ptr;
	}
	;
	assignment(const tree::node::ptr& p, const typename valuation_function<eval_type,valuation_type>::ptr& f_ptr) {
		if (variable_node* q = dynamic_cast<variable_node*>(p.get())){
				my_variable = tree::node::ptr(new variable_node(q->my_id));
				my_function_ptr = f_ptr;
		}
	}
	;
	assignment() {
	}
	;
	virtual ~assignment() {
	}
	;

	string get_assigned_variable() {
		return my_variable;
	}
	;
	const typename valuation_function<eval_type,valuation_type>::ptr& get_function() {
		return my_function_ptr;
	}
	;

};*/

}

#endif

