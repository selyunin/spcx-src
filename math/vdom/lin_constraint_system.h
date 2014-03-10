#ifndef LIN_CONSTRAINT_SYSTEM_H_
#define LIN_CONSTRAINT_SYSTEM_H_

#include <list>
#include <boost/shared_ptr.hpp>
#include "math/vdom/lin_constraint.h"
#include "math/vdom/lin_constraint_system_visitor.h"
#include "utility/printable.h"

namespace math {

template<typename scalar_type> class lin_constraint_system:
		public printable {
public:
	typedef typename boost::shared_ptr<lin_constraint_system<scalar_type> > ptr;

	typedef typename boost::shared_ptr<const lin_constraint_system<scalar_type> >
			const_ptr;

	typedef typename std::list<
			typename math::lin_constraint<scalar_type> > my_list_type;

	typedef typename my_list_type::iterator
			iterator;

	typedef typename my_list_type::const_iterator
			const_iterator;

	/** Constructor that builds an empty list. */
	lin_constraint_system<scalar_type> () : up_to_date(true),my_unified(true) {
	}

	/** Constructor that makes a copy of an already existing list. */
	lin_constraint_system<scalar_type> (const my_list_type& cons) : my_list(cons),up_to_date(false),my_unified(false) {
	}
	;

	/** Copy constructor (deep copy). */
	lin_constraint_system<scalar_type> (const lin_constraint_system<scalar_type>& s) :
		my_list(s.my_list),my_vars(s.my_vars),up_to_date(s.up_to_date),my_unified(s.my_unified) {
	}
	;

	iterator begin() {
		up_to_date = false;
		my_unified = false;
		return my_list.begin();
	}
	;
	iterator end() {
		up_to_date = false;
		my_unified = false;
		return my_list.end();
	}
	;

	/** Erase constraint
	 *
	 * Returns an iterator pointing to the constraint that followed the erased constraint.
	 * If the list is empty, returns end().
	 */
	iterator erase(const iterator& it) {
		up_to_date = false;
		return my_list.erase(it);
	}
	;
	const_iterator begin() const {
		return my_list.begin();
	}
	;
	const_iterator end() const {
		return my_list.end();
	}
	;
	/** Insert the constraint p.
	 *
	 * This redirects to push_back. Its purpose is to provide an interface
	 * similar to other containers like sets. */
	void insert(const typename math::lin_constraint<scalar_type>& p) {
		push_back(p);
	}
	;
	/** Push back the constraint p. */
	void push_back(const typename math::lin_constraint<scalar_type>& p) {
		if (my_list.empty()) {
			my_vars = p.get_variable_ids();
			up_to_date = true;
			my_unified = true;
		} else {
			// compare to the last constraint
			if (my_list.rbegin()->get_l().get_index_to_variable_id_map()
					!= p.get_l().get_index_to_variable_id_map()) {
				up_to_date = false;
				my_unified = false;
			}
			if (!my_unified) {
				// to be safe, set up_to_date = false
				// @todo we could check my_vars against p.get_l().get_index_to_variable_id_map()
				// but that wouldn't be efficient. we should really store a domain instead
				// of my_vars
				up_to_date = false;
			}
		}
		my_list.push_back(p);
	}
	;
	/** Push back the constraints in p. */
	void push_back(const lin_constraint_system<scalar_type>& p) {
		if (up_to_date && p.up_to_date) {
			my_vars.insert(p.my_vars.begin(),p.my_vars.end());
			my_unified = my_unified && p.my_unified;
			if (!empty() && !p.empty()) {
				if (my_list.begin()->get_l().get_index_to_variable_id_map()
						!= p.my_list.begin()->get_l().get_index_to_variable_id_map()) {
					my_unified = false;
				}
			}
		} else {
			up_to_date = false;
			// @todo Strictly speaking, this seems conservative
			my_unified = false;
		}
		my_list.insert(my_list.end(),p.my_list.begin(),p.my_list.end());
	}
	;
	/** Push back the constraints in p. */
	void push_back(typename lin_constraint_system<scalar_type>::const_ptr p) {
		push_back(*p);
	}
	;
	/** Remove the last element. */
	void pop_back() {
		if (!my_unified) {
			up_to_date=false;
		}
		my_list->pop_back();
	}
	;
	/** Returns the number of constraints. */
	unsigned int size() const {
		return my_list.size();
	}
	;
	/** Returns whether the list of constraints is empty. */
	bool empty() const {
		return my_list.empty();
	}
	;
	/** Satisfiability check
	 *
	 * Returns false if the constraints are unsatisfiable, and true otherwise. */
	math::tribool is_satisfiable() const {
		math::tribool res=true;
		for (const_iterator it = begin(); it != end(); ++it) {
			res = res && it->is_satisfiable();
			if (res==false)
				return false;
		}
		return res;
	}
	;

	/** Returns whether the vector v satisfies all constraints. */
	math::tribool is_satisfied(const vdom_vector<scalar_type>& v) const {
		math::tribool sat(true);
		for(const_iterator it=begin();it!=end();++it){
			sat = sat && it->is_satisfied(v);
		}
		return sat;
	}
	;


	/** Bring all constraints in the system to the given domain. */
	void reorder(const positional_vdomain& dom) {
		if (!empty()) {
			for (typename lin_constraint_system<scalar_type>::iterator it =
					begin(); it != end(); ++it) {
				it->get_l().reorder(dom);
			}
		}
		my_vars = dom.get_variable_ids();
		up_to_date = true;
		my_unified = true;
	}
	;

	/** Bring all constraints in the system to a common domain and return the domain. */
	positional_vdomain unify_domains() {
		if (empty())
			return positional_vdomain();
		if (!my_unified) {
			//make sure my_vars is set
			update_cache();

			positional_vdomain dom(my_vars);
			reorder(dom);
			return dom;
		} else {
			return this->begin()->get_l().domain();
		}
	}
	;

	/** Bring all constraints in the system to a common domain. */
	bool has_unified_domains() const {
		return my_unified;
	}
	;

	/** Update the internal cache. */
	void update_cache() const {
		if (!up_to_date) {
			lin_constraint_system<scalar_type>* nonconst_this =
					const_cast<lin_constraint_system<scalar_type>*> (this);
			if (my_unified && !my_list.empty()) {
				nonconst_this->my_vars = my_list.begin()->get_l().get_variable_ids();
			} else {
				nonconst_this->my_vars = variable_id_set();
				if (!empty()) {
					// start with the first one
					typename lin_constraint_system<scalar_type>::const_iterator
							it = this->begin();
					nonconst_this->my_vars = it->get_variable_ids();
					positional_vdomain current_dom = it->get_l().domain();
					for (; it != this->end(); ++it) {
						if (current_dom
								!= it->get_l().domain()) {
							const variable_id_set& v = it->get_variable_ids();
							nonconst_this->my_vars.insert(v.begin(), v.end());
							current_dom = it->get_l().domain();
						}
					}
				}
			}
			nonconst_this->up_to_date = true;
		}
	}
	;

	/** Returns the variables in the constraints. */
	 const variable_id_set& get_variable_ids() const {
		update_cache();
		return my_vars;
	}
	;

	/** Remove variables from constraints.
	 *
	 * Remove variables and the associated coefficients from
	 * the constraints. Equivalent to setting the coefficients to zero,
	 * and removing the variables from the domain.
	 *
	 * @attention This is not existential quantification!
	 */
	 void remove_variables(const variable_id_set& vis) {
		if (up_to_date) {
			//my_vars.erase(vis.begin(),vis.end());
			set_difference_assign(my_vars,vis);
		}
		for (iterator it = begin(); it != end(); ++it) {
			it->remove_variables(vis);
		}
		// this could now possibly be unified
	}
	;

	/** Convert to a list of constraints of type result_type. */
	template<typename result_type> lin_constraint_system<result_type> convert_to() const {
		lin_constraint_system<result_type> res;
		for (typename lin_constraint_system<scalar_type>::const_iterator it =
				this->begin(); it != this->end(); ++it) {
			res.push_back(it->template convert_to<result_type> ());
		}
		return res;
	}
	;

	/** Convert equalities to two inequality constraints. */
	void expand_equalities() {
		for (typename my_list_type::iterator it = my_list.begin(); it
				!= my_list.end(); ++it) {
			if (it->is_equality()) {
				typename math::lin_constraint<scalar_type> con = *it;
				con.set_sign(LE);
				it->set_sign(GE);
				my_list.push_front(con);
			}
		}
	}
	/**
	 * Collapses inequalities to equality constraint when possible.
	 */
	void collapse_inequalities() {
		my_list_type collapsed_list;
		comparison_operator cons_sign;

		for (typename my_list_type::iterator it = my_list.begin(); it
			!= my_list.end(); ++it){
			cons_sign = it->get_canonic_sign();
			if(cons_sign != LE){
				//typename math::lin_constraint<scalar_type> con = *it;
				collapsed_list.push_back(*it);
				continue;
			}
			else{
				typename my_list_type::iterator iter = it;
				for(++iter;iter!=my_list.end();++iter){
					if(iter->get_canonic_sign() == LE &&
							it->get_canonic_l() == -iter->get_canonic_l()){
						iter = my_list.erase(iter);
						it->set_sign(EQ);
						break;
					}
				}
				collapsed_list.push_back(*it);
			}
		}
		my_list = collapsed_list;
		// in case we lost variables from domains,
		// let's set up_to_date to false unless
		// they all had the same domain
		if (!my_unified)
			up_to_date = false;
	}
	;

	/** Process the visitor v. */
	void accept(lin_constraint_system_visitor<scalar_type>& v) const {
		v.visit(*this);
	}
	;

	void print(std::ostream& os) const {
		for (typename my_list_type::const_iterator it = my_list.begin(); it
				!= my_list.end(); ++it) {
			if (it != my_list.begin())
				os << " & ";
			os << *it;
		}
	}
	;

private:
	my_list_type my_list;
	variable_id_set my_vars;
	bool up_to_date; // true if my_vars = union of the variables in constraints
	bool my_unified; // true if all constraints are over the same domain
};

}

#endif /*LIN_CONSTRAINT_SYSTEM_H_*/
