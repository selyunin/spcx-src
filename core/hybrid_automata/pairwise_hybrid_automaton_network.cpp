/*
 * pairwise_hybrid_automaton_network.cpp
 *
 *  Created on: Aug 26, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/pairwise_hybrid_automaton_network.h"

#include <stdexcept>

#include "utility/basic_exception.h"
#include "utility/shared_ptr_output.h"
#include "utility/stl_helper_functions.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/location_constraint_set.h"
#include "core/hybrid_automata/hybrid_automaton_pair.h"
#include "core/predicates/dot_context_lookup.h"

namespace hybrid_automata {

/*
class pw_network_node {
public:
	virtual ~pw_network_node() {
	}
	;
	virtual pw_network_node* clone() = 0;
	virtual location_id
			get_location_id(const location_constraint_set& lcons) const = 0;
};

class pw_network_nonpair_node: public pw_network_node {
public:
	pw_network_nonpair_node(const hybrid_automaton::const_ptr& p) :
		my_aut(p) {
	}
	;
	pw_network_node* clone() {
		return new pw_network_nonpair_node(my_aut);
	}
	;
	// Return the id of the complete set lcons
	location_id get_location_id(const location_constraint_set& lcons) const {
		location_constraint_set::const_iterator_pair itp =
				lcons.get_constraints(my_aut->get_id());
		if (itp.first == itp.second) {
			std::cerr << lcons;
			throw std::runtime_error("no constraint found for aut "
					+ int2string(my_aut->get_id()));
		}
		return itp.first->second.get_id();
	}
	;
private:
	hybrid_automaton::const_ptr my_aut;
};

class pw_network_pair_node: public pw_network_node {
public:
	pw_network_pair_node(const hybrid_automaton_pair::const_ptr& p,
			pw_network_node* c1, pw_network_node* c2) :
		my_aut(p), child1(c1), child2(c2) {
		assert(child1);
		assert(child2);
	}
	;
	~pw_network_pair_node() {
		delete child1;
		delete child2;
	}
	;
	// Copy constructor needs to make a deep copy.
	pw_network_pair_node(const pw_network_pair_node& n) {
		my_aut = n.my_aut;
		child1 = n.child1->clone();
		child2 = n.child2->clone();
	}
	;
	pw_network_node* clone() {
		return new pw_network_pair_node(my_aut, child1, child2);
	}
	;

	location_id get_location_id(const location_constraint_set& lcons) const {
		location_id lid;
		location_id l1 = child1->get_location_id(lcons);
		location_id l2 = child2->get_location_id(lcons);
		//		if (!my_aut->find_location_id(lid, l1, l2))
		//			std::cerr << lcons;
		//			std::cerr << " resolved to " << "(" << l1 << "," << l2 << ")" << std::endl;
		//			my_aut->print(std::cerr);
		//			throw std::runtime_error("could not resolve location constraint in pairwise network");
		hybrid_automaton_pair* p =
				const_cast<hybrid_automaton_pair*> (my_aut.get());
		lid = p->get_or_add_location_id(l1, l2);
		return lid;
	}
	;

private:
	hybrid_automaton_pair::const_ptr my_aut;
	pw_network_node* child1;
	pw_network_node* child2;
};
*/

pairwise_hybrid_automaton_network::pairwise_hybrid_automaton_network() :
	my_impl(hybrid_automaton::ptr()) //, my_root(NULL)
{
}

pairwise_hybrid_automaton_network::~pairwise_hybrid_automaton_network() {
	//delete my_root;
}

pairwise_hybrid_automaton_network::ptr pairwise_hybrid_automaton_network::get_ptr() {
	return boost::static_pointer_cast<pairwise_hybrid_automaton_network>(
			hybrid_automaton::get_ptr());
}

pairwise_hybrid_automaton_network::const_ptr pairwise_hybrid_automaton_network::get_const_ptr() const {
	return boost::static_pointer_cast<const pairwise_hybrid_automaton_network>(
			hybrid_automaton::get_const_ptr());
}

void pairwise_hybrid_automaton_network::set_impl(
		const hybrid_automaton::ptr& aut) {
	my_impl = aut;
	my_impl_automata[aut->get_id()] = aut;
}

void pairwise_hybrid_automaton_network::add_base_automata(automaton_id aut_id) {
	hybrid_automaton::ptr aut = hybrid_automaton_cache::get_automaton(aut_id);
	if (!aut) {
		throw basic_exception(
				"pairwise_hybrid_automaton_network::add_base_automata: Automaton "
						+ to_string(aut_id) + " not found in automaton cache.");
	}
	// test if this is a base automaton
	hybrid_automaton_network::ptr net_aut = boost::dynamic_pointer_cast<
			hybrid_automaton_network>(aut);
	if (net_aut) {
		// add the base automata from net_aut
		automaton_id_set child_ids = net_aut->get_automata();
		for (automaton_id_set::const_iterator it = child_ids.begin(); it
				!= child_ids.end(); ++it) {
			add_base_automata(*it);
		}
	} else {
		my_base_automata[aut->get_id()] = aut;
	}
}

hybrid_automaton_network::ptr pairwise_hybrid_automaton_network::compute_or_assign_composition(
		hybrid_automaton::ptr aut) {
//	std::cout << "network id " << get_id() << " adding implementation " << aut->get_id();
	my_automata.insert(aut->get_id());
	// make sure the automaton is in the cache, otherwise the rest won't work
	hybrid_automaton_cache::add_automaton(aut);
	add_base_automata(aut->get_id());

	//pw_network_nonpair_node* added_leaf_node = new pw_network_nonpair_node(aut);
	if (!my_impl) {
		// It is the only automaton in the composition. The composition is identical to aut.
		set_impl(aut);
		//delete my_root;
		//my_root = added_leaf_node;
	} else {
		// create new automaton consisting of my_impl and aut
		// we want a new impl automaton, not a new network,
		hybrid_automaton::ptr new_impl(my_impl->create());
		hybrid_automaton_pair::ptr comp_aut = hybrid_automaton_pair::ptr(
				new hybrid_automaton_pair());
		comp_aut->init(new_impl, my_impl, aut);
		set_impl(comp_aut);
	}
	return get_ptr();
}

automaton_id_set pairwise_hybrid_automaton_network::get_automata() const {
	return my_automata;
}

location_id_set pairwise_hybrid_automaton_network::get_locations(
		const location_constraint_set& lcons) const {
	assert(my_impl);

	location_id_set locset;
	location_constraint_set fullcons; // empty set

	automaton_id this_id=get_id();
	automaton_id impl_id=my_impl->get_id();

	// shortcut if it's a single constraint that defines the top implementation
	location_constraint_set cons=lcons;
	cons.map(this_id,impl_id);

	locset = my_impl->get_locations(cons);

	//std::cout << lcons << " ----> " << locset << std::endl;
	return locset;
}

/*
location_id_set pairwise_hybrid_automaton_network::get_locations(
		const location_constraint_set& lcons) const {
	assert(my_impl);

	location_id_set locset;
	location_constraint_set fullcons; // empty set

	std::cout << "my id " << get_id() << " imp id " << my_impl->get_id()
			<< std::endl;

	// shortcut if it's a single constraint that defines the top implementation
	bool treated = false;
	if (lcons.size() == 1) {
		const automaton_id& aut_id = lcons.begin()->first;
		const location_constraint* con = &(lcons.begin()->second);
		if (aut_id == get_id()) {
			location_constraint_set lcons2;
			lcons2.add_constraint(my_impl->get_id(), *con);
			locset = my_impl->get_locations(lcons2);
			treated = true;
		} else if (aut_id == my_impl->get_id()) {
			locset = my_impl->get_locations(lcons);
			treated = true;
		}
	}
	if (!treated) {
		// @todo treat the case of mixed constraints over implementations and
		// originals
		//
		// first check if an implementation is referenced
		// for every implementation, remove the children from the list
		// of automata that needs expansion
		//
		// For now: assure it's just originals
		//assert(set_contains(my_automata_base,lcons.get_automata()));

		if (!lcons.is_empty())
			add_location_constraint(my_base_automata.begin(), my_base_automata.end(),
					lcons, fullcons, locset);
	}
	//std::cout << lcons << " ----> " << locset << std::endl;
	return locset;
}
*/

bool pairwise_hybrid_automaton_network::canonicalize_location_constraint(
		const automaton_id& aut_id, const location_constraint& con,
		location_constraint_set& lcons) const {
	return my_impl->canonicalize_location_constraint(my_impl->get_id(), con, lcons);
}

//location_id pairwise_hybrid_automaton_network::get_location_id(
//		const location_constraint_set& lcons) const {
//	//return my_root->get_location_id(lcons);
//}

// ----------------------------------------------------------------------------

hybrid_automaton* pairwise_hybrid_automaton_network::create() const {
	assert(my_impl);
	return my_impl->create();
}

pairwise_hybrid_automaton_network* pairwise_hybrid_automaton_network::create_network() const {
	return new pairwise_hybrid_automaton_network();
}

const symbolic_state_collection_ptr& pairwise_hybrid_automaton_network::get_initial_states() const {
	return my_impl->get_initial_states();
}

void pairwise_hybrid_automaton_network::set_initial_states(
		const symbolic_state_collection_ptr& sstate_set) {
	my_impl->set_initial_states(sstate_set);
}

std::pair<hybrid_automaton_pair::transition_const_iterator,
		hybrid_automaton_pair::transition_const_iterator> pairwise_hybrid_automaton_network::get_outgoing_transitions(
		location_id l, label_id a) const {
	return my_impl->get_outgoing_transitions(l, a);
}

transition_ptr pairwise_hybrid_automaton_network::get_transition(
		const transition_id& id) const {
	return my_impl->get_transition(id);
}

location_ptr pairwise_hybrid_automaton_network::get_location(
		const location_id& id) const {
	return my_impl->get_location(id);
}

location_id pairwise_hybrid_automaton_network::get_location_id(
		std::string loc_name) const {
	return my_impl->get_location_id(loc_name);
}

std::pair<hybrid_automaton_pair::location_const_iterator,
		hybrid_automaton_pair::location_const_iterator> pairwise_hybrid_automaton_network::get_locations() const {
	return my_impl->get_locations();
}

transition_id pairwise_hybrid_automaton_network::add_transition(
		const transition_ptr& trans, bool check_emptiness) {
	return my_impl->add_transition(trans, check_emptiness);
}

location_id pairwise_hybrid_automaton_network::add_location(
		const location_ptr& loc) {
	return my_impl->add_location(loc);
}

const label_id_set& pairwise_hybrid_automaton_network::get_labels() const {
	return my_impl->get_labels();
}

void pairwise_hybrid_automaton_network::add_label(const label_id& lab) {
	my_impl->add_label(lab);
}

void pairwise_hybrid_automaton_network::accept(hybrid_automaton_visitor& v) {
	std::string impl_name = my_impl->get_name();
	my_impl->set_name(get_name());

	try {
		my_impl->accept(v);
	} catch (std::exception& e) {
		if (my_automata.size()>1)
			throw basic_exception("Failure inside automaton network "+dot_context::context_free_name(get_name())+".",e);
		else {
			// If size 1, the implementation is an actual automaton (not an intermediary
			// so let's report it's name, too.
			throw basic_exception("Failure inside component "+dot_context::context_free_name(impl_name)+" of network "+dot_context::context_free_name(get_name())+".",e);
		}
	}

	my_impl->set_name(impl_name);
	// note : the visitor is passed on to all automata in
	// the network via the hybrid_automaton_pair accept method
}


/** Output as a stream of characters.
 */
void pairwise_hybrid_automaton_network::print(std::ostream& os) const {
	if (my_impl) {
		std::string impl_name = my_impl->get_name();
		my_impl->set_name(get_name());
		my_impl->print(os);
		my_impl->set_name(impl_name);
	}
}

const std::string& pairwise_hybrid_automaton_network::get_name() const {
	//return my_impl->get_name();
	return hybrid_automaton::get_name();
}

void pairwise_hybrid_automaton_network::set_name(std::string s) {
	//my_impl->set_name(s+"_pairwise_hybrid_automaton_network_implementation");
	hybrid_automaton::set_name(s);
}

void pairwise_hybrid_automaton_network::add_variable(const variable_id& vid,
		bool is_input, bool is_constant) {
	my_impl->add_variable(vid, is_input, is_constant);
}
void pairwise_hybrid_automaton_network::add_variables(
		const variable_id_set& vars, const variable_id_set& inp_vars,
		const variable_id_set& const_vars) {
	my_impl->add_variables(vars, inp_vars, const_vars);
}
const variable_id_set& pairwise_hybrid_automaton_network::get_variable_ids() const {
	return my_impl->get_variable_ids();
}
const variable_id_set& pairwise_hybrid_automaton_network::get_input_variables() const {
	return my_impl->get_input_variables();
}
const variable_id_set& pairwise_hybrid_automaton_network::get_const_variables() const {
	return my_impl->get_const_variables();
}
}

