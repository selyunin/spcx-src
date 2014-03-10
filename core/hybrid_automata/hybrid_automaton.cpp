#include "core/hybrid_automata/hybrid_automaton.h"

#include "core/discrete/discrete_set.h"
#include "core/discrete/singleton_set.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/transition.h"

namespace hybrid_automata {

hybrid_automaton::ptr hybrid_automaton::get_ptr() {
	hybrid_automaton::ptr p= boost::enable_shared_from_this<hybrid_automaton>::shared_from_this();
	return p;
}

hybrid_automaton::const_ptr hybrid_automaton::get_const_ptr() const {
	return boost::enable_shared_from_this<hybrid_automaton>::shared_from_this();
}

void hybrid_automaton::add_variables(const variable_id_set& vars,
		const variable_id_set& inp_vars,
		const variable_id_set& const_vars) {
	for (variable_id_set::const_iterator it = vars.begin(); it != vars.end(); ++it) {
		add_variable(*it, inp_vars.find(*it) != inp_vars.end(),
				const_vars.find(*it) != const_vars.end());
	}
}

variable_id_set hybrid_automaton::get_controlled_variables() const {
	variable_id_set contr_vars=get_variable_ids();
	const variable_id_set& input_vars=get_input_variables();
	set_difference_assign(contr_vars,input_vars);
	return contr_vars;
}

void hybrid_automaton::add_labels(const label_id_set& labs) {
	for (label_id_set::const_iterator it = labs.begin(); it != labs.end(); ++it) {
		add_label(*it);
	}
}

void hybrid_automaton::add_label(const std::string& lab_name) {
//	if (lab_name == named_label::silent_name())
//		throw basic_exception("silent label should not be added to alphabet");
	add_label(named_label::get_or_add_label_id(lab_name));
}

hybrid_automaton::hybrid_automaton() :
	my_id(get_new_id()), my_name("") {
}

hybrid_automaton::hybrid_automaton(std::string new_name) :
	my_id(get_new_id()) {
	set_name(new_name);
	//hybrid_automaton_cache::add_automaton(get_ptr());
}

hybrid_automaton::~hybrid_automaton() {
}

std::list<transition_id> hybrid_automaton::get_outgoing_transitions(
		const discrete::discrete_set::const_ptr d, label_id a) const {
	std::list<transition_id> out_list;
	std::pair<transition_const_iterator, transition_const_iterator> pr;
	for (discrete::discrete_set::const_iterator d_it = d->begin(); d_it != d->end(); ++d_it) {
		location_id_set locs = get_locations(*d_it);
		for (location_id_set::const_iterator l_it = locs.begin(); l_it != locs.end(); ++l_it) {
			pr = get_outgoing_transitions(*l_it, a);
			copy(pr.first, pr.second, std::back_inserter(out_list));
		}
	}
	return out_list;
}

std::list<discrete::discrete_set::ptr> hybrid_automaton::get_time_equiv_locations(
		const discrete::discrete_set::const_ptr& d,
		const continuous::continuous_set::const_ptr& cset) const {
	std::list<discrete::discrete_set::ptr> d_list;
	for (discrete::discrete_set::const_iterator d_it = d->begin(); d_it != d->end(); ++d_it) {
		location_id_set locs = get_locations(*d_it);
		for (location_id_set::const_iterator l_it = locs.begin(); l_it != locs.end(); ++l_it) {
			discrete::discrete_set::ptr new_d = discrete::discrete_set::ptr(
					new discrete::singleton_set(location_constraint_set(get_id(), *l_it)));
			d_list.push_back(new_d);
		}
	}
	return d_list;
}

label_id hybrid_automaton::get_label(const transition_id& id) const {
	return get_transition(id)->get_label();
}

const location_id& hybrid_automaton::get_source(const transition_id& id) const {
	return get_transition(id)->get_source();
}

const location_id& hybrid_automaton::get_target(const transition_id& id) const {
	return get_transition(id)->get_target();
}

const jump_constraints& hybrid_automaton::get_jump_constraints(const transition_id& id) const {
	return get_transition(id)->get_jump_constraints();
}

const std::string& hybrid_automaton::get_name() const {
	return my_name;
}

void hybrid_automaton::set_name(std::string s) {
	// Check whether the name is already in the cache
	//if (hybrid_automaton_cache::has_automaton(s))
	//	throw std::runtime_error("Automaton with name '"+s+"' already exists.");
	my_name=s;
}

const automaton_id& hybrid_automaton::get_id() const {
	return my_id;
}

automaton_id hybrid_automaton::get_new_id() {
	static automaton_id highest_id = 0;
	return ++highest_id;
}

}
