#include "core/pwl/PLWL_sstate_collection_stl_list.h"
//#include "symbolic_state_collection_stl_list.cpp"

#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection_stl_list.h"

namespace hybrid_automata {

PLWL_sstate_collection_stl_list::PLWL_sstate_collection_stl_list() {
	my_merge_passed_with_new = false;
	passed_list = symbolic_state_collection_stl_list::ptr(
			new symbolic_state_collection_stl_list);
	waiting_list = symbolic_state_collection_stl_list::ptr(
			new symbolic_state_collection_stl_list);
}

PLWL_sstate_collection_stl_list* PLWL_sstate_collection_stl_list::create() const {
	return new PLWL_sstate_collection_stl_list();
}

PLWL_sstate_collection_stl_list* PLWL_sstate_collection_stl_list::clone() const {

	PLWL_sstate_collection_stl_list* new_list = new PLWL_sstate_collection_stl_list(*this);

	new_list->passed_list = symbolic_state_collection_stl_list::ptr(
			passed_list->clone());
	new_list->waiting_list = symbolic_state_collection_stl_list::ptr(
			waiting_list->clone());

	return new_list;
}

symbolic_state_collection::ptr PLWL_sstate_collection_stl_list::get_passed_list() const {
	return passed_list;
}

symbolic_state_collection::ptr PLWL_sstate_collection_stl_list::get_waiting_list() const {
	//return symbolic_state_collection::ptr(const_cast<symbolic_state_collection_stl_list*>(&waiting_list));
	return waiting_list;
}

PLWL_sstate_collection_stl_list::info PLWL_sstate_collection_stl_list::get_info() const {
	info ret_info = { passed_list->size(), waiting_list->size() };
	return ret_info;
}

void PLWL_sstate_collection_stl_list::set_passed_list(
		const symbolic_state_collection::ptr& new_pl) {
	passed_list->copy(new_pl);
}

void PLWL_sstate_collection_stl_list::set_waiting_list(
		const symbolic_state_collection::ptr& new_wl) {
	waiting_list->copy(new_wl);
}

bool PLWL_sstate_collection_stl_list::is_waiting_list_empty() const {
	return waiting_list->is_empty();
}

symbolic_state::ptr PLWL_sstate_collection_stl_list::pick() {
	return waiting_list->erase_front();
}

void PLWL_sstate_collection_stl_list::add_passed(symbolic_state::ptr sstate,
		const hybrid_automaton_ptr H) {

	// Add state to passed list and retrieve the part that's new
	if (my_merge_passed_with_new) {
//std::cout << "adding to passed with merging:" << sstate << std::endl;
		passed_list->add_and_return_new_with_merging(sstate); // This returns in sstate only the new states
	} else {
//		passed_list->add(sstate); // this does no checking whatsoever
//		std::cout << "adding to passed without merging:" << sstate << std::endl;
		passed_list->add_and_return_new(sstate); // This returns in sstate only the new states*
	}
//	std::cout << "passed list:" << passed_list << std::endl;
}

void PLWL_sstate_collection_stl_list::add_waiting(symbolic_state::ptr sstate,
		const hybrid_automaton_ptr H) {
	// compute the states that are new
	passed_list->compute_new(sstate);
//std::cout << "not on passed list:" << sstate << std::endl;

	// Add remainder (the part that's new) to waiting list
//	if (my_merge_passed_with_new) {
//		/* The same kind of add_new() as above can be applied to waiting
//		 * list too. The overhead of add_new may be more than having
//		 * redundant states in WL although.
//		 */
//		waiting_list->compute_new(sstate);
//		waiting_list->add(sstate);
//	} else
//		waiting_list->add_and_return_new_with_merging(sstate); // This returns in sstate only the new states

	if (my_merge_passed_with_new) {
//std::cout << "adding to waiting with merging:" << sstate << std::endl;
		waiting_list->add_and_return_new_with_merging(sstate); // This returns in sstate only the new states
	} else {
//std::cout << "adding to waiting without merging:" << sstate << std::endl;
		waiting_list->add_and_return_new(sstate);
	}
//	std::cout << "waiting list:" << waiting_list << std::endl;
}

void PLWL_sstate_collection_stl_list::clear_passed_list() {
	passed_list->clear();
}

void PLWL_sstate_collection_stl_list::clear_waiting_list() {
	waiting_list->clear();
}

/**
 * \brief Destructor
 */
PLWL_sstate_collection_stl_list::~PLWL_sstate_collection_stl_list() {
}

void PLWL_sstate_collection_stl_list::set_option(std::string opt) {
	if (opt == "MERGEP_ON") {
		my_merge_passed_with_new = true;
	} else if (opt == "MERGEP_OFF") {
		my_merge_passed_with_new = false;
	} else {
		throw std::runtime_error(
				"PLWL_sstate_collection_stl_list: unknown option " + opt + ".");
	}
}

}
