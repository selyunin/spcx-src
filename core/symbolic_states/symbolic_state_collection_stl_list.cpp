#include "core/symbolic_states/symbolic_state_collection_stl_list.h"

#include "core/continuous/continuous_set_operators.h"
#include "core/continuous/continuous_set.h"
#include "core/discrete/discrete_set.h"
#include "core/discrete/discrete_set_operators.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/hybrid_automata/adapt_automaton_visitor.h"

namespace hybrid_automata {

symbolic_state_collection_stl_list::~symbolic_state_collection_stl_list() {
}

symbolic_state_collection_stl_list::symbolic_state_collection_stl_list() {
}

symbolic_state_collection_stl_list* symbolic_state_collection_stl_list::create() const {
	return new symbolic_state_collection_stl_list;
}

symbolic_state_collection_stl_list* symbolic_state_collection_stl_list::clone() const {
	symbolic_state_collection_stl_list* stl_p = new symbolic_state_collection_stl_list();
	// Need to make a deep copy
	for (std::list<symbolic_state::ptr>::const_iterator it = mysymbolic_state_list.begin(); it
			!= mysymbolic_state_list.end(); ++it) {
		stl_p->add((*it)->clone());
	}
	return stl_p;
}

std::size_t symbolic_state_collection_stl_list::get_memory() const {
	//return symbolic_states_list.size()*sizeof(symbolic_state);
	// Need to see later.
	return 0;
}

void symbolic_state_collection_stl_list::add(const symbolic_state::ptr& sstate) {
	mysymbolic_state_list.push_back(sstate);
}

void symbolic_state_collection_stl_list::compute_new(const symbolic_state::ptr& sstate) {
	using namespace continuous;
	// Take away the states that are already in symbolic_state_list
	for (std::list<symbolic_state::ptr>::const_iterator it = mysymbolic_state_list.begin(); it
			!= mysymbolic_state_list.end(); ++it) {
		// If the discrete states of sstate are a subset of those in it, then we can substract from the continuous states of sstate
		if ((*it)->get_discrete_set()->contains(sstate->get_discrete_set())) {
			//			sstate->get_continuous_set()->cheap_difference_assign(it->get()->get_continuous_set());
			continuous_set_ptr new_part = compute_or_assign_cheap_difference(sstate->get_continuous_set(), (*it)->get_continuous_set());
			sstate->set_continuous_set(new_part);
			//std::cout << "difference:" << sstate << std::endl << std::flush;
			if (math::definitely(new_part->is_empty())) {
				return;
			}
		}
	}
}

void symbolic_state_collection_stl_list::add_and_return_new(const symbolic_state::ptr& sstate) {

	compute_new(sstate);
	// Add whatever is left unless it's empty
	if (!math::definitely(sstate->get_continuous_set()->is_empty()))
		add(sstate);
}

void symbolic_state_collection_stl_list::add_and_return_new_with_merging(
		const symbolic_state::ptr& sstate) {
	bool merged = false;

	// Note: If an accurate difference operator is used, merging should be performed before subtracting anything from sstate.
	// Here we use cheap difference (A\B returns either A or the empty set), so this seems less important

	bool is_empty = math::definitely(sstate->get_continuous_set()->is_empty());
	// Take away the states that are already in symbolic_state_list
	for (std::list<symbolic_state::ptr>::iterator it = mysymbolic_state_list.begin(); it
			!= mysymbolic_state_list.end() && !is_empty;) {
		// If the discrete states of sstate are a subset of those in it, then we can substract from the continuous states of sstate
		if ((*it)->get_discrete_set()->contains(sstate->get_discrete_set())) {
			//			sstate->get_continuous_set()->cheap_difference_assign(it->get()->get_continuous_set());
			sstate->set_continuous_set(compute_or_assign_cheap_difference(
					sstate->get_continuous_set(), (*it)->get_continuous_set()));
			//std::cout << "difference:" << sstate << std::endl << std::flush;
			is_empty = math::definitely(sstate->get_continuous_set()->is_empty());
		}

		bool to_delete = false;
		// Merging : if the new state contains an old one, replace the old by the new
		if (!is_empty && sstate->get_discrete_set()->contains((*it)->get_discrete_set())) {
			if (sstate->get_continuous_set()->contains((*it)->get_continuous_set())) {
				// replace the first one and delete all the other ones
				if (!merged) {
					(*it)->set_discrete_set(sstate->get_discrete_set());
					(*it)->set_continuous_set(sstate->get_continuous_set());
					merged = true;
				} else {
					to_delete = true;
				}
			}
		}

		if (!to_delete)
			++it;
		else
			it = mysymbolic_state_list.erase(it);
	}
	// Add whatever is left unless it's empty or already merged
	if (!merged && !is_empty)
		add(sstate);
}

bool symbolic_state_collection_stl_list::contains(const symbolic_state_collection::ptr& sstate_set) const {
	symbolic_state_collection::ptr new_sstate_set(sstate_set->clone());
	new_sstate_set->difference_assign(get_const_ptr());
	if (new_sstate_set->is_empty())
		return true;
	else
		return false;
}
void symbolic_state_collection_stl_list::cheap_difference_assign(const symbolic_state::ptr& sstate) {
	symbolic_state_collection::ptr new_list = symbolic_state_collection::ptr(
			new symbolic_state_collection_stl_list);

	for (std::list<symbolic_state::ptr>::iterator it1 = mysymbolic_state_list.begin(); it1
			!= mysymbolic_state_list.end(); ++it1) {
		discrete::discrete_set::ptr d1 = (*it1)->get_discrete_set();
		continuous::continuous_set::ptr c1 = (*it1)->get_continuous_set();

		discrete::discrete_set::ptr d2 = sstate->get_discrete_set();
		continuous::continuous_set::ptr c2 = sstate->get_continuous_set();

		discrete::discrete_set::ptr d = discrete::compute_difference(d1, d2);

		if (!d->is_empty()) {
			continuous::continuous_set::ptr c1_clone(c1->clone());
			symbolic_state::ptr s = symbolic_state::ptr(new symbolic_state(d, c1_clone));
			new_list->add(s);
		}

		d = discrete::compute_intersection(d1, d2);
		if (!d->is_empty()) {
			continuous::continuous_set::ptr c = compute_cheap_difference(c1, c2);

			if (!math::definitely(c->is_empty())) {
				symbolic_state::ptr s = symbolic_state::ptr(new symbolic_state(d, c));
				new_list->add(s);
			}
		}
	}
	swap(new_list);
}
void symbolic_state_collection_stl_list::cheap_difference_assign(
		const symbolic_state_collection::ptr& sstate_set) {
	for (symbolic_state_collection::const_iterator it2 = sstate_set->begin(); it2
			!= sstate_set->end(); ++it2) {
		cheap_difference_assign(*it2);
	}
}
void symbolic_state_collection_stl_list::difference_assign(const symbolic_state::const_ptr& sstate) {
	symbolic_state_collection::ptr new_list = symbolic_state_collection::ptr(
			new symbolic_state_collection_stl_list);
	new_list->clear();
	for (std::list<symbolic_state::ptr>::iterator it1 = mysymbolic_state_list.begin(); it1
			!= mysymbolic_state_list.end(); ++it1) {
		discrete::discrete_set::ptr d1 = (*it1)->get_discrete_set();
		continuous::continuous_set::ptr c1 = (*it1)->get_continuous_set();

		discrete::discrete_set::ptr d2 = sstate->get_discrete_set();
		continuous::continuous_set::ptr c2 = sstate->get_continuous_set();

		discrete::discrete_set::ptr d = discrete::compute_difference(d1, d2);

		if (!d->is_empty()) {
			symbolic_state::ptr s = symbolic_state::ptr(new symbolic_state(d, c1));
			new_list->add(s);
		}

		d = discrete::compute_intersection(d1, d2);
		if (!d->is_empty()) {
			continuous::continuous_set::ptr c = compute_cheap_difference(c1, c2);
			if (!math::definitely(c->is_empty())) {
				symbolic_state::ptr s = symbolic_state::ptr(new symbolic_state(d, c));
				new_list->add(s);
			}
		}
	}
	swap(new_list);
}

void symbolic_state_collection_stl_list::difference_assign(
		const symbolic_state_collection::const_ptr& sstate_set) {
	for (symbolic_state_collection::const_iterator it2 = sstate_set->begin(); it2
			!= sstate_set->end(); ++it2) {
		difference_assign(*it2);
	}
}

bool symbolic_state_collection_stl_list::accept(adapt_discrete_set_visitor& v) {
	bool success=true;
	for (std::list<symbolic_state::ptr>::const_iterator it = mysymbolic_state_list.begin(); it
			!= mysymbolic_state_list.end(); ++it) {
		v.reset();
		(*it)->get_discrete_set()->accept(v);
		if (v.get_success()) {
			(*it)->set_discrete_set(v.get_discrete_set());
		} else {
			success=false;
		}
	}
	return success;
}

bool symbolic_state_collection_stl_list::accept(adapt_continuous_set_visitor& v) {
	bool success=true;
	for (std::list<symbolic_state::ptr>::const_iterator it = mysymbolic_state_list.begin(); it
			!= mysymbolic_state_list.end(); ++it) {
		v.reset();
		(*it)->get_continuous_set()->accept(v);
		if (v.get_success()) {
			(*it)->set_continuous_set(v.get_continuous_set());
		} else {
			success=false;
		}
	}
	return success;
}

void symbolic_state_collection_stl_list::intersection_assign(
		const symbolic_state_collection::ptr& sstate_set) {
	assert(boost::dynamic_pointer_cast<symbolic_state_collection_stl_list>(sstate_set));

	std::list<symbolic_state_ptr> new_list;

	for (std::list<symbolic_state::ptr>::const_iterator it =
			mysymbolic_state_list.begin(); it != mysymbolic_state_list.end(); ++it) {
		for (symbolic_state_collection::const_iterator it1 =
				sstate_set->begin(); it1 != sstate_set->end(); ++it1) {
			discrete::discrete_set::ptr dset_ptr = compute_intersection(
					(*it)->get_discrete_set(),
					(*it1)->get_discrete_set());

			if (!dset_ptr->is_empty()) {
				continuous::continuous_set::ptr cset_ptr =
						compute_intersection((*it)->get_continuous_set(),
								(*it1)->get_continuous_set());

				if (!math::definitely(cset_ptr->is_empty())) {
					symbolic_state::ptr sstate_ptr =
							(symbolic_state::ptr) new symbolic_state;
					sstate_ptr->set_discrete_set(dset_ptr);
					sstate_ptr->set_continuous_set(cset_ptr);
					new_list.push_back(sstate_ptr);
				}
			}
		}
	}
	mysymbolic_state_list.swap(new_list);
}

symbolic_state::ptr symbolic_state_collection_stl_list::erase_front() {
	symbolic_state::ptr s = symbolic_state::ptr();
	if (!mysymbolic_state_list.empty()) {
		s = mysymbolic_state_list.front();
		mysymbolic_state_list.pop_front();
	}
	return s;
}

bool symbolic_state_collection_stl_list::is_empty() const {
	// it's empty unless there's a non-empty element in the list
	for (std::list<symbolic_state::ptr>::const_iterator it =
			mysymbolic_state_list.begin(); it != mysymbolic_state_list.end(); ++it) {
		if (!(*it)->is_empty())
			return false;
	}
	return true;
}

unsigned int symbolic_state_collection_stl_list::size() const {
	return mysymbolic_state_list.size();
}

void symbolic_state_collection_stl_list::swap(symbolic_state_collection::ptr sstate_set) {
	assert(boost::dynamic_pointer_cast<symbolic_state_collection_stl_list>(sstate_set));

	mysymbolic_state_list.swap(
			((symbolic_state_collection_stl_list *) sstate_set.get())->mysymbolic_state_list);
}

void symbolic_state_collection_stl_list::union_assign(
		const symbolic_state_collection::ptr& sstate_set) {
	/**
	 * Simply join the 2 lists.
	 */
	assert(boost::dynamic_pointer_cast<symbolic_state_collection_stl_list>(sstate_set));

	mysymbolic_state_list.splice(mysymbolic_state_list.end(),
			((symbolic_state_collection_stl_list *) sstate_set.get())->mysymbolic_state_list);

}

}
