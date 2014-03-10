#ifndef EXPLICIT_POST_H_
#define EXPLICIT_POST_H_

virtual void add_post_discrete(
		symbolic_state_collection::ptr result_set,
		symbolic_state::ptr sstate) const;
virtual void add_post_time_infinity(
		symbolic_state_collection::ptr result_set,
		symbolic_state::ptr sstate) const;
virtual void add_post_time_bounded(
		symbolic_state_collection::ptr result_set,
		symbolic_state::ptr sstate, double delta) const;

void explicit_automaton::add_post_discrete(
		symbolic_state_collection::ptr result_set,
		symbolic_state::ptr sstate) const {
	discrete_set_stl_set_ptr d= static_pointer_cast<discrete_set_stl_set,
			discrete_set>(sstate->get_discrete_set());
	continuous_set_ptr c=sstate->get_continuous_set();
	// iterate through the locations in the set
	for (discrete_set_stl_set::const_iterator loc_it=d->begin(); loc_it
			!=d->end(); ++loc_it) {
		// iterate through all the labels
		for (label_id_set::const_iterator lab_it=my_labels.begin(); lab_it
				!=my_labels.end(); ++lab_it) {
			// iterate through all the outgoing transitions of the location
			std::pair<transition_const_iterator,transition_const_iterator>
					transp= get_outgoing_transitions(*loc_it, *lab_it);
			for (transition_const_iterator trans_it=transp.first; trans_it
					!=transp.second; ++trans_it) {
				// apply the post of the transition to the continuous states
				continuous_set_ptr c_res=(*trans_it)->post_discrete(c);

				// create a discrete_set consisting of the target location
				discrete_set_stl_set_ptr d_res=
						discrete_set_stl_set_ptr(new discrete_set_stl_set);
				d_res->add((*trans_it)->get_target());
				// add both to the result
				symbolic_state::ptr s=symbolic_state::ptr(new symbolic_state(d_res,c_res));
				result_set->add(s);
			}
		}
	}
}
void explicit_automaton::add_post_time_infinity(
		symbolic_state_collection::ptr result_set,
		symbolic_state::ptr sstate) const {

	discrete_set_stl_set_ptr d= static_pointer_cast<discrete_set_stl_set,
			discrete_set>(sstate->get_discrete_set());
	continuous_set_ptr c=sstate->get_continuous_set();
	// iterate through the locations in the set
	for (discrete_set_stl_set::const_iterator loc_it=d->begin(); loc_it
			!=d->end(); ++loc_it) {
		// apply the post of the location to the continuous states
		location::ptr lptr=get_location(*loc_it);
		continuous_set_ptr c_res=lptr->post_time_infinity(c);

		// create a discrete_set consisting of the current location
		discrete_set_stl_set_ptr d_res=
				discrete_set_stl_set_ptr(new discrete_set_stl_set);
		d_res->add(*loc_it);

		// add both to the result
		symbolic_state::ptr s=symbolic_state::ptr(new symbolic_state(d_res,c_res));
		result_set->add(s);
	}
}
void explicit_automaton::add_post_time_bounded(
		symbolic_state_collection::ptr result_set,
		symbolic_state::ptr sstate, double delta) const {

	discrete_set_stl_set_ptr d= static_pointer_cast<discrete_set_stl_set,
			discrete_set>(sstate->get_discrete_set());
	continuous_set_ptr c=sstate->get_continuous_set();
	// iterate through the locations in the set
	for (discrete_set_stl_set::const_iterator loc_it=d->begin(); loc_it
			!=d->end(); ++loc_it) {
		// apply the post of the location to the continuous states
		location::ptr lptr=get_location(*loc_it);
		continuous_set_ptr c_res=lptr->post_time_bounded(c, delta);

		// create a discrete_set consisting of the current location
		discrete_set_stl_set_ptr d_res=
				discrete_set_stl_set_ptr(new discrete_set_stl_set);
		d_res->add(*loc_it);

		// add both to the result
		symbolic_state::ptr s=symbolic_state::ptr(new symbolic_state(d_res,c_res));
		result_set->add(s);
	}
}

#endif /*EXPLICIT_POST_H_*/
