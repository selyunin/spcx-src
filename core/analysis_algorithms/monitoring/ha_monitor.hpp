#include "ha_monitor.h"
//#include "core/analysis_algorithms/reachability_algorithm.h"

#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/continuous_dynamics/continuous_dynamics.h"
#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"
#include "core/analysis_algorithms/reachability_algorithm.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/hybrid_automaton_pair.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/explicit_transition.h"
#include "core/discrete/singleton_set.h"
#include "core/discrete/discrete_set.h"
#include "core/hybrid_automata/location_constraint_set.h"
#include "math/vdom/lin_constraint.h"
#include "core/scenarios/reachability_scenario.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection_stl_list.h"
#include "core/scenarios/support_fun_scenario.h"
#include "io/common_output/output_options.h"
#include "core/continuous/continuous_dynamics/continuous_dynamics.h"
#include "core/hybrid_automata/hybrid_automaton_network.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "io/GEN_format/GEN_formatter.h"
#include "core/hybrid_automata/hybrid_automaton_utility.h"
#include "core/post_operators/discrete_post.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"
#include "math/scalar_types/scalar_with_infinity.h"


namespace hybrid_automata {

template<typename scalar_type>
std::string ha_monitor<scalar_type>::monitor_name() {
	return "MONITOR";
}

template<typename scalar_type>
const variable_id& ha_monitor<scalar_type>::get_timer_id() {
	/* static initialization, see also
	 * http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.15
	 * http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.16
	 */
	static variable_id my_id(variable::get_or_add_variable_id(get_timer_name()));
	return my_id;
}

template<typename scalar_type>
const std::string& ha_monitor<scalar_type>::get_timer_name() {
	/* static initialization, see also
	 * http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.15
	 * http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.16
	 */
	static std::string my_timer_name(monitor_name()+".TIMER");
	return my_timer_name;
}

template<typename scalar_type>
ha_monitor<scalar_type>::ha_monitor(hybrid_automaton::ptr ha_ptr,
		hybrid_automata::symbolic_state_collection::ptr ini_states,
		hybrid_automata::reachability_scenario scen,
		std::vector<io::output_formatter*> output_formatters) {
	using namespace math;
	using namespace continuous;
	my_ha_ptr = ha_ptr;
	my_scen = scen;
	// create the constraint TIMER==0
	lin_expression<scalar_type> v;
	v.set_coeff_with_id(get_timer_id(), scalar_type(1));
	v.set_inh_coeff(scalar_type(0));
	lin_constraint<scalar_type> inv_cons;
	inv_cons = lin_constraint<scalar_type> (v, EQ);
	typename constr_polyhedron<scalar_type>::ptr init_poly(
			new constr_polyhedron<scalar_type> ());
	init_poly->add_constraint(inv_cons);
	// intersect continuous sets to add TIMER==0
	symbolic_state_collection_stl_list::ptr start_states(new symbolic_state_collection_stl_list());
	// = my_scen.create_symbolic_state_collection();
	for (symbolic_state_collection::const_iterator it =
			ini_states->begin(); it != ini_states->end(); ++it) {
		discrete::discrete_set_ptr dset = (*it)->get_discrete_set();
		continuous_set_ptr cset = (*it)->get_continuous_set();
		continuous_set_ptr newset = compute_intersection(cset, init_poly);
		start_states->add_and_return_new_with_merging(symbolic_state::ptr(new symbolic_state(dset, newset)));
	}
	hybrid_automata::canonicalize_location_constraints(*start_states);
	current_states_set = start_states;

	// by default, don't output states
	output_states = false;
	of = output_formatters;
}

template <typename scalar_type>
bool ha_monitor<scalar_type>::reach_unobs(const scalar_type& t1, const scalar_type& t2, const std::vector<std::string>& unobs_events){
	using namespace math;
	using namespace continuous;

	assert(t2 >= t1);

	/*step 1: make a monitor automata */

	hybrid_automata::hybrid_automaton::ptr monitor_ptr = hybrid_automata::hybrid_automaton::ptr(new hybrid_automata::explicit_automaton(monitor_name()));
	hybrid_automata::hybrid_automaton_cache::add_automaton(monitor_ptr);
	lin_constraint<scalar_type> inv_cons;

	typename constr_polyhedron<scalar_type>::ptr init_poly, inv_poly_start, inv_poly_stop, guard_poly, dyn_poly;

	monitor_ptr->add_variable(get_timer_id());
	monitor_ptr->add_label(monitor_name()+".TAU");


	inv_poly_start = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	inv_poly_stop = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());

	// add constraint: TIMER >= 0
	lin_expression<scalar_type> v;
	v.set_coeff_with_id(get_timer_id(), scalar_type(1));
	v.set_inh_coeff(scalar_type(0));
	inv_cons = lin_constraint<scalar_type> (v,GE);

	inv_poly_start->add_constraint(inv_cons);

	// add another constraint: TIMER <= t1
	v.set_coeff_with_id(get_timer_id(),scalar_type(1));
	v.set_inh_coeff(-t1);
	inv_cons = lin_constraint<scalar_type> (v,LE);
	inv_poly_start->add_constraint(inv_cons);

	// add another constraint: TIMER >= t1
	v.set_coeff_with_id(get_timer_id(),scalar_type(1));
	v.set_inh_coeff(-t1);
	inv_cons = lin_constraint<scalar_type> (v,GE);
	inv_poly_stop->add_constraint(inv_cons);

	// add another constraint: TIMER <= t2
	v.set_coeff_with_id(get_timer_id(),scalar_type(1));
	v.set_inh_coeff(-t2);
	inv_cons = lin_constraint<scalar_type> (v,LE);
	inv_poly_stop->add_constraint(inv_cons);

	//continuous_set_ptr inv = inv_poly;

	// create the guard poly (TIMER >= t1)
	v.set_coeff_with_id(get_timer_id(),scalar_type(1));
	v.set_inh_coeff(-t1);
	inv_cons = lin_constraint<scalar_type>(v,GE);
	guard_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	guard_poly->add_constraint(inv_cons);


	//continuous_set_ptr deriv = dyn_poly;
	math::vdom_matrix<scalar_type> zero_mat(v.domain(),v.domain());
	math::vdom_vector<scalar_type> b(v.domain(),scalar_type(1));

	math::affine_map<scalar_type> dyn_map(zero_mat,b);

	continuous_dynamics::ptr dyn = continuous_dynamics::ptr(new ode_affine_dynamics<scalar_type>(dyn_map));

	location::ptr loc = location::ptr(new location("START", time_constraints(inv_poly_start, dyn)));
	monitor_ptr->add_location(loc);

	// Add another location

	//------------------

	dyn = continuous_dynamics::ptr(new ode_affine_dynamics<scalar_type>(dyn_map));
	loc = location::ptr(new location("STOP", time_constraints(inv_poly_stop, dyn)));
	monitor_ptr->add_location(loc);

	// add a transition

	math::affine_map<scalar_type> map(v.domain(),scalar_type(0)); // Make an identity affine map.
	continuous_set_transform::const_ptr transf = continuous_set_transform::const_ptr
			(new reset_affine_transform<scalar_type>(map));

	// get the ids of source and target locations

	location_id sloc = monitor_ptr->get_location_id("START");
	location_id tloc = monitor_ptr->get_location_id("STOP");

	transition::ptr t = transition::ptr(new explicit_transition(sloc, monitor_name()+".TAU", jump_constraints(guard_poly,transf), tloc));
	monitor_ptr->add_transition(t);


	symbolic_state_collection::ptr start_states = my_scen.create_symbolic_state_collection();

	// let's create some initial states
//	v.set_coeff_with_id(get_timer_id(),scalar_type(1));
//	v.set_inh_coeff(scalar_type(0));
//	inv_cons = lin_constraint<scalar_type> (v,EQ);
//	init_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
//	init_poly->add_constraint(inv_cons);

	// the initial states are the monitor in location START;
	//
	init_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	location_constraint_set lcs = location_constraint_set(monitor_ptr->get_id(),monitor_ptr->get_location_id("START"));
	discrete::discrete_set::ptr ini_d = discrete::discrete_set::ptr(new discrete::singleton_set(lcs));
	symbolic_state::ptr ini_ss = symbolic_state::ptr(new symbolic_state(ini_d, init_poly));
	start_states->add(ini_ss);

	//monitor_ptr->set_initial_states(start_states);
	/* Add Labels */
	for(std::vector<std::string>::const_iterator it = unobs_events.begin();it!=unobs_events.end(); it++)
		monitor_ptr->add_label(*it);

	//debug
	//std::cout << "Monitor Automata is:" << monitor_ptr << std::endl;
/*
	std::cout << "Monitor labels are:";
	label_id_set lid_set = monitor_ptr->get_labels();

	for(label_id_set::iterator it = lid_set.begin(); it != lid_set.end();it++)
		std::cout << *it << std::endl;
*/

	/* Monitor automata created */

	/* reset the initial set of the system automat to the current_states_set.*/
	//my_ha_ptr->set_initial_states(current_states_set);
	/* Make a composition of the monitor automata and my_ha*/

	//my_ha_ptr->set_initial_states(current_states_set); // Ensures that the system moves on from the current state set.

	hybrid_automaton_pair::ptr compose_ptr = hybrid_automaton_pair::ptr(new hybrid_automaton_pair());
	compose_ptr->set_name("compose");
	//my_ha_ptr->set_initial_states();
	hybrid_automaton::ptr composed_ha_ptr =  hybrid_automata::hybrid_automaton::ptr(my_ha_ptr->create());
	hybrid_automata::hybrid_automaton_cache::add_automaton(composed_ha_ptr);
	composed_ha_ptr->set_name("composed_ha");

	//Initialise the start states of my_ha_ptr with current_states-set
	start_states->intersection_assign(current_states_set);
	//my_ha_ptr->set_initial_states(current_states_set);

	compose_ptr->init(composed_ha_ptr,my_ha_ptr,monitor_ptr);
	// Note: compose_ptr should now have taken the place of composed_ha_ptr,
	// i.e., taken the name and id that composed_ha_ptr had previously.

/*
	hybrid_automaton_cache::clear();
	hybrid_automaton_network::ptr aut_net = scen.create_hybrid_automaton_network();
	hybrid_automaton_cache::add_automaton(aut_net);
	aut_net->set_name("sys");

	hybrid_automaton_cache::add_automaton(my_ha_ptr);
	hybrid_automaton_cache::add_automaton(monitor_ptr);

	aut_net = aut_net->compute_or_assign_composition(my_ha_ptr);
	aut_net = aut_net->compute_or_assign_composition(monitor_ptr);
*/

	/* Change the symbolic states to contain poly collection and not sfms*/
	//std::cout << "states before clustering =" << states << std::endl;
	start_states = cluster_sym_states(compose_ptr, start_states);
	//std::cout << "states after clustering = " << states << std::endl;

	symbolic_state_collection::ptr states, bad_states;
	reachability_algorithm alg(my_scen);
//	std::cout << "starting from states:" << start_states;
	states = alg.reach(compose_ptr,start_states);
	//std::cout << "Composed system:" << compose_ptr;

	hybrid_automata::canonicalize_location_constraints(*states);

	/* For printing */
	if (output_states) {
		LOGGERSW(LOW, "reach_unobs", "Output of reachable states");
		for (unsigned int i = 0; i < of.size(); ++i) {
			of[i]->output(*states);
		}
	}

/*
	std::cout << "system automata:" << my_ha_ptr <<std::endl;
	std::cout << "system automata initial states:" << my_ha_ptr->get_initial_states() << std::endl;
	std::cout << "monitor_automata:" << monitor_ptr << std::endl;
	std::cout << "Monitor automata initial states:" << monitor_ptr->get_initial_states() << std::endl;
	std::cout << "The composed automata is:\n" << compose_ptr << std::endl;
*/

	// get discrete post op from scenario
	discrete_post_ptr my_post_d;
	my_post_d = my_scen.create_discrete_post_operator();


	/*Check for emptiness of the reach set of the composed automata for time horizon of t2 */

	/* Create the Bad states to compute Reach(t1,t2) */
	lcs = location_constraint_set(monitor_ptr->get_id(),monitor_ptr->get_location_id("STOP"));
	ini_d = discrete::discrete_set::ptr(new discrete::singleton_set(lcs));
	init_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	ini_ss = symbolic_state::ptr(new symbolic_state(ini_d, init_poly));
	bad_states = my_scen.create_symbolic_state_collection();
	bad_states->add(ini_ss);

	if (bad_states) {
			states->intersection_assign(bad_states);
	}

	hybrid_automata::canonicalize_location_constraints(*states);

	{
		LOGGERSW(MEDIUM, "Monitoring::take_transitions", "Projecting and merging states");
		/*
		 * Existentially Quantify the monitor states
		 */
		variable_id_set monitor_ids;
		monitor_ids.insert(get_timer_id());
		//std::cout << "Variables Ids of Composed HA:" << compose_ptr->get_variable_ids() << std::endl;
		//std::cout << "Monitor Id=" << get_timer_id() << std::endl;
		symbolic_state_collection_stl_list::ptr merged_states(
				new symbolic_state_collection_stl_list());
		for (symbolic_state_collection::const_iterator it = states->begin(); it
				!= states->end(); ++it) {
			// quantify away the monitor automaton location (STOP)
			discrete::discrete_set_ptr dset = (*it)->get_discrete_set();
			dset->existentially_quantify(monitor_ptr->get_id());

			continuous_set_ptr cset = (*it)->get_continuous_set();
			// cset->existentially_quantify_variables(monitor_ids);
			//std::cout << "cset: " << cset << std::endl;
			merged_states->add_and_return_new_with_merging(
					symbolic_state::ptr(new symbolic_state(dset, cset)));
		}
		//std::cout << "current states :" << current_states_set;
		current_states_set = merged_states;
	}

/*
		std::ofstream my_file;
		my_file.open("/tmp/reach_unobs.txt");
		io::GEN_formatter formatter(my_file,0.0);
		formatter.output(*states);
		my_file.close();

		std::system("graph -TX -C -B /tmp/reach_unobs.txt");
*/

	//my_ha_ptr->set_initial_states(current_states_set);

	/* For printing */
	if (output_states) {
		LOGGERSW(LOW, "reach_unobs", "Output of reachable states");
		for (unsigned int i = 0; i < of.size(); ++i) {
			of[i]->output(*current_states_set);
		}
	}

	hybrid_automata::hybrid_automaton_cache::remove_automaton(monitor_ptr->get_id());
	hybrid_automata::hybrid_automaton_cache::remove_automaton(composed_ha_ptr->get_id());

	return math::definitely(current_states_set->is_empty());

}

template <typename scalar_type>
bool ha_monitor<scalar_type>::take_transitions(const std::string& event){
	using namespace math;
	using namespace continuous;

	/*step 1: make a monitor automata */

	hybrid_automata::hybrid_automaton::ptr monitor_ptr = hybrid_automata::hybrid_automaton::ptr(new hybrid_automata::explicit_automaton("monitor"));
	hybrid_automata::hybrid_automaton_cache::add_automaton(monitor_ptr);
	lin_constraint<scalar_type> inv_cons;

	typename constr_polyhedron<scalar_type>::ptr init_poly, inv_poly_start, inv_poly_stop, guard_poly, dyn_poly;

	monitor_ptr->add_variable(get_timer_id());
	positional_vdomain dom(monitor_ptr->get_variable_ids());

	inv_poly_start = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	math::vdom_matrix<scalar_type> zero_mat(dom,dom);
	math::vdom_vector<scalar_type> b(dom,scalar_type(1));
	math::affine_map<scalar_type> dyn_map(zero_mat,b);
	continuous_dynamics::ptr dyn = continuous_dynamics::ptr(new ode_affine_dynamics<scalar_type>(dyn_map));

	location::ptr loc = location::ptr(new location("START", time_constraints(inv_poly_start, dyn)));
	monitor_ptr->add_location(loc);

	/*step 2: compose system with monitor */
	hybrid_automaton_pair::ptr compose_ptr = hybrid_automaton_pair::ptr(new hybrid_automaton_pair());
	compose_ptr->set_name("compose");
	//my_ha_ptr->set_initial_states();
	hybrid_automaton::ptr composed_ha_ptr =  hybrid_automata::hybrid_automaton::ptr(my_ha_ptr->create());
	hybrid_automata::hybrid_automaton_cache::add_automaton(composed_ha_ptr);
	composed_ha_ptr->set_name("composed_ha");

	compose_ptr->init(composed_ha_ptr,my_ha_ptr,monitor_ptr);

/*  removed so things work with the phaver scenario also
 *
	typename support_function::sfm_cont_set<scalar_type>::const_ptr sfm_set_ptr;
	for(symbolic_state_collection::const_iterator it = current_states_set->begin();
			it!=current_states_set->end(); ++it){

		sfm_set_ptr =
				boost::dynamic_pointer_cast<const support_function::sfm_cont_set<scalar_type> >((*it)->get_continuous_set());
		if (!(sfm_set_ptr)) {
			std::cout << "MONITORING: take_transitions: Continuous Set is not an SFM in a sym_state of the current_states_set.\n "
					"May be because no postc operation done before call to take_transitions.\n";
			return false;
		}
	}
*/

	std::cout << "Computing Image of Transitions with Label " << event << " ....." << std::endl;
	symbolic_state_collection::ptr interm_states,result_states;

	discrete_post_ptr postd_ptr = my_scen.get_discrete_post_operator();

	result_states = my_scen.create_symbolic_state_collection();

	label_id event_label_id = named_label::get_or_add_label_id(event);

	{
		LOGGERSWNC(MEDIUM, "Monitoring::take_transitions", "Time to compute the reach state");
		unsigned int count = 0;
		for (symbolic_state_collection::const_iterator it =
				current_states_set->begin(); it != current_states_set->end(); ++it, count++) {
			postd_ptr->add_post_states(compose_ptr, result_states, result_states,
					event_label_id, (*it));
		}
		//std::cout << "number of sym_states inside current_states_set:" << count << std::endl;
	}


	/*	LOGGERSWNC(MEDIUM, "Monitoring::take_transitions", "Time to Plot");*/

	{
		LOGGERSW(MEDIUM, "Monitoring::take_transitions", "Projecting and merging states");
		hybrid_automata::canonicalize_location_constraints(*result_states);
		symbolic_state_collection_stl_list::ptr merged_states(
				new symbolic_state_collection_stl_list());
		for (symbolic_state_collection::const_iterator it =
				result_states->begin(); it != result_states->end(); ++it) {
			// quantify away the monitor automaton location (STOP)
			discrete::discrete_set_ptr dset = (*it)->get_discrete_set();
			dset->existentially_quantify(monitor_ptr->get_id());

			continuous_set_ptr cset = (*it)->get_continuous_set();
			// cset->existentially_quantify_variables(monitor_ids);
			//std::cout << "cset: " << cset << std::endl;
			merged_states->add_and_return_new_with_merging(
					symbolic_state::ptr(new symbolic_state(dset, cset)));
		}
		//std::cout << "current states :" << current_states_set;
		current_states_set = merged_states;
	}

	/* For printing */
	if (output_states) {
		LOGGERSW(LOW, "reach_unobs", "Output of reachable states");
		for (unsigned int i = 0; i < of.size(); ++i) {
			of[i]->output(*current_states_set);
		}
	}

	hybrid_automata::hybrid_automaton_cache::remove_automaton(monitor_ptr->get_id());
	hybrid_automata::hybrid_automaton_cache::remove_automaton(composed_ha_ptr->get_id());

	return (math::definitely(current_states_set->is_empty()));
}

template<typename scalar_type>
bool ha_monitor<scalar_type>::filter_states(const std::string& variable,
		const scalar_type& low,
		const scalar_type& high){
	using namespace math;
	using namespace continuous;

	// need to get variable using context dependent lookup
	std::string var_name =
			valuation_functions::variable_node_creator::context_lookup(variable,
					my_ha_ptr->get_name());
	variable_id var_id = variable::get_or_add_variable_id(var_name, 1);

	// Make constraint: variable <= high
	lin_expression<scalar_type> lin_exp;
	math::lin_constraint<scalar_type> filter_cons1, filter_cons2;

	lin_exp.set_coeff_with_id(var_id, scalar_type(1));
	lin_exp.set_inh_coeff(-scalar_type(high));
	filter_cons1 = math::lin_constraint<scalar_type> (lin_exp,LE);

	// Make constraint: -variable<= -low
	lin_exp.set_coeff_with_id(var_id, scalar_type(-1));
	lin_exp.set_inh_coeff(scalar_type(low));
	filter_cons2 = math::lin_constraint<scalar_type> (lin_exp,LE);

	// We need to filter out states that become empty.
	// So let's copy everything to a new list
	symbolic_state_collection_stl_list::ptr result_states(new symbolic_state_collection_stl_list());
	//result_states = my_scen.create_symbolic_state_collection();

	// Add the filtering constraint
	for(symbolic_state_collection::const_iterator it = current_states_set->begin();
			it!=current_states_set->end(); ++it){

		discrete::discrete_set::ptr dset = (*it)->get_discrete_set();
		continuous_set_ptr c_set_ptr = (*it)->get_continuous_set();

		bool res_empty = false;

		typename support_function::sfm_cont_set<scalar_type>::ptr sfm_set_ptr =
				boost::dynamic_pointer_cast<support_function::sfm_cont_set<scalar_type> >(c_set_ptr);
		if(sfm_set_ptr){
			sfm_set_ptr->intersection_with_constraint(filter_cons1);
			sfm_set_ptr->intersection_with_constraint(filter_cons2);
			// intersection is not carried out explicitly
			// need to check emptiness properly
			res_empty = math::definitely(sfm_set_ptr->is_empty_outer());
		}
		else{
			typename polyhedron<scalar_type>::ptr poly_set_ptr =
							boost::dynamic_pointer_cast<polyhedron<scalar_type> >(c_set_ptr);
			if(poly_set_ptr){
				poly_set_ptr->add_constraint(filter_cons1);
				poly_set_ptr->add_constraint(filter_cons2);
			}
			else{
				typename polyhedron_collection<scalar_type>::ptr poly_collection_ptr =
						boost::dynamic_pointer_cast<polyhedron_collection<scalar_type> >(c_set_ptr);
				if(poly_collection_ptr){
					poly_collection_ptr->add_constraint(filter_cons1);
					poly_collection_ptr->add_constraint(filter_cons2);
				}
				else{
					throw std::runtime_error("MONITORING:filter_states: Continuous set type not supported.");
				}
			}
			res_empty = math::definitely(c_set_ptr->is_empty());
		}

		if (!res_empty) {
//std::cout << "adding " << c_set_ptr << " to " << 	result_states << std::endl;
			symbolic_state::ptr new_state(new symbolic_state(dset, c_set_ptr));
			result_states->add_and_return_new_with_merging(new_state);
//std::cout << "got back " << new_state << std::endl;
		}
	}

	current_states_set = result_states;

	/* For printing */
	if (output_states) {
		LOGGERSW(LOW, "reach_unobs", "Output of reachable states");
		for (unsigned int i = 0; i < of.size(); ++i) {
			of[i]->output(*current_states_set);
		}
	}

	return (math::definitely(current_states_set->is_empty()));
}


template<typename scalar_type>
typename ha_monitor<scalar_type>::label_intv_map ha_monitor<scalar_type>::get_transitions(const scalar_type& T) const {
	using namespace math;
	using namespace continuous;

	label_intv_map my_map;
	std::list<transition_id> out_trans;
	const label_id_set& lab_set = my_ha_ptr->get_labels();
//	double T = double(time_horizon);
//	my_scen.set_time_horizon(time_horizon);

	/* Make a monitor automata */

	/*step 1: make a monitor automata */

	hybrid_automata::hybrid_automaton::ptr monitor_ptr = hybrid_automata::hybrid_automaton::ptr(new hybrid_automata::explicit_automaton("monitor"));
	hybrid_automata::hybrid_automaton_cache::add_automaton(monitor_ptr);
	lin_constraint<scalar_type> inv_cons;

	typename constr_polyhedron<scalar_type>::ptr inv_poly,dyn_poly;

	monitor_ptr->add_variable(get_timer_id());
	variable_id z_id = variable::get_or_add_variable_id("monitor.z");
	monitor_ptr->add_variable(z_id);

	monitor_ptr->add_label("monitor.tau");

	// add constraint: TIMER >= 0
	lin_expression<scalar_type> v;
	v.set_coeff_with_id(get_timer_id(), scalar_type(1));
	v.set_coeff_with_id(z_id, scalar_type(0));
	v.set_inh_coeff(scalar_type(0));
	inv_cons = lin_constraint<scalar_type> (v,GE);
	inv_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	inv_poly->add_constraint(inv_cons);

	// add another constraint: TIMER <=T
	v.set_coeff_with_id(get_timer_id(),scalar_type(1));
	v.set_coeff_with_id(z_id, scalar_type(0));
	v.set_inh_coeff(scalar_type(-1) * T);
	inv_cons = lin_constraint<scalar_type> (v,LE);
	inv_poly->add_constraint(inv_cons);

	continuous_set_ptr inv = inv_poly;

	// Create Dynamics
	math::vdom_matrix<scalar_type> zero_mat(v.domain(),v.domain());
	math::vdom_vector<scalar_type> b(v.domain(),scalar_type(1));
	b.set_coeff_with_id(z_id,scalar_type(0));

	math::affine_map<scalar_type> dyn_map(zero_mat,b); // timer'=1 and z'=0

	continuous_dynamics::ptr dyn = continuous_dynamics::ptr(new ode_affine_dynamics<scalar_type>(dyn_map));

	location::ptr loc = location::ptr(new location("START", time_constraints(inv, dyn)));
	monitor_ptr->add_location(loc);

	// Add Labels
	for(label_id_set::const_iterator it = lab_set.begin(); it!= lab_set.end(); ++it)
		monitor_ptr->add_label(*it);

	// Add a stop location for each system label and a transition from START to this location
	location_id sloc,tloc;
	std::string state_name;
	math::vdom_vector<scalar_type> b_stop;
	math::affine_map<scalar_type> map;
	continuous_set_transform::const_ptr transf;
	transition::ptr t;

	for(label_id_set::const_iterator it = lab_set.begin(); it != lab_set.end(); ++it){
		state_name = "STOP_" + hybrid_automata::named_label::get_name(*it);

		// TIMER <= T
		v.set_coeff_with_id(get_timer_id(),scalar_type(1));
		v.set_coeff_with_id(z_id,scalar_type(0));
		v.set_inh_coeff(scalar_type(-1)* T);
		inv_cons = lin_constraint<scalar_type> (v,LE);
		inv_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
		inv_poly->add_constraint(inv_cons);

		// add constraint: TIMER >= 0
		v.set_coeff_with_id(get_timer_id(), scalar_type(1));
		v.set_coeff_with_id(z_id, scalar_type(0));
		v.set_inh_coeff(scalar_type(0));
		inv_cons = lin_constraint<scalar_type> (v,GE);
		inv_poly->add_constraint(inv_cons);

		// z = 0
		v.set_coeff_with_id(get_timer_id(),scalar_type(0));
		v.set_coeff_with_id(z_id,scalar_type(1));
		v.set_inh_coeff(scalar_type(0));
		inv_cons = lin_constraint<scalar_type> (v,EQ);
		inv_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
		inv_poly->add_constraint(inv_cons);

		inv = inv_poly;

		b_stop = math::vdom_vector<scalar_type>(v.domain(),scalar_type(1));

		dyn_map = math::affine_map<scalar_type>(zero_mat,b_stop);

		dyn = continuous_dynamics::ptr(new ode_affine_dynamics<scalar_type>(dyn_map));
		loc = location::ptr(new location(state_name, time_constraints(inv, dyn)));
		monitor_ptr->add_location(loc);

		// add a transition

		map = math::affine_map<scalar_type>(v.domain(),v.domain()); // Make an identity affine map.
		transf = continuous_set_transform::const_ptr
				(new reset_affine_transform<scalar_type>(map));

		// get the ids of source and target locations

		sloc = monitor_ptr->get_location_id("START");
		tloc = monitor_ptr->get_location_id(state_name);

		t = transition::ptr(new explicit_transition(sloc,hybrid_automata::named_label::get_name(*it), jump_constraints(transf), tloc));
		monitor_ptr->add_transition(t);

		// add another transition
		map = math::affine_map<scalar_type>(v.domain(),v.domain()); // Make an identity affine map.
		transf = continuous_set_transform::const_ptr
				(new reset_affine_transform<scalar_type>(map));

		sloc = monitor_ptr->get_location_id(state_name);
		tloc = monitor_ptr->get_location_id("START");

		t = transition::ptr(new explicit_transition(sloc,"monitor.tau", jump_constraints(transf), tloc));
		monitor_ptr->add_transition(t);

	}



	symbolic_state_collection::ptr start_states = my_scen.create_symbolic_state_collection();

	// let's create some initial states
	inv_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());

	// constraint TIMER==0
//	v.set_coeff_with_id(get_timer_id(),scalar_type(1));
//	v.set_coeff_with_id(z_id,scalar_type(0));
//	v.set_inh_coeff(0);
//	inv_cons = lin_constraint<scalar_type> (v,EQ);
//	inv_poly->add_constraint(inv_cons);

	// constraint z==0
	v.set_coeff_with_id(get_timer_id(),scalar_type(0));
	v.set_coeff_with_id(z_id,scalar_type(1));
	v.set_inh_coeff(0);
	inv_cons = lin_constraint<scalar_type> (v,EQ);
	inv_poly->add_constraint(inv_cons);


	continuous_set::ptr ini_c  = inv_poly;


	location_constraint_set lcs = location_constraint_set(monitor_ptr->get_id(),monitor_ptr->get_location_id("START"));
	discrete::discrete_set::ptr ini_d = discrete::discrete_set::ptr(new discrete::singleton_set(lcs));
	symbolic_state::ptr ini_ss = symbolic_state::ptr(new symbolic_state(ini_d, ini_c));
	start_states->add(ini_ss);

	monitor_ptr->set_initial_states(start_states);
	//DEBUG
	//std::cout << "Monitor Automata:\n" << monitor_ptr  << std::endl;
	//std::cout << "System Automata:\n" << my_ha_ptr << std::endl;

	hybrid_automaton_pair::ptr compose_ptr = hybrid_automaton_pair::ptr(new hybrid_automaton_pair());
	compose_ptr->set_name("compose");
	hybrid_automaton::ptr composed_ha_ptr =  hybrid_automata::hybrid_automaton::ptr(my_ha_ptr->create());
	hybrid_automata::hybrid_automaton_cache::add_automaton(composed_ha_ptr);
	composed_ha_ptr->set_name("composed_ha");

	compose_ptr->init(composed_ha_ptr,my_ha_ptr,monitor_ptr);
	// Note: compose_ptr should now have taken the place of composed_ha_ptr,
	// i.e., taken the name and id that composed_ha_ptr had previously.

	symbolic_state_collection::ptr states,bad_states,alt_states;
	reachability_algorithm alg(my_scen);

	// the start states are the monitor states intersected with the current states
	start_states->intersection_assign(current_states_set);

	states = alg.reach(compose_ptr,start_states);


	/* For printing */
	if (output_states) {
		LOGGERSW(LOW, "get_transitions", "Output of reachable states");
		for (unsigned int i = 0; i < of.size(); ++i) {
			of[i]->output(*states);
		}
	}

	scalar_type t_s, t_e;
	bool is_empty, is_bounded;
	math::vdom_vector<scalar_type> l_s,l_e,support_vec;

	l_e = math::vdom_vector<scalar_type>(v.domain());
	l_e.set_coeff_with_id(get_timer_id(),scalar_type(1));
	l_s = math::vdom_vector<scalar_type>(v.domain());
	l_s.set_coeff_with_id(get_timer_id(),scalar_type(-1));

	continuous::continuous_set_ptr cont_set_ptr;
	time_interval time_intv;
	std::pair<std::string,time_interval> my_pair;

	my_map = label_intv_map();
	for(label_id_set::const_iterator it = lab_set.begin(); it != lab_set.end(); ++it){
		std::string label_name  = hybrid_automata::named_label::get_name(*it);
		LOGGER_OS(DEBUG7,"ha_monitor<scalar_type>::get_transitions") << "Computing transition times for label " << label_name;
		state_name = "STOP_" + label_name;
		lcs = location_constraint_set(monitor_ptr->get_id(),monitor_ptr->get_location_id(state_name));
		ini_d = discrete::discrete_set::ptr(new discrete::singleton_set(lcs));
		ini_c = continuous_set::ptr(inv_poly->create_universe());
		ini_ss = symbolic_state::ptr(new symbolic_state(ini_d, ini_c));
		bad_states = my_scen.create_symbolic_state_collection();
		bad_states->add(ini_ss);

		alt_states = my_scen.create_symbolic_state_collection();
		alt_states->copy(states);

		hybrid_automata::canonicalize_location_constraints(*alt_states);

		alt_states->intersection_assign(bad_states);

		for(symbolic_state_collection::const_iterator iter = alt_states->begin();
				iter!=alt_states->end(); ++iter){
			cont_set_ptr = (*iter)->get_continuous_set();
			typename continuous::support_function_provider::const_ptr sf_provider_ptr =
					boost::dynamic_pointer_cast<const continuous::support_function_provider>(cont_set_ptr);

			if (!sf_provider_ptr)
				std::runtime_error("Monitoring:get_transitions:state not a support function provider!\n");

			sf_provider_ptr->compute_support(l_s,t_s,support_vec,is_empty,is_bounded);
			if(!is_bounded){
				std::runtime_error("Monitoring:get_transitions:No lower bound found with timer variable!\n");
			}

			sf_provider_ptr->compute_support(l_e,t_e,support_vec,is_empty,is_bounded);
			if(!is_bounded){
				std::runtime_error("Monitoring:get_transitions:No upper bound found with timer variable!\n");
			}

			if(!is_bounded){
				time_intv = time_interval(scalar_with_infinity<scalar_type>(scalar_type(-1) * t_s));
			}
			else{
				time_intv = time_interval(scalar_with_infinity<scalar_type>(scalar_type(-1) * t_s),scalar_with_infinity<scalar_type>( t_e));
			}
			my_pair = std::pair<std::string, time_interval>(label_name,time_intv);
			my_map.insert(my_pair);
		}
	}

	/* To Visualize
	std::ofstream my_file;

	my_file.open("/tmp/reach_unobs.txt");
	io::GEN_formatter formatter(my_file,0.0);


	variable_id_list vid_list;
	vid_list.push_back(variable::get_or_add_variable_id("system.x"));
	vid_list.push_back(variable::get_or_add_variable_id("system.v"));


	formatter.set_output_variables(vid_list);
	formatter.output(*states);
	my_file.close();

	std::system("graph -TX -C -B /tmp/reach_unobs.txt"); */

	hybrid_automata::hybrid_automaton_cache::remove_automaton(monitor_ptr->get_id());
	hybrid_automata::hybrid_automaton_cache::remove_automaton(composed_ha_ptr->get_id());

	return my_map;

}
template<typename scalar_type>
typename ha_monitor<scalar_type>::interval_set ha_monitor<scalar_type>::get_variable_bounds(const scalar_type& t_s, const scalar_type& t_e,
		const std::string& label, const std::string& variable) const {

	interval_set my_set;
	return my_set;
}

template<typename scalar_type>
hybrid_automata::symbolic_state_collection::const_ptr ha_monitor<scalar_type>::get_current_states() const {
	return current_states_set;
}
/**
 * Clusters a collection of symbolic states obtained from post_c
 * computation according to the scenario clustering options.
 *
 * @param sym_states
 * @return a symbolic state
 */

template<typename scalar_type>
symbolic_state_collection::ptr ha_monitor<scalar_type>::cluster_sym_states(
		hybrid_automaton::ptr compose_ha_ptr, symbolic_state_collection::ptr sym_states){
	using namespace continuous;

	continuous_set::ptr cluster_states_ptr;
	continuous_set::const_ptr source_inv, target_inv;


	symbolic_state_collection_stl_list::ptr result(new symbolic_state_collection_stl_list());
	//= my_scen.create_symbolic_state_collection();

	const variable_id_set& v_ids = compose_ha_ptr->get_variable_ids();
	positional_vdomain vdom = positional_vdomain(v_ids);

	math::affine_map<scalar_type> map(vdom,scalar_type(0)); // Make an identity affine map.
	typename continuous_set_transform::const_ptr transf = typename continuous_set_transform::const_ptr(new reset_affine_transform<scalar_type>(map));

	typename constr_polyhedron<scalar_type>::ptr guard_poly,src_inv_poly, tgt_inv_poly;
	guard_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	guard_poly = typename constr_polyhedron<scalar_type>::ptr(guard_poly->create_universe());

	const jump_constraints& j = jump_constraints(guard_poly, transf);

	src_inv_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	tgt_inv_poly = typename constr_polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	src_inv_poly = typename constr_polyhedron<scalar_type>::ptr(src_inv_poly->create_universe());
	tgt_inv_poly = typename constr_polyhedron<scalar_type>::ptr(tgt_inv_poly->create_universe());

	discrete_post_ptr postd = my_scen.get_discrete_post_operator();
	continuous_set_collection cset_coll;
	for(symbolic_state_collection::const_iterator it = sym_states->begin();it!=sym_states->end();++it){
		discrete::discrete_set::ptr dset = (*it)->get_discrete_set();
		continuous_set::const_ptr cset = (*it)->get_continuous_set();

		// If it's an sfm we need to turn it into polyhedra
		if (boost::dynamic_pointer_cast<const support_function::sfm_cont_set<scalar_type> >(cset)) {
			cset_coll = postd->post(j, src_inv_poly, tgt_inv_poly, cset);

			// add the continuous sets to the list individually
			for (continuous::continuous_set_collection::iterator it =
					cset_coll.begin(); it != cset_coll.end(); ++it) {
				continuous::continuous_set::ptr cset = *it;
				if (cset && !math::definitely(cset->is_empty())) {
					symbolic_state::ptr new_sym_state = symbolic_state::ptr(
							new symbolic_state(dset, cset));
					//std::cout << "sym_state before adding to list:" << new_sym_state << std::endl;
					//result->add_and_return_new_with_merging(new_sym_state);
					// Note: not merging here because we'll do it later
					result->add(new_sym_state);
				}
			}
		} else {
			result->add(*it);
		}

	}
	//std::cout << "sym_states_coll size = " << result->size() << std::endl;
	//std::cout << "result sym states:" << result << std::endl;
	return result;

}

} // end of hybrid_automata namespace
