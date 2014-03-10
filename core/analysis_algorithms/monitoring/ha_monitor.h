
#ifndef _HA_MONITOR_
#define _HA_MONITOR_

#include "core/hybrid_automata/explicit_automaton.h"
#include "core/scenarios/scenario_chooser.h"
#include "math/vdom/lin_expression.h"
#include "utility/logger_stopwatch.h"
#include "math/numeric/interval.h"

namespace hybrid_automata {

template <typename scalar_type>
class ha_monitor {
public:

	typedef math::numeric::interval<scalar_type> time_interval;
	typedef std::multimap<std::string, time_interval> label_intv_map;
	typedef std::set<time_interval> interval_set;

	/** Create a monitoring setup for the hybrid automaton ha_ptr, which has initial states ini_states, and
	 * using the scenario scen.
	 */
	ha_monitor(hybrid_automaton::ptr ha_ptr,
			hybrid_automata::symbolic_state_collection::ptr ini_states,
			hybrid_automata::reachability_scenario scen,
			std::vector<io::output_formatter*> output_formatters);

	/** destructor*/
	virtual ~ha_monitor(){};

	/**
	 * Computes all the reachable states between time instants t1 and t2, considering delay and all possible
	 * jumps except the ones with labels in the events list.
	 *
	 * @param t1 Start time instant
	 * @param t2 End time instant
	 * @param events List of events
	 * @return true if the computed current state set is empty, false otherwise
	 */
	bool reach_unobs(const scalar_type& t1, const scalar_type& t2, const std::vector<std::string>& events);

	/**
	 * Triggers the discrete transitions with the passed label on the current state set, without applying any delay.
	 * @return true if the computed current state set is empty, false otherwise.
	 */

	bool take_transitions(const std::string& label);

	/**
	 * Filter the current_states_set with the constraint low <= variable <= high,
	 * where variable is passed with is variable id as parameter.
	 *
	 * @param var_id The identifier of the variable of the constraint
	 * @param low Lower bound on the variable
	 * @param high Upper bpund on the variable
	 * @return true if the computed current state set is empty, false otherwise.
	 */
	bool filter_states(const std::string& variable, const scalar_type& low, const scalar_type& high);

	/**
	 * Computes the list of available labeled transition(s) with the time interval
	 * when they are enabled. The list is stored as a std::multimap since there
	 * could be more than 1 interval associated with a label for a given time horizon.
	 *
	 * @param time_horizon Upper bound on the time.
	 */
	label_intv_map get_transitions(const scalar_type& time_horizon) const;
	/**
	 * @todo Dont understand what to do here!
	 * @param t_s
	 * @param t_e
	 * @param label
	 * @param variable
	 * @return
	 */
	interval_set get_variable_bounds(const scalar_type& t_s, const scalar_type& t_e, const std::string& label,const std::string& variable) const;

	/** Returns the current states
	 *
	 */
	symbolic_state_collection::const_ptr get_current_states() const;
	/**
	 * Clusters a symbolic_states_collection accosrding to the user clustering
	 * options options. The sfm_set on the symbolic states are converted to polytope
	 * collections.
	 *
	 * @param compose_ha_ptr
	 * @param postc_states
	 * @return Clustered symbolic state collection with sfm changed to poly collection
	 */
	symbolic_state_collection::ptr cluster_sym_states(hybrid_automata::hybrid_automaton::ptr compose_ha_ptr,
			symbolic_state_collection::ptr postc_states);

	/** Returns the name of the monitor automaton */
	static std::string monitor_name();

	/** Bool for activating output of reachable states */
	bool output_states;

	/** Get the id of the variable used as global timer */
	static const variable_id& get_timer_id();

	/** Get the id of the variable used as global timer */
	static const std::string& get_timer_name();

private:
	/** forbid use of default constructor  */
	ha_monitor();

	hybrid_automata::reachability_scenario my_scen;
	hybrid_automata::symbolic_state_collection::ptr current_states_set; /*the reachable state set at a time instant */
	hybrid_automata::hybrid_automaton::ptr my_ha_ptr; /* The HA to monitor */
	std::vector<io::output_formatter*> of; // vector of output formatters
};


} // end of hybrid automata namespace

#include "ha_monitor.hpp"

#endif

