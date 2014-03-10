#ifndef REACHABILITY_SCENARIO_H_
#define REACHABILITY_SCENARIO_H_

#include "boost/shared_ptr.hpp"
#include "application/options.h"
#include "core/hybrid_automata/adapt_automaton_visitor.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
class hybrid_automaton_network;
typedef boost::shared_ptr<hybrid_automaton_network>
		hybrid_automaton_network_ptr;
typedef boost::shared_ptr<const hybrid_automaton_network>
		hybrid_automaton_network_const_ptr;
class passed_and_waiting_list;
typedef boost::shared_ptr<passed_and_waiting_list> passed_and_waiting_list_ptr;
typedef boost::shared_ptr<const passed_and_waiting_list>
		passed_and_waiting_list_const_ptr;
class symbolic_state;
typedef boost::shared_ptr<symbolic_state> symbolic_state_ptr;
typedef boost::shared_ptr<const symbolic_state> symbolic_state_const_ptr;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection>
		symbolic_state_collection_ptr;
typedef boost::shared_ptr<const symbolic_state_collection>
		symbolic_state_collection_const_ptr;
class post_operator;
typedef boost::shared_ptr<post_operator> post_operator_ptr;
typedef boost::shared_ptr<const post_operator> post_operator_const_ptr;
class discrete_post;
typedef boost::shared_ptr<discrete_post> discrete_post_ptr;
typedef boost::shared_ptr<const discrete_post> discrete_post_const_ptr;
class continuous_post;
typedef boost::shared_ptr<continuous_post> continuous_post_ptr;
typedef boost::shared_ptr<const continuous_post>
		continuous_post_const_ptr;
//class adapt_automaton_visitor;
//typedef boost::shared_ptr<adapt_automaton_visitor> adapt_automaton_visitor_ptr;
}

namespace hybrid_automata {

/** A reachability_scenario provides continuous and discrete post operators,
 * and a PWL list.
 *
 * All semantic information (how a scenario does things) is passed by objects
 * (post operators etc.). This is so that new scenarios can be built by
 * combining copying old ones and mixing/combining their elements.  */

class reachability_scenario {
public:
	typedef adapt_automaton_visitor::coefficient_type coefficient_type;
	typedef void (*option_handler)(reachability_scenario*,options::options_processor::variables_map& vmap);

	reachability_scenario();
	hybrid_automaton_ptr create_hybrid_automaton() const;
	hybrid_automaton_network_ptr create_hybrid_automaton_network() const;
	passed_and_waiting_list_ptr create_passed_and_waiting_list() const;
	symbolic_state_collection_ptr create_symbolic_state_collection() const;
	continuous_post_ptr create_continuous_post_operator() const;
	discrete_post_ptr create_discrete_post_operator() const;
	adapt_automaton_visitor_ptr create_adapt_automaton_visitor() const;
	coefficient_type get_bool_type() const;
	coefficient_type get_number_type() const;
	continuous_post_ptr get_continuous_post_operator();
	discrete_post_ptr get_discrete_post_operator();
	passed_and_waiting_list_ptr get_passed_and_waiting_list();
	symbolic_state_collection_ptr get_symbolic_state_collection();
	hybrid_automaton_ptr get_hybrid_automaton();
	hybrid_automaton_network_ptr get_hybrid_automaton_network();
	adapt_automaton_visitor_ptr get_adapt_automaton_visitor();
	void set_continuous_post_operator(continuous_post_ptr p);
	void set_discrete_post_operator(discrete_post_ptr p);
	void set_passed_and_waiting_list(passed_and_waiting_list_ptr p);
	void set_symbolic_state_collection(symbolic_state_collection_ptr p);
	void set_hybrid_automaton(hybrid_automaton_ptr p);
	void set_hybrid_automaton_network(hybrid_automaton_network_ptr p);
	void set_adapt_automaton_visitor(adapt_automaton_visitor_ptr p);

	void apply_options(options::options_processor::variables_map& vmap);
	void set_option_handler(option_handler f);

	/** Set time horizon to th.
	 *
	 * A negative value denotes unlimited time. */
	void set_time_horizon(double th);
	double get_time_horizon() const;

	/** Set sampling time to ts.
	 *
	 * A negative value denotes no sampling. */
	void set_sampling_time(double ts);
	double get_sampling_time() const;

	/** Set number of iterations of the algorithm to k.
	 *
	 * A negative value denotes no limit. */
	void set_iter_max(int k);
	int get_iter_max() const;

	/** Set relative error at intersection operations.
	 *
	 * A negative value denotes no guaranteed bound. */
	void set_intersection_error(double e);
	double get_intersection_error() const;

	/**
	 * Sets the intersection algorithm type to be carried out.
	 * @param algo_type Could be lb_search or lb_search_bb
	 */
	void set_intersection_type(std::string algo_type);
	std::string get_intersection_type() const;
	/**
	 * Sets the minima bracketing algorithm.
	 * @param minbrak_type
	 */
	void set_minbrak_type(std::string minbrak_type);
	std::string get_minbrak_type() const;

	/**
	 * Sets the sfm interval size to chull before intersecting with
	 * the guard constraint.
	 *
	 * @param split_size
	 */
	void set_split_size(size_t split_size);
	size_t get_split_size() const ;

	/** Returns the name of the scenario. */
	std::string get_name() const;

	/** Sets the name of the scenario. */
	void set_name(const std::string& n);

private:
	continuous_post_ptr my_continuous_post_operator;
	discrete_post_ptr my_discrete_post_operator;
	passed_and_waiting_list_ptr my_pwl;
	symbolic_state_collection_ptr my_symbolic_state_collection;
	hybrid_automaton_ptr my_hybrid_automaton;
	hybrid_automaton_network_ptr my_hybrid_automaton_network;
	adapt_automaton_visitor_ptr my_adapt_automaton_visitor;
	option_handler my_option_handler;

	double my_time_horizon;
	double my_sampling_time;
	int my_iter_max;
	double my_intersection_error;
	std::string my_intersection_type;
	std::string my_minbrak_type;
	size_t my_split_size;
	std::string my_name;
};

}

#endif /*REACHABILITY_SCENARIO_H_*/
