/*
 * monitoring_options.cpp
 *
 *  Created on: Sep 28, 2009
 *      Author: frehse
 */

#include "core/analysis_algorithms/monitoring/monitoring_options.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/analysis_algorithms/monitoring/ha_monitor.h"

#include <iostream>
#include <iterator>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "application/options.h"

namespace options {

using namespace std;

namespace {
/**
 * Tokenizes the option and puts the first option in the method_name parameter.
 * Also puts the parsed tokens in the 'tokens' object and tok_cnt has the number of
 * tokens in the 'tokens' object.
 *
 * @Returns An iterator to the next token.
 **/

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

tokenizer tokenize_events(std::string command, std::vector<std::string>& toks, unsigned int& tok_cnt) {
	unsigned int cnt=0;
	//Tokenise the command
	boost::char_separator<char> sep(" ,\r");
	tokenizer tokens(command, sep);
	for(tokenizer::iterator tok_counter = tokens.begin();tok_counter!=tokens.end();tok_counter++){
		toks.push_back(*tok_counter);
		cnt++;
	}
	tok_cnt = cnt;
	return tokens;
}
}

/*
 * UPPAAL TRON Interface
 * @author Rajarshi
 */

//void process_monitoring_options(const string& command, hybrid_automata::reachability_scenario& scen) {
void process_monitoring_options(hybrid_automata::hybrid_automaton::ptr aut_net,
		hybrid_automata::symbolic_state_collection::ptr ini_states,
		const hybrid_automata::reachability_scenario& scen,
		std::vector<io::output_formatter*> of) {

	unsigned int tok_cnt;
	hybrid_automata::ha_monitor<double> my_ha_monitor(aut_net, ini_states, scen, of);

	hybrid_automata::symbolic_state_collection::const_ptr previous_states(ini_states->clone());
	//logger::copyfmt_to(std::cout);

	bool ret_type = false; // assume the set is not empty

	do{
		std::cout << "Monitoring command options are:\n"
				<< "[reach_unobs]\n[take_transitions]\n[filter_states]\n[get_transitions]\n[show_current_states]\n[show_previous_states]\n[trace]\nType EXIT to quit Monitoring console.\n>>";
		//std::cout << "\nEnter the next monitoring command:\n";

		std::string command;
		cin >> command;

		if (command.compare("reach_unobs")==0) {
	/*
			if(tok_cnt < 2){
				throw std::runtime_error("MONITORING: Illegal arguments passed\n");
			}
			unsigned int t1, t2;
			vector<char> unobs_events;
			t1 = boost::lexical_cast<unsigned int>(*tok_iter);
			tok_iter++;
			t2 = boost::lexical_cast<unsigned int>(*tok_iter);
			tok_iter++;
			while (tok_iter != tokens.end()) {
				unobs_events.push_back(boost::lexical_cast<char>(*tok_iter));
				tok_iter++;
			}
	*/
			std::cout << "Enter Time Instant t1, t2, followed by the observable events, which are NOT executed:\n";
			double t1,t2;
			//std::string events;
			string events;
			std::vector<std::string> unobs_events;

			std::cin >> t1 >> t2;
			std::getline(cin,events);
			tokenize_events(events,unobs_events,tok_cnt);
/*
			while(e != '\n'){
				unobs_events.push_back(e);
				e = std::cin.get();
			};
*/

			//call the required function to compute reach_unobs(t1,t2,unobs_events)

/*
			cout << "t1=" << t1 << ",t2=" << t2 << "\n";
			cout << "events are:\n";

			for (vector<std::string>::iterator it = unobs_events.begin(); it
					!= unobs_events.end(); it++)
				std::cout << *it << std::endl;
*/

			/* Call the function that computes the reach_unobs */
			previous_states = hybrid_automata::symbolic_state_collection::const_ptr(my_ha_monitor.get_current_states()->clone());
			ret_type = my_ha_monitor.reach_unobs(t1, t2, unobs_events);
			(ret_type == true) ? std::cout << "\nResulting Reach Set Empty\n": std::cout << "\nResulting Reach Set Not Empty" << " (" << my_ha_monitor.get_current_states()->size() << " symbolic states)\n";

		} else if (command.compare("take_transitions")==0) {
			std::cout << "Enter the Transition Event:\n";
/*
			unsigned int t1, t2;
			std::cin >> t1 >> t2 >> event;
*/
			std::string event;
			std::cin.clear();
			std::cin >> event;

			previous_states = hybrid_automata::symbolic_state_collection::const_ptr(my_ha_monitor.get_current_states()->clone());
			ret_type = my_ha_monitor.take_transitions(event);
			(ret_type == true) ? std::cout << "\nResulting Reach Set Empty\n": std::cout << "\nResulting Reach Set Not Empty" << " (" << my_ha_monitor.get_current_states()->size() << " symbolic states)\n";
			//std::cout << " Number of arguments in take_transitions:" << tok_cnt;

/*
			else if (tok_cnt == 3){
				t1 = boost::lexical_cast<unsigned int>(*tok_iter);
				tok_iter++;
				t2 = boost::lexical_cast<unsigned int>(*tok_iter);
				event = boost::lexical_cast<char>(*++tok_iter);
				//call the function take_transition with the arguments
				cout << "t1=" << t1 << ",t2=" << t2 << "\n";
				cout << "Transition event is:" << event << endl;
			}
*/
			//else{ throw std::runtime_error("MONITORING: illegal arguments passed\n");}
		} else if (command.compare("filter_states")==0) {
			double low,high;
			std::string var;
			std::cout << "Enter the filter variable followed by the lower and upper bound" << std::endl;
			std::cin >> var >> low >> high;
			std::cout << "Reach Set Filtering in progress...\n";

			previous_states = hybrid_automata::symbolic_state_collection::const_ptr(my_ha_monitor.get_current_states()->clone());
			ret_type = my_ha_monitor.filter_states(var,low,high);
			(ret_type == true) ? std::cout << "\nResulting Reach Set Empty\n": std::cout << "\nResulting Reach Set Not Empty" << " (" << my_ha_monitor.get_current_states()->size() << " symbolic states)\n";

		} else if (command.compare("get_transitions")==0){
			std::cout << "Enter the Time Horizon:" << std::endl;
			double T;
			std::cin >> T;
			typedef hybrid_automata::ha_monitor<double>::label_intv_map interval_map_type;
			interval_map_type my_intv_map = my_ha_monitor.get_transitions(T);
			/* Debug Purpose */
			//std::cout << "get_transitions map:\n";
			for(interval_map_type::const_iterator it = my_intv_map.begin(); it!=my_intv_map.end();++it){
				std::cout <<"Label:" << it->first << std::endl;
				std::cout <<"Interval:" << it->second << std::endl << std::endl;
			}

		} else if (command.compare("show_current_states")==0){
			std::cout << "Current set of states (interval bounds):" << std::endl;
			io::INTV_formatter of(std::cout);
			of.set_context(aut_net->get_name());
			//LOGGER_OS(LOW, "monitoring", "Forbidden states are not reachable.");
			of.output(*my_ha_monitor.get_current_states());
			//std::cout << *my_ha_monitor.get_current_states();
		} else if (command.compare("trace") == 0) {
			std::string param;
			std::cin >> param;
			if (param=="on") {
				my_ha_monitor.output_states = true;
//			} else if (param=="time") {
//				// add time to the output variables
//				for (unsigned int i = 0; i < of.size(); ++i) {
//					variable_id_list vars = of[i]->get_output_variables();
//					if (find(vars.begin(), vars.end(),
//							my_ha_monitor.get_timer_id()) == vars.end()) {
//						vars.push_front(my_ha_monitor.get_timer_id());
//						of[i]->set_output_variables(vars);
//					}
//				}
//				my_ha_monitor.output_states = true;
			} else if (param=="off") {
				my_ha_monitor.output_states = false;
			} else {
				std::cout << "Monitoring: Unknown trace option: " << param;
			}
		} else if (command.compare("EXIT") == 0) {
			return;
		} else if (command.compare("show_previous_states")==0){
			// will be treated below; include here so that it's not considered an unknown option
		}
		else {
			std::cout << "Monitoring: Unknown Option:\n";
		}

		// check whether the current states are empty
		if (command.compare("show_previous_states")==0 || ret_type) {
			std::cout << "Last non-empty set of states (interval bounds):" << std::endl;
			io::INTV_formatter of(std::cout);
			of.set_context(aut_net->get_name());
			of.output(*previous_states);
		}
	} while(!ret_type); // end of do while
}
;

void process_testing_options(const string& command) {
/*
	cout << "Testing Interface" << endl;
	cout << "Options are: " << command << endl;
	vector<char> events;
	string method_name;
	unsigned int tok_cnt;
	tokenizer tokens = tokenize_command(command, method_name,tok_cnt);
	tokenizer::iterator tok_iter = tokens.begin();
	tok_iter++;

	if (method_name == "get_intervals") {
		vector<char> events;
		while (tok_iter != tokens.end()) {
			events.push_back(boost::lexical_cast<char>(*tok_iter));
			tok_iter++;
		}
		for (vector<char>::iterator it = events.begin(); it != events.end(); it++)
			cout << *it << endl;
		// call required function to compute the intervals in which these events will be enabled.
	} else if (method_name == "take_trans_var") {
		unsigned int t1, t2;
		t1 = boost::lexical_cast<unsigned int>(*tok_iter);
		tok_iter++;
		t2 = boost::lexical_cast<unsigned int>(*tok_iter);
		char event = boost::lexical_cast<char>(*++tok_iter);
		char variable = boost::lexical_cast<char>(*++tok_iter);
		;
		// call required function with the parsed tokens
		cout << "t1=" << t1 << ",t2=" << t2 << "\n";
		cout << "Transition event is:" << event << endl;
		cout << "Variable is: " << variable << "\n";
		//call required function with parsed arguments.
	} else {
		runtime_error("Unknown Option");
	}
*/

}
;

void add_monitoring_options() {
	/*@author Rajarshi
	 * Options for Talking with UPPAAL TRON
	 */
//	options::options_processor::config.add_options()("MONITORING",
//			boost::program_options::value<string>()->notifier(
//					&process_monitoring_options),
//			"Commands to Interface for Monitoring");
	options::options_processor::config.add_options()("MONITORING",
			"Enter the Commands for Monitoring:\n");

	options::options_processor::config.add_options()("TESTING",
			boost::program_options::value<string>()->notifier(
					&process_testing_options), "Commands for Testing");
}

bool check_monitoring_options(options::options_processor::variables_map& vmap) {
	return true;
}

bool apply_monitoring_options(options::options_processor::variables_map& vmap) {
	return true;
}

bool use_monitoring(options::options_processor::variables_map& vmap) {
	 if (vmap.count("MONITORING")) {
		 return true;
	 } else return false;
}

void monitor(hybrid_automata::hybrid_automaton::ptr aut_net,
		hybrid_automata::symbolic_state_collection::ptr ini_states,
		options::options_processor::variables_map& vmap, std::vector<io::output_formatter*> of) {
	const hybrid_automata::reachability_scenario& scen =
			hybrid_automata::scenario_chooser::get_scenario();
	process_monitoring_options(aut_net, ini_states, scen, of);

}

}
