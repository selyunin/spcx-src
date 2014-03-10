/*
 * reachability_scenario.cpp
 *
 *  Created on: Jun 16, 2009
 *      Author: frehse
 */

#include "core/scenarios/reachability_scenario.h"

#include "core/post_operators/post_operator.h"
#include "core/post_operators/continuous_post.h"
#include "core/post_operators/discrete_post.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/hybrid_automaton_network.h"
#include "core/hybrid_automata/adapt_automaton_visitor.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/pwl/passed_and_waiting_list.h"

namespace hybrid_automata {

reachability_scenario::reachability_scenario() :
  my_continuous_post_operator(continuous_post_ptr()),
  my_discrete_post_operator(discrete_post_ptr()), my_pwl(
							 passed_and_waiting_list_ptr()),
  my_symbolic_state_collection(symbolic_state_collection_ptr()),
  my_hybrid_automaton(hybrid_automaton_ptr()),
  my_hybrid_automaton_network(hybrid_automaton_network_ptr()),
  my_option_handler(0),
  my_time_horizon(-1), my_sampling_time(1), my_iter_max(-1) {
}

hybrid_automaton_ptr reachability_scenario::create_hybrid_automaton() const {
	return hybrid_automaton_ptr(my_hybrid_automaton->create());
}

hybrid_automaton_network_ptr reachability_scenario::create_hybrid_automaton_network() const {
	return hybrid_automaton_network_ptr(my_hybrid_automaton_network->create_network());
}

passed_and_waiting_list_ptr reachability_scenario::create_passed_and_waiting_list() const {
	return passed_and_waiting_list_ptr(my_pwl->clone());
}

symbolic_state_collection_ptr reachability_scenario::create_symbolic_state_collection() const {
	return symbolic_state_collection_ptr(my_symbolic_state_collection->create());
}

continuous_post_ptr reachability_scenario::create_continuous_post_operator() const {
	continuous_post_ptr p=continuous_post_ptr(my_continuous_post_operator->clone());
	p->set_time_horizon(get_time_horizon());
	p->set_sampling_time(get_sampling_time());
	return p;
}

discrete_post_ptr reachability_scenario::create_discrete_post_operator() const {
	return discrete_post_ptr(my_discrete_post_operator->clone());
}

adapt_automaton_visitor_ptr reachability_scenario::create_adapt_automaton_visitor() const {
	return adapt_automaton_visitor_ptr(my_adapt_automaton_visitor->clone());
}

reachability_scenario::coefficient_type reachability_scenario::get_bool_type() const {
	assert(my_adapt_automaton_visitor);
	return my_adapt_automaton_visitor->get_bool_type();
}
reachability_scenario::coefficient_type reachability_scenario::get_number_type() const {
	assert(my_adapt_automaton_visitor);
	return my_adapt_automaton_visitor->get_number_type();
}

continuous_post_ptr reachability_scenario::get_continuous_post_operator() {
	return my_continuous_post_operator;
}

discrete_post_ptr reachability_scenario::get_discrete_post_operator() {
	return my_discrete_post_operator;
}

passed_and_waiting_list_ptr reachability_scenario::get_passed_and_waiting_list() {
	return my_pwl;
}

symbolic_state_collection_ptr reachability_scenario::get_symbolic_state_collection() {
	return my_symbolic_state_collection;
}

hybrid_automaton_ptr reachability_scenario::get_hybrid_automaton() {
	return my_hybrid_automaton;
}

hybrid_automaton_network_ptr reachability_scenario::get_hybrid_automaton_network() {
	return my_hybrid_automaton_network;
}

adapt_automaton_visitor_ptr reachability_scenario::get_adapt_automaton_visitor() {
	return my_adapt_automaton_visitor;
}

void reachability_scenario::set_continuous_post_operator(continuous_post_ptr p) {
	my_continuous_post_operator = p;
}

void reachability_scenario::set_discrete_post_operator(discrete_post_ptr p) {
	my_discrete_post_operator = p;
}

void reachability_scenario::set_passed_and_waiting_list(
		passed_and_waiting_list_ptr p) {
	my_pwl = p;
}

void reachability_scenario::set_symbolic_state_collection(
		symbolic_state_collection_ptr p) {
	my_symbolic_state_collection = p;
}

void reachability_scenario::set_hybrid_automaton(hybrid_automaton_ptr p) {
	my_hybrid_automaton = p;
}

void reachability_scenario::set_hybrid_automaton_network(
		hybrid_automaton_network_ptr p) {
	my_hybrid_automaton_network = p;
}

void reachability_scenario::set_adapt_automaton_visitor(
		adapt_automaton_visitor_ptr p) {
	my_adapt_automaton_visitor = p;
}

void reachability_scenario::apply_options(options::options_processor::variables_map& vmap) {
	if (my_option_handler) {
		my_option_handler(this,vmap);
	}
}

void reachability_scenario::set_option_handler(option_handler f) {
	my_option_handler=f;
}

void reachability_scenario::set_time_horizon(double th) {
	my_time_horizon = th;
}

double reachability_scenario::get_time_horizon() const {
	return my_time_horizon;
}

void reachability_scenario::set_sampling_time(double ts) {
	my_sampling_time = ts;
}

double reachability_scenario::get_sampling_time() const {
	return my_sampling_time;
}

void reachability_scenario::set_iter_max(int k) {
	my_iter_max = k;
}

int reachability_scenario::get_iter_max() const {
	return my_iter_max;
}

double reachability_scenario::get_intersection_error() const {
	return my_intersection_error;
}

void reachability_scenario::set_intersection_error(double e) {
	my_intersection_error = e;
}

std::string reachability_scenario::get_intersection_type() const {
	return my_intersection_type;
}
void reachability_scenario::set_intersection_type(std::string algo_type){
	my_intersection_type = algo_type;
}

std::string reachability_scenario::get_minbrak_type() const {
	return my_minbrak_type;
}
void reachability_scenario::set_minbrak_type(std::string minbrak_type){
	my_minbrak_type = minbrak_type;
}
size_t reachability_scenario::get_split_size() const {
	return my_split_size;
}
void reachability_scenario::set_split_size(size_t split_size) {
	assert(split_size >= 0);
	my_split_size = split_size;
}
std::string reachability_scenario::get_name() const {
	return my_name;
}

void reachability_scenario::set_name(const std::string& n) {
	my_name=n;
}

}
