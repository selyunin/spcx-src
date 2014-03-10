/*
 * automaton_cache.cpp
 *
 *  Created on: Aug 25, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/automaton_cache.h"

#include <stdexcept>
#include "utility/stl_helper_functions.h"
//#include"../../utility/shared_ptr_output.h"
#include "core/hybrid_automata/hybrid_automaton.h"

namespace hybrid_automata {

void hybrid_automaton_cache::add_automaton(hybrid_automaton_ptr p) {
	automaton_id_to_ptr_map[p->get_id()] = p;
}

void hybrid_automaton_cache::remove_automaton(automaton_id id) {
	automaton_id_to_ptr_map.erase(id);
}

void hybrid_automaton_cache::clear() {
	automaton_id_to_ptr_map.clear();
}

hybrid_automaton_ptr hybrid_automaton_cache::get_automaton(automaton_id id) {
	container_type::const_iterator it = automaton_id_to_ptr_map.find(id);
	if (it != automaton_id_to_ptr_map.end())
		return it->second;
	else
		return hybrid_automaton_ptr();
}

hybrid_automaton_ptr hybrid_automaton_cache::get_automaton(const std::string& name) {
	for (container_type::const_iterator it = automaton_id_to_ptr_map.begin(); it
			!= automaton_id_to_ptr_map.end(); ++it) {
		if (it->second->get_name() == name)
			return it->second;
	}
	return hybrid_automaton_ptr();
}

automaton_id hybrid_automaton_cache::get_automaton_id(const std::string& name) {
	hybrid_automaton::ptr p=get_automaton(name);
	if (!p)
		throw std::runtime_error("Could not find automaton with name '" + name + "'.");
	return p->get_id();
}

bool hybrid_automaton_cache::has_automaton(automaton_id id) {
	return automaton_id_to_ptr_map.find(id)!=automaton_id_to_ptr_map.end();
}

bool hybrid_automaton_cache::has_automaton(const std::string& name) {
	for (container_type::const_iterator it = automaton_id_to_ptr_map.begin(); it
			!= automaton_id_to_ptr_map.end(); ++it) {
		if (it->second->get_name() == name)
			return true;
	}
	return false;
}

bool hybrid_automaton_cache::has_automaton(const hybrid_automaton_ptr& p) {
	for (container_type::const_iterator it = automaton_id_to_ptr_map.begin(); it
			!= automaton_id_to_ptr_map.end(); ++it) {
		if (it->second == p)
			return true;
	}
	return false;
}

void hybrid_automaton_cache::swap_identity(hybrid_automaton_ptr p1, hybrid_automaton_ptr p2) {
	assert(p1);
	assert(p2);

//	print(std::cerr);

	// Note: don't put a new automaton in the cache.
	automaton_id id1=p1->my_id;
	automaton_id id2=p2->my_id;

	bool has_p1=has_automaton(id1);
	bool has_p2=has_automaton(id2);

	std::swap(p1->my_id, p2->my_id);
	std::swap(p1->my_name, p2->my_name);
	// Fix the cache
	// 1) remove the old associations
	remove_automaton(id1);
	remove_automaton(id2);
	// 2) add the new ones
	if (has_p2)
	add_automaton(p1);
	if (has_p1)
	add_automaton(p2);

//	std::cerr << "after" << std::endl;
//	print(std::cerr);
}

void hybrid_automaton_cache::print(std::ostream& os) {
	//os << automaton_id_to_ptr_map;
	for (container_type::const_iterator it = automaton_id_to_ptr_map.begin(); it
			!= automaton_id_to_ptr_map.end(); ++it) {
		if (it != automaton_id_to_ptr_map.begin())
			os << std::endl;
		os << it->first << " --> " << it->second->get_name();
	}
}

void hybrid_automaton_cache::print_all_automata(std::ostream& os) {
	//os << automaton_id_to_ptr_map;
	for (container_type::const_iterator it = automaton_id_to_ptr_map.begin(); it
			!= automaton_id_to_ptr_map.end(); ++it) {
		if (it != automaton_id_to_ptr_map.begin()) {
			os << std::endl;
			os << std::endl;
			os << std::endl;
		}
		os << it->second;
	}
}

std::set<automaton_id> hybrid_automaton_cache::get_automata() {
	//os << automaton_id_to_ptr_map;
	std::set<automaton_id> aut_set;
	for (container_type::const_iterator it = automaton_id_to_ptr_map.begin(); it
			!= automaton_id_to_ptr_map.end(); ++it) {
		aut_set.insert(it->first);
	}
	return aut_set;
}

std::set<std::string> hybrid_automaton_cache::get_automaton_names() {
	std::set<std::string> aut_set;
	for (container_type::const_iterator it = automaton_id_to_ptr_map.begin(); it
			!= automaton_id_to_ptr_map.end(); ++it) {
		aut_set.insert(it->second->get_name());
	}
	return aut_set;
}

std::map<automaton_id, hybrid_automaton_ptr> hybrid_automaton_cache::automaton_id_to_ptr_map;

}
