/*
 * automaton_cache.h
 *
 *  Created on: Aug 25, 2009
 *      Author: frehse
 */

#ifndef AUTOMATON_CACHE_H_
#define AUTOMATON_CACHE_H_

#include <string>
#include <map>
#include <set>
#include <stdexcept>
#include "boost/shared_ptr.hpp"
#include "core/hybrid_automata/automaton_id.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
//typedef boost::shared_ptr<const discrete_set> const_ptr;
}

namespace hybrid_automata {

/** Stores automata globally.
 */

class hybrid_automaton_cache {

public:
	/** Add an automaton to the cache. */
	static void add_automaton(hybrid_automaton_ptr p);

	/** Remove an automaton from the cache. */
	static void remove_automaton(automaton_id id);

	/** Remove all automata from the cache. */
	static void clear();

	/** Get a pointer to the automaton corresponding to id. Returns automaton::ptr() if the id doesn't exist. */
	static hybrid_automaton_ptr get_automaton(automaton_id id);

	/** Get a pointer to the automaton corresponding to name. Returns automaton::ptr() if the id doesn't exist.
	 *\todo check name is single
	 * */
	static hybrid_automaton_ptr get_automaton(const std::string& name);

	/** Get the id of the automaton corresponding to name. Throws if the name doesn't exist.
	 * */
	static automaton_id get_automaton_id(const std::string& name);

	/** Returns true if an automaton with id \p id exists in the cache.
	 * */
	static bool has_automaton(automaton_id id);

	/** Returns true if the automaton *p exists in the cache.
	 * */
	static bool has_automaton(const hybrid_automaton_ptr& p);

	/** Returns true if an automaton corresponding to name exists in the cache.
	 * */
	static bool has_automaton(const std::string& name);

	/** Swap the identities of automata p1 and p2, i.e., exchange their
	 * names and ids.
	 * All future calls to p1 that pass via ids will so be redirected to p2 and
	 * vice versa.
	 *
	 * @note The swap function takes pointers as arguments because it needs
	 * to work even if the automata are not in the cache.
	 * (This needn't be a member function, actually.)
	 */
	static void swap_identity(hybrid_automaton_ptr p1, hybrid_automaton_ptr p2);

	/** Output as a stream of characters.
	 */
	static void print(std::ostream& os);

	/** Output as a stream of characters.
	 */
	static void print_all_automata(std::ostream& os);

	/** Retrieve the set of all automaton_id in the cache. */
	static std::set<automaton_id> get_automata();

	/** Retrieve the set of all names in the cache. */
	static std::set<std::string> get_automaton_names();

private:
	typedef std::map<automaton_id, hybrid_automaton_ptr> container_type;
	static container_type automaton_id_to_ptr_map;
};

}

#endif /* AUTOMATON_CACHE_H_ */
