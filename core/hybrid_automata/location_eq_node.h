/*
 * location_eq_node.h
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#ifndef LOCATION_EQ_NODE_H_
#define LOCATION_EQ_NODE_H_

#include "core/hybrid_automata/location_id.h"
#include "core/hybrid_automata/automaton_id.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/predicates/dot_context_lookup.h"
#include "utility/tree_node.h"

//#include "../utility/tree_node_visitor.h"
//#include "../valuation_function/valuation_function_tree_nodes.h"

namespace hybrid_automata {

/**
 * \todo create base class for loc_node
 */
class location_eq_node: public tree::node {
public:
	location_eq_node(const std::string& aut, const std::string& loc, bool equal);
	location_eq_node(automaton_id aut, location_id loc, bool equal);
	virtual ~location_eq_node();

	const automaton_id& get_automaton_id() const;

	const location_id& get_location_id() const;
	const bool& get_equal() const;
	virtual void accept(tree::node_visitor& v) const;
private:
	automaton_id my_automaton;
	location_id my_location;
	bool my_equal;
};

/** A class that controls the creation/lookup of automaton names.
 *
 * Having it as a class gives us the ability to add
 * static memory, e.g., for locking the creation of variables.
 *
 * When parsing a model file, variables are already created
 * when the symbol table is created, so we need lookup only.
 * A variable that is not found should generate an exception.
 *
 * However, when we want to parse a simple string inside
 * a tester or elsewhere, we might not want to bother
 * with a symbol table etc. For this the creation
 * should be unlocked.
 * */
class location_node_creator {
public:
	static location_eq_node* create(const std::string& aut_name,
			const std::string& context, const std::string& loc_name, bool equal) {
		std::string namec = name_in_context(aut_name, context);
		if (locked) {
			namec = context_lookup(aut_name, context);
		}
		return new location_eq_node(namec, loc_name, equal);
	}
	;
	static bool is_locked() {
		return locked;
	}
	static void set_locked(bool l) {
		locked = l;
	}
	/** Return the name within the context */
	static std::string name_in_context(const std::string& name,
			const std::string& context) {
		if (context.empty()) {
			return name;
		} else {
			return context + "." + name;
		}
	}
private:
	static std::string context_lookup(const std::string& name,
			const std::string& context) {
		std::string namec = name_in_context(name, context);
		// try lookup inside context
		if (!name.empty() && !hybrid_automaton_cache::has_automaton(namec)) {
			std::set<std::string> known =
					hybrid_automaton_cache::get_automaton_names();
			std::set<std::string> candidates = dot_context::lookup(name,
					context, known.begin(), known.end());
			if (candidates.size() == 1) {
				namec = name_in_context(*candidates.begin(), context);
			} else if (candidates.size() == 0) {
				throw basic_exception("Could not find component " + name
						+ " in component " + context + ".");
			} else {
				std::string s;
				for (std::set<std::string>::const_iterator it =
						candidates.begin(); it != candidates.end(); ++it) {
					if (it != candidates.begin()) {
						s += ", ";
					}
					s += *it;
				}
				throw basic_exception("Could not find component " + name
						+ " in component " + context + ". Found instead " + s
						+ ". Please use a more explicit identifier.");
			}
		}
		return namec;
	}

	static bool locked;
};

}

#endif /* LOCATION_EQ_NODE_H_ */
