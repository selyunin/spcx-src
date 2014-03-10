/*
 * dot_context_lookup.h
 *
 *  Created on: Oct 19, 2010
 *      Author: frehse
 */

#ifndef DOT_CONTEXT_LOOKUP_H_
#define DOT_CONTEXT_LOOKUP_H_

#include <iostream>
#include <string>
#include <set>
#include "utility/stl_helper_functions.h"
#include "math/vdom/variable.h"

namespace dot_context {

/** Returns true if name begins with the context. */
bool starts_with_context(const std::string& name, const std::string context);

/** Returns the name within the the context.
 *
 * If the context was not found, returns the unmodified string. */
std::string in_context_name(const std::string& name, const std::string context);

/** Returns the name stripped of any context. */
std::string context_free_name(const std::string& name);

/** Returns the names in the container that match the context.
 *
 * @attention This function assumes that context.name is not in the store.
 * The caller should verify that this is the case, and use
 * in_context_name otherwise.
 *
 * A string matches name and context if it begins with the context
 * and ends with "." followed by the name.
 */
template<class const_iterator> std::set<std::string> lookup(
		const std::string& name, const std::string context,
		const_iterator ibeg, const const_iterator& iend) {
	std::set<std::string> res;
	std::string dotname = "." + name;
	size_t dns = dotname.size();
	size_t cs = context.size();
	for (; ibeg != iend; ++ibeg) {
		if (starts_with_context(*ibeg, context) && ibeg->size()>=dns && ibeg->substr(ibeg->size()
				- dns, dns) == dotname) {
			res.insert(ibeg->substr(cs + 1)); // return the name without the context
		}
	}
	return res;
}
;

/** Strip the name of any context.
 *
 * @attention This function assumes that name.context is not in the store.
 * The caller should verify that this is the case, and use
 * in_context_name otherwise.
 *
 */
template<class const_iterator>
std::string smallest_unique_name(const std::string& name,
		const std::string& context, const_iterator ibeg,
		const const_iterator& iend) {
	std::string good_name = name;
	std::string candidate_name = good_name;
	std::set<std::string> s;
	size_t cs;
	bool found_better=false;
	do {
		found_better=false;
		cs = candidate_name.find(".");
		if (cs != std::string::npos && cs+1<candidate_name.size()) {
			candidate_name = candidate_name.substr(cs + 1);
			s = lookup(candidate_name, context, ibeg, iend);
			if (s.size() == 1) {
				good_name = candidate_name;
				found_better=true;
			}
			cs = candidate_name.find(".");
		}
	} while (found_better);
	return good_name;
}

/** Check if the set of variables v1 contains those in v2.
 *
 * If not, throw a basic_exception with the message prologue+variables+epilogue.
 * The variables are reported context_free.
 */
void throw_if_not_contains_context_free(const variable_id_set& v1,
		const variable_id_set& v2, std::string prologue, std::string epilogue);

}

#endif /* DOT_CONTEXT_LOOKUP_H_ */
