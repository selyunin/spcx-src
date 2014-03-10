#ifndef named_label_H_
#define named_label_H_

#include <set>
#include <vector>
#include <map>
#include <string>
#include <iostream>
//#include <iostream>

namespace hybrid_automata {

typedef unsigned int label_id;
typedef std::set<label_id> label_id_set;

/** A class for attributing label names (string) to identifiers of type \p label_id.
 * The association is cached in a bidirectional map.
 *
 * There is a special label called the silent label, with a fixed id different from 0.
 * Its string representation is the empty string.
 *
 * The id 0 can be used similar to a null pointer.
 */
class named_label {

public:
	/** Obtain a unique identifier for the label name \p new_name.
	 * If the identifier is not found in the cache, then a new one is created.
	 */
	static label_id get_or_add_label_id(std::string new_name);

	/** Returns the name associated with the identifier \p id.
	 */
	static std::string get_name(label_id id);

	/** Returns the number of labels in the cache.
	 */
	static unsigned int get_label_count();

	/** Output the list of names and associated ids. */
	static void print_named_label_cache(std::ostream& o);

	/** Clear the cache of all labels. */
	static void reset();

	/** Get the id of the silent label. */
	static const label_id& silent_id();

	/** Get the name of the silent label. */
	static const std::string& silent_name();

	/** If silent label, return false */
	static bool is_not_silent(label_id id);

private:
	named_label();
	/** Helper function that returns an id that is guaranteed to be new. */
	static label_id get_new_id();

	/** Init string_to_label_id_chache */
	static std::map<std::string,label_id> init_string_to_label_id_cache();

	/** Init label_id_to_string_cache */
	static std::vector<std::string> init_label_id_to_string_cache();

	/** Init label_count */
	static unsigned int init_label_count();

	static std::map<std::string,label_id> string_to_label_id_chache;
	static std::vector<std::string> label_id_to_string_cache;
};

}

#endif /*named_label_H_*/
