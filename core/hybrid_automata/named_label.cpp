#include "core/hybrid_automata/named_label.h"

#include <stdexcept>
#include "utility/stl_helper_functions.h"

namespace hybrid_automata {

using namespace std;
std::map<string, label_id> named_label::string_to_label_id_chache =
		named_label::init_string_to_label_id_cache();
std::vector<string> named_label::label_id_to_string_cache =
		named_label::init_label_id_to_string_cache();

named_label::named_label() {
}

const label_id& named_label::silent_id() {
	static label_id sid = 1;
	return sid;
}

const std::string& named_label::silent_name() {
	static std::string sstr = "";
	return sstr;
}

label_id named_label::get_or_add_label_id(std::string new_name) {
	std::map<string, label_id>::const_iterator pos =
			string_to_label_id_chache.find(new_name);
	if (pos == string_to_label_id_chache.end()) {
		// create a new id
		label_id new_id = get_new_id();
		string_to_label_id_chache.insert(make_pair(new_name, new_id));
		label_id_to_string_cache.resize(get_label_count()+1);
		label_id_to_string_cache[new_id] = new_name;
		return new_id;
	} else // already in the cache
	{
		return pos->second;
	}
}

bool named_label::is_not_silent(label_id id)
{
	if (id > get_label_count())
		throw std::out_of_range("label id " + to_string(id) + " out of range");
	else if (id) {
		if(label_id_to_string_cache[id].compare("") == 0)
			return false;
		else
			return true;
	}
	else
		throw std::out_of_range("label id 0 has no name");
}

string named_label::get_name(label_id id) {
	// we use the fact that ids are given contiguously:
	// id is in the cache iff id<get_label_count()
	if (id > get_label_count())
		throw std::out_of_range("label id " + to_string(id) + " out of range");
	else if (id)
		return label_id_to_string_cache[id];
	else
		throw std::out_of_range("label id 0 has no name");
}

unsigned int named_label::get_label_count() {
	return label_id_to_string_cache.size();
}

void named_label::print_named_label_cache(std::ostream& o) {
	for (std::map<std::string, label_id>::const_iterator it =
			string_to_label_id_chache.begin(); it
			!= string_to_label_id_chache.end(); ++it) {
		o << it->first << " <-> " << it->second << endl;
	}
}

label_id named_label::get_new_id() {
	return label_id_to_string_cache.size();
}

void named_label::reset() {
	string_to_label_id_chache = init_string_to_label_id_cache();
	label_id_to_string_cache = init_label_id_to_string_cache();
}

std::map<std::string, label_id> named_label::init_string_to_label_id_cache() {
	// add the silent label
	std::map<std::string, label_id> m;
	m.insert(make_pair(silent_name(), silent_id()));
	return m;
}

std::vector<std::string> named_label::init_label_id_to_string_cache() {
	// add the silent label
	std::vector<std::string> v(silent_id() + 1);
	v[silent_id()] = silent_name();
	return v;
}

}
