#include "math/vdom/index_to_variable_id_map.h"

#include <set>
#include <stdexcept>
#include "utility/stl_helper_functions.h"
//#include "shared_ptr_user.h"
//#include "../set_implementations/ppl_polyhedron/general.h"
#include "math/vdom/variable.h"
#include "math/vdom/index_to_index_bimap.h"

using namespace std;

// Defining static members
std::map<std::vector<variable_id>, index_to_variable_id_map_ptr>
		index_to_variable_id_map::index_to_variable_id_map_cache;

//std::map<std::vector<variable_id>,index_to_variable_id_map_ptr> index_to_variable_id_map_cache=std::map<std::vector<variable_id>,index_to_variable_id_map_ptr>();



index_to_variable_id_map_ptr index_to_variable_id_map::get_map_with_ids(
		const variable_id_set& new_ids) {
	// construct the new vector
	vector<variable_id> new_v(new_ids.size());
	vector<variable_id>::iterator it = new_v.begin();
	copy(new_ids.begin(), new_ids.end(), it);

	return get_index_to_variable_id_map_ptr(new_v);
}

index_to_variable_id_map_ptr index_to_variable_id_map::get_map_with_id_added(
		const variable_id& id) const {
	// don't add an id that already exists
	map<variable_id, index_type>::const_iterator pos = my_id_to_index_map.find(
			id);
	if (pos == my_id_to_index_map.end()) {
		vector<variable_id> new_v(dimensions() + 1);
		vector<variable_id>::iterator i = new_v.begin();
		copy(my_id_vector.begin(), my_id_vector.end(), i);
		(*new_v.rbegin()) = id; // assign id as the last element
		return get_index_to_variable_id_map_ptr(new_v);
	} else {
		return shared_from_this();
	}
}

index_to_variable_id_map_ptr index_to_variable_id_map::get_map_with_ids_added(
		const variable_id_set& id_set) const {
	// don't add ids that already exist
	variable_id_set new_ids;
	variable_id_set::const_iterator it = id_set.begin();
	while (it != id_set.end()) {
		map<variable_id, index_type>::const_iterator pos =
				my_id_to_index_map.find(*it);
		if (pos == my_id_to_index_map.end()) {
			new_ids.insert(*it);
		}
		++it;
	}

	if (new_ids.size() > 0) {
		//cout << "new ids:" << new_ids;
		// construct the new vector
		vector<variable_id> new_v(dimensions() + new_ids.size());
		vector<variable_id>::iterator i = new_v.begin();
		copy(my_id_vector.begin(), my_id_vector.end(), i);
		//cout << "before new:" << new_v;
		copy(new_ids.begin(), new_ids.end(), i + dimensions());
		//cout << "after new:" << new_v;

		return get_index_to_variable_id_map_ptr(new_v);
	} else {
		// there's nothing to do because there are no new id's

		// return a pointer to *this
		// see: http://www.boost.org/doc/libs/1_35_0/libs/smart_ptr/sp_techniques.html#from_this
		return shared_from_this();
	}
}

index_to_variable_id_map_ptr index_to_variable_id_map::get_map_with_ids_removed(
		const variable_id_set& id_set) const {
	// check if id_set is not empty
	if (id_set.empty()) {
		return shared_from_this();
	} else {
		// count how many ids will be removed
		unsigned int n_removed_ids = 0;
		bool found;
		for (variable_id_set::const_iterator it = id_set.begin(); it
				!= id_set.end(); ++it) {
			check_for_index(*it, found);
			if (found) {
				++n_removed_ids;
			}
		}

		// create the new vector
		std::vector<variable_id> new_v(dimensions() - n_removed_ids);
		index_type j = 0;
		for (index_type i = 0; i < my_id_vector.size() && j < new_v.size(); ++i) {
			if (id_set.find(get_id(i)) == id_set.end()) {
				new_v[j] = get_id(i);
				++j;
			}
		}
		return get_index_to_variable_id_map_ptr(new_v);
	}
}

index_to_variable_id_map_ptr index_to_variable_id_map::get_map_with_primedness_reassigned(
		unsigned int d, unsigned int p) const {
	// create the new vector
	std::vector<variable_id> new_v(my_id_vector.size());
	for (index_type i = 0; i < my_id_vector.size(); ++i) {
		if (variable::get_prime_count(my_id_vector[i]) == d)
			new_v[i] = variable::get_primed_id(my_id_vector[i], p);
	}
	return get_index_to_variable_id_map_ptr(new_v);
}

index_to_variable_id_map_ptr index_to_variable_id_map::get_map_with_primedness_decreased(
		unsigned int d) const {
	// create the new vector
	std::vector<variable_id> new_v(my_id_vector.size());
	for (index_type i = 0; i < my_id_vector.size(); ++i) {
		if (d == 0 || variable::get_prime_count(my_id_vector[i]) == d)
			new_v[i] = variable::get_id_primedness_decreased(
					my_id_vector[i]);
	}
	return get_index_to_variable_id_map_ptr(new_v);
}

index_to_variable_id_map_ptr index_to_variable_id_map::get_map_with_primedness_increased(
		unsigned int d) const {
	// create the new vector
	std::vector<variable_id> new_v(my_id_vector.size());
	for (index_type i = 0; i < my_id_vector.size(); ++i) {
		if (d == 0 || variable::get_prime_count(my_id_vector[i]) == d)
			new_v[i] = variable::get_id_primedness_increased(
					my_id_vector[i]);
	}
	return get_index_to_variable_id_map_ptr(new_v);
}

index_to_variable_id_map_ptr index_to_variable_id_map::get_index_to_variable_id_map_ptr(
		const std::vector<variable_id>& new_v) {

	// look in the cache if it's already there
	std::map<std::vector<variable_id>, index_to_variable_id_map_ptr>::const_iterator
			p = index_to_variable_id_map_cache.find(new_v);
	if (p == index_to_variable_id_map_cache.end()) {
		// generate new id_to_index_map
		std::map<variable_id, index_type> new_id_to_index_map;

		for (index_type i = 0; i < new_v.size(); ++i) {
			new_id_to_index_map.insert(make_pair(new_v[i], i));
		}

		//cout << "generated id_to_index_map:" << new_id_to_index_map;
		index_to_variable_id_map_ptr new_ptr = index_to_variable_id_map_ptr(
				new index_to_variable_id_map(new_v, new_id_to_index_map, compute_id_set(new_v)));
		//cout << "new iim in cache:" << new_ptr << "!" << endl;
		index_to_variable_id_map_cache.insert(make_pair(new_v, new_ptr));
		return new_ptr;
		// The following returns the exact same pointer, but this is a little paranoid. Equality of shared_ptr is decided based on the pointer contained, so it doesn't matter.
		//pair<std::map<std::vector<variable_id>,index_to_variable_id_map_ptr>::iterator,bool> ret;
		//ret=index_to_variable_id_map_cache.insert(make_pair(new_v,new_ptr));
		//return ret.first->second; // return exactly the pointer, so map equality can be decided based on pointer equality
	} else {
		return p->second;
	}
}

variable_id_set index_to_variable_id_map::compute_id_set(const std::vector<variable_id>& id_vec) {
	variable_id_set vis;
	for (std::vector<variable_id>::const_iterator it = id_vec.begin(); it
			!= id_vec.end(); ++it) {
		vis.insert(*it);
	}
	return vis;
}

bool index_to_variable_id_map::contains_variables(
		const index_to_variable_id_map_ptr& iimap) const {
	for (std::vector<variable_id>::const_iterator it =
			iimap->my_id_vector.begin(); it != iimap->my_id_vector.end(); ++it) {
		if (!has_id(*it)) {
			return false;
		}
	}
	return true;
}

void index_to_variable_id_map::print(std::ostream& os) const {
	for (std::map<variable_id, index_type>::const_iterator it =
			my_id_to_index_map.begin(); it != my_id_to_index_map.end(); ++it) {
		os << it->second << "(index)" << " <-> " << it->first << endl;
	}
	//os << my_id_vector;
	/* For printing the cache. (doesn't work yet)
	 for (std::map<std::vector<variable_id>,index_to_variable_id_map_ptr>::const_iterator
	 p = index_to_variable_id_map_cache.begin(); p!= index_to_variable_id_map_cache.end(); ++p)
	 {
	 os << p->first << " <-> " << p->second << "(id)" << endl;
	 }
	 */
}

void get_common_map(const index_to_variable_id_map_ptr& M1,
		const index_to_variable_id_map_ptr& M2,
		index_to_variable_id_map_ptr& M, index_type& newdim,
		index_to_index_bimap& iimap2) {

	// We construct the new id_vector, and instantiate M from that.
	// First, we build the new iimap2 to get newdim

	// Construct iimap2
	iimap2 = index_to_index_bimap();

	newdim = M1->dimensions();
	variable_id new_id;
	index_type i1, new_index;
	bool found;
	set<index_type> new_from_M2;
	// Make a first pass to get the new dimension
	// Count the ids that aren't already there
	new_index = newdim; // the index at which new variables will be inserted
	for (dimension_t i = 0; i < M2->dimensions(); ++i) {
		new_id = M2->get_id(i);
		i1 = M1->check_for_index(new_id, found);
		if (found) {
			iimap2.insert(i, i1);
		} else {
			new_from_M2.insert(i);
			iimap2.insert(i, new_index);
			++newdim;
			new_index = newdim;
		}
	}
	iimap2.fill_up_to(newdim); // add the remaining indices up to newdim-1

	// Construct M
	vector<variable_id> new_v(newdim);
	copy(M1->my_id_vector.begin(), M1->my_id_vector.end(), new_v.begin()); // fill with M1
	new_index = M1->dimensions(); // the index at which new variables will be inserted
	for (set<index_type>::const_iterator it = new_from_M2.begin(); it
			!= new_from_M2.end(); ++it) {
		new_v[iimap2.get_map(*it)] = M2->get_id(*it);
	}

	M = index_to_variable_id_map::get_index_to_variable_id_map_ptr(new_v);
}

index_to_index_bimap get_index_to_index_mapping(
		const index_to_variable_id_map_ptr& to_M,
		const index_to_variable_id_map_ptr& from_M) {

	// Construct iimap
	index_to_index_bimap iimap = index_to_index_bimap();
	variable_id vid;
	index_type to_index;
	bool found;
	for (dimension_t i = 0; i < from_M->dimensions(); ++i) {
		vid = from_M->get_id(i);
		to_index = to_M->check_for_index(vid, found);
		if (found)
			iimap.insert(i, to_index);
	}
	return iimap;
}

std::ostream& operator<<(std::ostream& os,
		const index_to_variable_id_map_ptr dsp) {
	dsp->print(os);
	return os;
}
