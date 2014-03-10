/*
 * index_to_variable_id_map_provider.hpp
 *
 *  Created on: May 12, 2010
 *      Author: frehse
 */

#ifndef INDEX_TO_VARIABLE_ID_MAP_PROVIDER_HPP_
#define INDEX_TO_VARIABLE_ID_MAP_PROVIDER_HPP_

inline index_to_variable_id_map_provider::index_to_variable_id_map_provider() {
}

/** Constructor with given iimap. */
inline index_to_variable_id_map_provider::index_to_variable_id_map_provider(index_to_variable_id_map_ptr iimap) :
	my_domain(iimap) {
}

inline index_to_variable_id_map_provider::index_to_variable_id_map_provider(
		positional_vdomain dom) :
		my_domain(dom) {
}


inline index_to_variable_id_map_provider::~index_to_variable_id_map_provider() {
}

/** Returns the \p variable_id at index i.
 */
inline variable_id index_to_variable_id_map_provider::index_to_variable_id_map_provider::get_id(index_type i) const {
	return get_index_to_variable_id_map()->get_id(i);
}

/** Returns the index of \p variable_id id.
 */
inline bool index_to_variable_id_map_provider::index_to_variable_id_map_provider::has_id(variable_id id) const {
	bool has=false;
	get_index_to_variable_id_map()->check_for_index(id, has);
	return has;
}

/** Returns the index of \p variable_id id.
 */
inline variable_id index_to_variable_id_map_provider::index_to_variable_id_map_provider::get_index(variable_id id) const {
	return get_index_to_variable_id_map()->get_index(id);
}

inline index_type index_to_variable_id_map_provider::index_to_variable_id_map_provider::check_for_index(const variable_id& id, bool& has_id) const {
	return get_index_to_variable_id_map()->check_for_index(id,has_id);
}

/** Returns the ids of all variables in the index_to_variable_id_map. */
inline const variable_id_set& index_to_variable_id_map_provider::index_to_variable_id_map_provider::get_variable_ids() const {
	return my_domain.get_variable_ids();
}

/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
 */
inline void index_to_variable_id_map_provider::reassign_primedness(unsigned int d, unsigned int p) {
	my_domain
			= positional_vdomain(get_index_to_variable_id_map()->get_map_with_primedness_reassigned(d, p));
}

/** Increase the primedness of the variables with primedness of degree \p d by 1.
 * If d is 0, increase all. */
inline void index_to_variable_id_map_provider::increase_primedness(unsigned int d) {
	my_domain
				= positional_vdomain(get_index_to_variable_id_map()->get_map_with_primedness_increased(d));
}

/** Decrease the primedness of the variables with primedness of degree \p d by 1.
 * If d is 0, decrease all. */
inline void index_to_variable_id_map_provider::decrease_primedness(unsigned int d) {
	my_domain
				= positional_vdomain(get_index_to_variable_id_map()->get_map_with_primedness_decreased(d));
}

inline variable_id_set index_to_variable_id_map_provider::get_primed_variables(unsigned int prime_count) const {
	variable_id_set vis;
	variable_id id;
	// Check for all indices whether their corresponding id is of prime degree prime_count
	for (index_type i=0; i< my_domain.size(); ++i) {
		id = get_id(i);
		if (variable::get_prime_count(id)==prime_count) {
			vis.insert(id);
		}
	}
	return vis;
}

/*! Returns a pointer to the index_to_variable_id_map of \p *this. */
inline const index_to_variable_id_map_ptr& index_to_variable_id_map_provider::get_index_to_variable_id_map() const {
	return my_domain.get_index_to_variable_id_map();
}

/** Assigns pnew_map to the index_to_variable_id_map of \p *this.
 @attention No checks for consistency are applied.
 */
inline void index_to_variable_id_map_provider::set_index_to_variable_id_map(
		const index_to_variable_id_map_ptr& pnew_map) {
	my_domain = positional_vdomain(pnew_map);
	//std::cout << _index_to_variable_id_map_ptr;
}

inline
void index_to_variable_id_map_provider::set_domain(const positional_vdomain& dom) {
	my_domain=dom;
}

inline
const positional_vdomain& index_to_variable_id_map_provider::domain() const {
	return my_domain;
}

inline
std::set<variable> index_to_variable_id_map_provider::get_variables() const {
	return domain().get_variables();
}

inline
index_type index_to_variable_id_map_provider::pos(const variable& x) const {
	return get_index(x.get_id());
}

inline
variable index_to_variable_id_map_provider::get_variable(const index_type& i) const {
	return variable(get_id(i));
}

#endif /* INDEX_TO_VARIABLE_ID_MAP_PROVIDER_HPP_ */
