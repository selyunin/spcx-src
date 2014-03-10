/*
 * positional_vdomain.h
 *
 *  Created on: Mar 25, 2010
 *      Author: frehse
 */

#ifndef positional_vdomain_H_
#define positional_vdomain_H_

#include "index_to_variable_id_map.h"

/** Forward declarations. */
class index_to_index_bimap;
namespace math {
template<typename T> class affine_map;
template<typename T> class vdom_vector;
}

typedef index_to_index_bimap position_map;

/** A named positional domain keeps track which variables are associated
 * to which space dimension.
 *
 * Variables can be added, removed, or remapped to other dimensions.
 */
class positional_vdomain {
public:
	typedef unsigned int size_type;

	/** Create a domain without any variables. */
	positional_vdomain();

	/** Create a domain with variables at the respective positions
	 * as given in the vector. */
	positional_vdomain(const std::vector<variable>& vars);

	/** Create a domain with variables at the respective positions. */
	positional_vdomain(const variable_id_list& vars);

	/** Create a domain with variables at arbitrary positions. */
	positional_vdomain(const variable_id_set& vars);

	/** Add a variable to the domain, assigning the next higher index in sequence. */
	void add_variable(const variable& x);

	/** Add variables to the domain. */
	void add_variables(const std::set<variable>& x);

	/** Remove a variable from the domain. */
	void remove_variable(const variable& x);

	/** Remove variables from the domain. */
	void remove_variables(const std::set<variable>& x);

	/** Get the variables in the domain. */
	std::set<variable> get_variables() const;

	/** Get a vector variables in the domain, ordered by their index. */
	std::vector<variable> get_variable_vector() const;

	/** Get the variable ids in the domain. */
	const variable_id_set& get_variable_ids() const;

	/** Returns the number of variables in the domain. */
	size_type size() const;

	/** Returns the position of variable x. */
	size_type pos(const variable& x) const;

	/** Returns the variable at position i. */
	variable get_variable(const size_type& i) const;

	/** Returns whether the variable x is in this domain. */
	bool in_domain(const variable& x) const;

	/** Returns true if the variable x is in this domain, and
	 * assigns its position to pos. */
	bool in_domain(const variable& x, size_type& pos) const;

	/** Check if the variables in dom are contained in *this.
	 *
	 * The check is independent of the position of each variable. */
	bool contains_variables(const positional_vdomain& dom) const;

	/** Swap */
	void swap(positional_vdomain& D);

	friend void swap(positional_vdomain& D1, positional_vdomain& D2);

	/** Domain comparison */
	bool operator==(const positional_vdomain& D) const;

	/** Domain comparison */
	bool operator!=(const positional_vdomain& D) const;

	/** Output to stream */
	void print(std::ostream& os) const;

	/** Construct using index_to_variable_id_map. */
	positional_vdomain(index_to_variable_id_map_ptr iimap);

	/** Obtain index_to_variable_id_map. */
	const index_to_variable_id_map_ptr& get_index_to_variable_id_map() const;

	/** Returns the number representing an invalid position.
	 *
	 * Used for example in reorderng.
	 */
	static size_type invalid_pos();
private:


	friend class index_to_variable_id_map_provider;
	template<typename T> friend class math::affine_map;
	template<typename T> friend class math::vdom_vector;
	friend positional_vdomain compose(const positional_vdomain& D1,
			const positional_vdomain& D2, position_map& f1,
			position_map& f2);

	/** The index_to_variable_id_map is how the domain internally keeps track
	 * of things.
	 */
	index_to_variable_id_map_ptr my_iimap;
};

void swap(positional_vdomain& D1, positional_vdomain& D2);

/** Yields the position_map that reorders the elements of D1
 * to be elements of D2.
 *
 * If use_invalid_pos=true then invalid_pos is inserted if an element of D1
 * is not in D2. Otherwise throws.
 */
position_map compute_reordering(const positional_vdomain& D1,
		const positional_vdomain& D2, bool use_invalid_pos=false);

/** Join two domains D1 and D2 into a single one D.
 *
 * Returns D. Assigns to f1 the map that reorders elements of D1 into D,
 * and to f2 the reordering of D2 into D.
 */
positional_vdomain compose(const positional_vdomain& D1,
		const positional_vdomain& D2, position_map& f1, position_map& f2);

/** Join two domains D1 and D2 into a single one D.
 *
 * Returns D.
 */
positional_vdomain compose(const positional_vdomain& D1,
		const positional_vdomain& D2);

/**
 * Output as a stream of characters. Calls the print method.
 */
std::ostream& operator<<(std::ostream& os, const positional_vdomain& D);

/** Utility function to convert variable_id_sets to sets of variables. */
std::set<variable> create_variable_set(const variable_id_set& vis);

/** Utility function to convert sets of variables to variable_id_sets. */
variable_id_set create_variable_id_set(const std::set<variable>& vars);

#include "positional_vdomain.hpp"

#endif /* positional_vdomain_H_ */
