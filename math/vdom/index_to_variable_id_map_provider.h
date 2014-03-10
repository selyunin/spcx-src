#ifndef INDEX_TO_VARIABLE_ID_MAP_PROVIDER_H_
#define INDEX_TO_VARIABLE_ID_MAP_PROVIDER_H_

#include "math/vdom/index_to_variable_id_map.h"
#include "math/vdom/primed_variable_provider.h"
#include "positional_vdomain.h"

class index_to_variable_id_map_provider : public virtual primed_variable_provider {
public:
	/** Constructor with empty iimap. */
	index_to_variable_id_map_provider();

	/** Constructor with given iimap. */
	index_to_variable_id_map_provider(index_to_variable_id_map_ptr iimap);

	/** Constructor with given domain. */
	index_to_variable_id_map_provider(positional_vdomain dom);

	virtual ~index_to_variable_id_map_provider();

	// --------------------------------------------
	/** \name Non-modifying methods
	 *  \{ */
	// --------------------------------------------

	/** Returns the \p variable_id at index i. */
	variable_id get_id(index_type i) const;

	/** Returns the index of \p variable_id id. */
	bool has_id(variable_id id) const;

	/** Returns the index of \p variable_id id. */
	index_type get_index(variable_id id) const;

	/** Returns the index belonging to the id \p id if there is one.
	 * has_id is true if the id was found, and false otherwise.
	 * If has_id is false, then the returned index_type defaults to 0.
	 */
	index_type check_for_index(const variable_id& id, bool& has_id) const;

	/** Returns the ids of all variables in the index_to_variable_id_map. */
	const variable_id_set& get_variable_ids() const;

	/** Returns the ids of the variables that are primed to degree \p
	 * prime_count. */
	variable_id_set get_primed_variables(unsigned int prime_count) const;

	/* \} */
	// --------------------------------------------
	/** \name Methods using the domain interface
	 *  \{ */
	// --------------------------------------------

	/** Get the domain. */
	const positional_vdomain& domain() const;

	/** Get the variables in the domain. */
	std::set<variable> get_variables() const;

	/** Returns the position of variable x. */
	index_type pos(const variable& x) const;

	/** Returns the variable at position i. */
	variable get_variable(const index_type& i) const;

	/* \} */
	// --------------------------------------------
	/** \name Methods changing the primedness of variables
	 *  \{ */
	// --------------------------------------------

	/** Set the primedness of the variables with primedness of degree \p d to
	 * degree \p p. */
	virtual void reassign_primedness(unsigned int d, unsigned int p = 0);

	/** Increase the primedness of the variables with primedness of
	 * degree \p d by 1.
	 * If d is 0, increase all. */
	virtual void increase_primedness(unsigned int d = 0);

	/** Decrease the primedness of the variables with primedness of
	 * degree \p d by 1.
	 * If d is 0, decrease all. */
	virtual void decrease_primedness(unsigned int d = 0);

	/* \} */
	// --------------------------------------------
	/** \name Low level functions
	 *  \{ */
	// --------------------------------------------

	/** Returns a pointer to the index_to_variable_id_map of \p *this. */
	const index_to_variable_id_map_ptr& get_index_to_variable_id_map() const;

	/** Swap with another index_to_variable_id_map_provider. */
	void swap(index_to_variable_id_map_provider& v) {
		if (this != &v) {
			//swap(_index_to_variable_id_map_ptr,v._index_to_variable_id_map_ptr);
			my_domain.swap(v.my_domain);
		}
	}
	;

	/** Swap two index_to_variable_id_map_providers. */
	friend void swap(index_to_variable_id_map_provider &v1, index_to_variable_id_map_provider &v2) {
		v1.swap(v2);
	}
	;

	/* \} */
	// --------------------------------------------

protected:
	/** Assigns pnew_map to the index_to_variable_id_map of \p *this.
	 @attention No checks for consistency are applied. */
	void set_index_to_variable_id_map(const index_to_variable_id_map_ptr& pnew_map);

	/** Assigns dom to *this.
	 @attention No checks for consistency are applied. */
	void set_domain(const positional_vdomain& dom);

private:
	//index_to_variable_id_map_ptr _index_to_variable_id_map_ptr;
	positional_vdomain my_domain;
};

#include "index_to_variable_id_map_provider.hpp"

#endif /*INDEX_TO_VARIABLE_ID_MAP_PROVIDER_H_*/
