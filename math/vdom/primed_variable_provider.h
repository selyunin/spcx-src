#ifndef PRIMED_VARIABLE_PROVIDER_H_
#define PRIMED_VARIABLE_PROVIDER_H_

#include "math/vdom/variable.h"

/** Abstract class defining the interactions with "primed" variables.
 * */

class primed_variable_provider {
public:
	virtual ~primed_variable_provider() {
	}
	;

	/** Returns the ids of all variables over which the set is defined. */
	virtual const variable_id_set& get_variable_ids() const = 0;

	/** Returns the ids of the variables that are primed to degree \p prime_count. */
	virtual variable_id_set get_primed_variables(unsigned int prime_count) const {
		variable_id_set all_vars=get_variable_ids();
		variable_id_set vis;
		// Check for all indices whether their corresponding id is of prime degree prime_count
		for (variable_id_set::const_iterator it=all_vars.begin(); it!=all_vars.end(); ++it) {
			if (variable::get_prime_count(*it)==prime_count) {
				vis.insert(*it);
			}
		}
		return vis;
	}
	;

	// --------------------------------------------
	/** \name Methods changing the primedness of variables
	 *  \{ */
	// --------------------------------------------

	/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
	 */
	virtual void reassign_primedness(unsigned int d, unsigned int p = 0) = 0;

	/** Increase the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, increase all. */
	virtual void increase_primedness(unsigned int d = 0) = 0;

	/** Decrease the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, decrease all. */
	virtual void decrease_primedness(unsigned int d = 0) = 0;

	/* \} */
	// --------------------------------------------

};

#endif /*PRIMED_VARIABLE_PROVIDER_H_*/
