#ifndef LOCATION_H_
#define LOCATION_H_

#include <boost/any.hpp>
#include "core/continuous/continuous_dynamics/continuous_dynamics.h"

namespace hybrid_automata {

/** The time constraints define the continuous semantics of a location, i.e.,
 * how a state can evolve over time.
 *
 * @note They primarily consist of the invariant and dynamics (together often
 * referred to as the flow constraints), but future extensions might include
 * a time-can-pass set or flags such as "committed".
 */
class time_constraints {
public:
	time_constraints(continuous::continuous_set::const_ptr inv,
			continuous::continuous_dynamics::const_ptr dyn) :
		my_invariant(inv), my_dynamics(dyn) {
		assert(my_invariant);
		assert(my_dynamics);
	}
	;
	continuous::continuous_set::const_ptr get_invariant() const {
		return my_invariant;
	}
	;
	void set_invariant(continuous::continuous_set::const_ptr inv) {
		my_invariant = inv;
	}
	;
	continuous::continuous_dynamics::const_ptr get_dynamics() const {
		return my_dynamics;
	}
	;
	void set_dynamics(continuous::continuous_dynamics::const_ptr dyn) {
		my_dynamics = dyn;
	}
	;
protected:
	continuous::continuous_set::const_ptr my_invariant;
	continuous::continuous_dynamics::const_ptr my_dynamics;
};

/** A location class defined by dynamics and an invariant. */
class location {
public:
	typedef boost::shared_ptr<location> ptr;
	typedef boost::shared_ptr<const location> const_ptr;

	/** Construct a location with name \p loc_name, dynamics \p dyn and invariant \p inv. */
	location(std::string loc_name, time_constraints tcons) :
		my_name(loc_name), my_time_constraints(tcons) {
	}
	;
	virtual ~location();

	/** Get the time constraints. */
	virtual const time_constraints& get_time_constraints() const {
		return my_time_constraints;
	}
	;
	/** Get the time constraints. */
	virtual time_constraints& get_time_constraints() {
		return my_time_constraints;
	}
	;
	/** Set the time constraints.
	 *
	 * @note The time constraints contain pointers; the corresponding
	 * objects are adopted by *this. */
	virtual void set_time_constraints(const time_constraints& tcons) {
		my_time_constraints = tcons;
	}
	;
	virtual const std::string& get_name() const {
		return my_name;
	}
	;

	/** The following are for adding arbitrary tags */
	typedef std::string tag_key_type;
	typedef boost::any tag_value_type;
	/** Add a value with a tag */
	template<typename value_type>
	void add_tag(const tag_key_type& key, const value_type& value) {
		my_tag_map[key] = tag_value_type(value);
	};
	/** Get a tag
	 *
	 * If not present, return result of default constructor of value_type. */
	template<typename value_type>
	value_type get_tag(const tag_key_type& key) const {
		tag_map::const_iterator it = my_tag_map.find(key);
		if (it != my_tag_map.end()) {
			return boost::any_cast<value_type>(it->second);
		} else {
			return value_type();
		}
	};
	/** Check if a tag is present */
	bool has_tag(const tag_key_type& key) const {
		return my_tag_map.find(key)!=my_tag_map.end();
	};
private:
	std::string my_name;
	time_constraints my_time_constraints;
	typedef std::map<tag_key_type,tag_value_type> tag_map;
	tag_map my_tag_map;
};

}

#endif /*LOCATION_H_*/
