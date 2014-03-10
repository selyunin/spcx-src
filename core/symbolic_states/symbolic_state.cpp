#include "core/symbolic_states/symbolic_state.h"

#include "core/continuous/continuous_set.h"
#include "core/discrete/discrete_set.h"
#include "core/discrete/discrete_set_operators.h"

/** Forward declarations. */
namespace continuous {
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
continuous_set_ptr compute_intersection(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2);
continuous_set_ptr compute_or_assign_intersection(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2);
}

namespace hybrid_automata {

symbolic_state::symbolic_state() {
}

symbolic_state::symbolic_state(const discrete::discrete_set_ptr& ds,
		const continuous::continuous_set_ptr& cs) :
	d_set(ds), c_set(cs) {
}

void symbolic_state::clone(const discrete::discrete_set_ptr& ds,
		const continuous::continuous_set_ptr& cs)//:d_set(ds),c_set(cs) {};
{
	d_set = discrete::discrete_set_ptr(ds->clone());
	c_set = continuous::continuous_set_ptr(cs->clone());
}

symbolic_state::ptr symbolic_state::clone() {
	discrete::discrete_set_ptr new_d(d_set->clone());
	continuous::continuous_set_ptr new_c(c_set->clone());
	symbolic_state::ptr p = symbolic_state::ptr(new symbolic_state(new_d, new_c));
	return p;
}

bool symbolic_state::is_empty(void) const {
	return (d_set->is_empty() || c_set->is_empty()); // test discrete first, it's probably faster
}
void symbolic_state::set_continuous_set(const continuous::continuous_set::ptr& cs) {
	c_set = cs;
}
void symbolic_state::set_discrete_set(const discrete::discrete_set::ptr& ds) {
	d_set = ds;
}
const continuous::continuous_set::ptr& symbolic_state::get_continuous_set() const {
	return c_set;
}
const discrete::discrete_set::ptr& symbolic_state::get_discrete_set() const {
	return d_set;
}

void symbolic_state::intersection_assign(const symbolic_state::ptr& sstate) {
	d_set = discrete::compute_or_assign_intersection(d_set, sstate->get_discrete_set());
	if (!d_set->is_empty()) {
		c_set = compute_or_assign_intersection(c_set, sstate->get_continuous_set());
		if (c_set->is_empty()) {
			d_set=discrete::discrete_set::ptr(d_set->create_empty());
		} else {
			d_set=compute_or_assign_intersection(d_set, sstate->get_discrete_set());
		}
	} else {
		d_set=discrete::discrete_set::ptr(d_set->create_empty());
		c_set=continuous::continuous_set::ptr(c_set->create_empty());
	}
}
void symbolic_state::print(std::ostream& os) const {
	os << my_output_format.preamble;
	if (!my_output_format.skip_discrete)
		os << get_discrete_set();
	os << my_output_format.element_separator;
	if (!my_output_format.skip_continuous) {
		if (my_output_format.output_variables.empty())
		os << get_continuous_set();
		else {
			continuous::continuous_set::ptr cset(get_continuous_set()->clone());
			cset->project_to_variables(my_output_format.output_variables);
			os << cset;
		}
	}
	os << my_output_format.epilogue;
}

/* need to declare static members */
symbolic_state::output_format symbolic_state::my_output_format;

symbolic_state::~symbolic_state() {

}

}
