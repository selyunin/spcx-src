#include "core/continuous/predicate_continuous_set.h"

#include "utility/shared_ptr_output.h"

#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/predicates/valuation_function_tree_utility.h"

/** Forward declarations */
//std::ostream& operator<<(std::ostream& os, const tree::node& p);
//std::ostream& operator<<(std::ostream& os, const tree::node::ptr& p);
std::ostream& print_as_predicate(std::ostream& os, const tree::node::ptr& p);

namespace continuous {

/** Declaration of static member variables. */
predicate_continuous_set::output_format predicate_continuous_set::my_output_format =
		predicate_continuous_set::output_format();

/** Return a shared_ptr to *this.
 */
continuous_set::ptr predicate_continuous_set::get_ptr() {
	continuous_set::ptr p =
			boost::enable_shared_from_this<predicate_continuous_set>::shared_from_this();
	return p;
}

/** Return a shared_ptr to const *this.
 */
continuous_set::const_ptr predicate_continuous_set::get_const_ptr() const {
	continuous_set::const_ptr p =
			boost::enable_shared_from_this<predicate_continuous_set>::shared_from_this();
	return p;
}

/** Initialize with an element. */
predicate_continuous_set::predicate_continuous_set(predicate_type_ptr root_node) :
	my_predicate(root_node) {
}

predicate_continuous_set::~predicate_continuous_set() {
}

predicate_continuous_set::predicate_type_ptr predicate_continuous_set::get_predicate() const {
	return my_predicate;
}

void predicate_continuous_set::set_predicate(predicate_type_ptr new_pred) {
	my_predicate = new_pred;
}

predicate_continuous_set* predicate_continuous_set::create_universe() const {
	/** Create boolean node with value true. */
	predicate_type_ptr new_root = predicate_type_ptr(
			new valuation_functions::const_node<bool>(true));
	return new predicate_continuous_set(new_root);
}

predicate_continuous_set* predicate_continuous_set::create_empty() const {
	/** Create boolean node with value false. */
	predicate_type_ptr new_root = predicate_type_ptr(
			new valuation_functions::const_node<bool>(false));
	return new predicate_continuous_set(new_root);
}

predicate_continuous_set* predicate_continuous_set::clone() const {
	return new predicate_continuous_set(my_predicate);
}

int predicate_continuous_set::get_memory() const {
	throw std::runtime_error("missing implementation");
	return 0;
}

const variable_id_set& predicate_continuous_set::get_variable_ids() const {
	static variable_id_set vis;
	vis = valuation_functions::get_variable_ids(my_predicate);
	return vis;
}

unsigned int predicate_continuous_set::get_dim() const {
	throw std::runtime_error("missing implementation");
	return 0;
}

/** The collection is by construction empty iff the first element is empty. */
math::tribool predicate_continuous_set::is_empty() const {
	//throw std::runtime_error("missing implementation");
	return math::indeterminate();
}

void predicate_continuous_set::embed_variables(const variable_id_set& id_set) {
	throw std::runtime_error("missing implementation");
}

void predicate_continuous_set::existentially_quantify_variables(const variable_id_set& id_set) {
	throw std::runtime_error("missing implementation");
}

/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
 */
void predicate_continuous_set::reassign_primedness(unsigned int d, unsigned int p) {
	throw std::runtime_error("missing implementation");
}

/** Increase the primedness of the variables with primedness of degree \p d by 1.
 * If d is 0, increase all. */
void predicate_continuous_set::increase_primedness(unsigned int d) {
	throw std::runtime_error("missing implementation");
}

/** Decrease the primedness of the variables with primedness of degree \p d by 1.
 * If d is 0, decrease all. */
void predicate_continuous_set::decrease_primedness(unsigned int d) {
	throw std::runtime_error("missing implementation");
}

void predicate_continuous_set::intersection_assign(const predicate_continuous_set& cset) {
	/** Create boolean AND node from my_predicate and cset.my_predicate. */
	if (my_predicate && cset.my_predicate) {
		predicate_type_ptr new_root = predicate_type_ptr(
				new valuation_functions::boolean_node(AND, my_predicate,
						cset.my_predicate));
	my_predicate=new_root;
	}
	else if(cset.my_predicate)
		my_predicate = cset.my_predicate;

}

/** Accept a visitor. */
void predicate_continuous_set::accept(dispatching::dispatcher<continuous_set_typelist>& d) const {
	d.dispatch(this);
}

/** Output as a stream of characters. */
void predicate_continuous_set::print(std::ostream& os) const {
	output_format of = get_output_format();
	os << of.preamble;

	// GF: which printer to use?
	//throw std::runtime_error("missing implementation");
	if (my_predicate)
		print_as_predicate(os, my_predicate);
	else
		os << of.nullstring;

	os << of.epilogue;
}

predicate_continuous_set::predicate_continuous_set() {
}

}

