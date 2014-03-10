#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC.h"

#include "core/continuous/polyhedra/ppl_polyhedron/fp_interface.h"
//#include <ppl.hh>
#include "core/continuous/polyhedra/ppl_polyhedron/rat_linexpression.h"
#include "core/continuous/polyhedra/ppl_polyhedron/convert_matrix_math.h"
//#include "../../math/lin_constraint_system.h"

//#include "convert_to_ppl.h" // I couldn't get this to compile without including this here (circular declaration)
//#include "../../abstract_framework/continuous/continuous_set_operators.h"
#include "core/continuous/continuous_set_operators.h"
#include "core/continuous/polyhedra/polyhedron_output.h"

#include "core/continuous/continuous_set_transforms/constant_bound_time_elapse_transform.h"
#include "core/continuous/continuous_set_transforms/reset_function_transform.h"

/** Foward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
typedef boost::shared_ptr<const node> node_const_ptr;
}
namespace ppl_polyhedron {
Rational_Linear_Expression convert_to_Rational_Linear_Expression(const tree::node_const_ptr& p,
		const index_to_variable_id_map_ptr& iimap);
}

namespace ppl_polyhedron {

using namespace std;
using namespace Parma_Polyhedra_Library;
using namespace Parma_Polyhedra_Library::IO_Operators;
using namespace continuous;

//template<> continuous::polyhedron<Rational>::output_format
//		continuous::polyhedron<Rational>::my_output_format = TEXTUAL;

continuous_set_PPL_NNC::continuous_set_PPL_NNC() {
}

continuous_set_PPL_NNC::~continuous_set_PPL_NNC() {
}

continuous_set_PPL_NNC::continuous_set_PPL_NNC(const index_to_variable_id_map_ptr& pnew_map):index_to_variable_id_map_provider(pnew_map) {
	// Create a universe polyhedron of the right dimensions
	_mypoly = NNC_Polyhedron(pnew_map->dimensions());
}

continuous_set_PPL_NNC::continuous_set_PPL_NNC(const my_poly_type& poly,
		const index_to_variable_id_map_ptr& pnew_map):index_to_variable_id_map_provider(pnew_map),
	_mypoly(poly) {
}

continuous_set_PPL_NNC* continuous_set_PPL_NNC::clone() const {
	return new continuous_set_PPL_NNC(_mypoly, get_index_to_variable_id_map());
}

continuous_set_PPL_NNC* continuous_set_PPL_NNC::create_universe(
		const index_to_variable_id_map_ptr pnew_map) const {
	continuous_set_PPL_NNC* p = new continuous_set_PPL_NNC(NNC_Polyhedron(
			get_index_to_variable_id_map()->dimensions()), pnew_map);
	return p;
}

continuous_set_PPL_NNC* continuous_set_PPL_NNC::create_empty(const index_to_variable_id_map_ptr pnew_map) const {
	continuous_set_PPL_NNC* p = new continuous_set_PPL_NNC(NNC_Polyhedron(get_index_to_variable_id_map()->dimensions(),
					Parma_Polyhedra_Library::EMPTY), pnew_map);
	return p;
}

continuous_set_PPL_NNC* continuous_set_PPL_NNC::create_universe() const {
	return create_universe(get_index_to_variable_id_map());
}

continuous_set_PPL_NNC* continuous_set_PPL_NNC::create_empty() const {
	return create_empty(get_index_to_variable_id_map());
}

int continuous_set_PPL_NNC::get_memory() const {
	return (int) _mypoly.total_memory_in_bytes();
}

//void continuous_set_PPL_NNC::print_textual(std::ostream& os) const {
//	// construct vnvec
//	variable_id_map vnvec;
//	for (index_type i = 0; i < get_dim(); ++i) {
//		vnvec.insert(i, variable(get_id(i)));
//		//std::cout << vnvec << std::flush << std::endl;
//	}
//
//	if (_mypoly.is_empty()) {
//		os << "false";
//	} else if (_mypoly.is_universe()) {
//		os << "true";
//	} else {
//		Constraint_System cs;
//		cs = _mypoly.constraints();
//		Constraint_System::const_iterator i;
//		for (i = cs.begin(); i != cs.end(); i++) {
//			if (i != cs.begin())
//				os << " & ";
//
//			ppl_polyhedron::print_constraint(os, *i, vnvec);
//		}
//	}
//}
//
//void continuous_set_PPL_NNC::print_double_constraints(std::ostream& os) const {
//	double_point_list pl;
//	add_ccvs_consys_to_double_point_list(_mypoly, pl);
//	print_fp_raw(os, pl);
//}
//
//void continuous_set_PPL_NNC::print_double_generators(std::ostream& os) const {
//	double_point_list pl;
//	add_ccvs_to_double_point_list(_mypoly, pl);
//	print_fp_raw(os, pl);
//}

dimension_t continuous_set_PPL_NNC::get_dim() const {
	return _mypoly.space_dimension();
}

math::tribool continuous_set_PPL_NNC::is_empty() const {
	return _mypoly.is_empty();
}

math::tribool continuous_set_PPL_NNC::is_universe() const {
	return _mypoly.is_universe();
}

bool continuous_set_PPL_NNC::computes_support_vector() const {
	return true;
}
;

void continuous_set_PPL_NNC::compute_support(const math::vdom_vector<Rational>& l,
		Rational& max_value, math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	Integer denominator;
	Linear_Expression le;
	bool has_other_variables;
	has_other_variables = convert_to_Linear_Expression(math::lin_expression<Rational>(l), le, denominator,
			get_index_to_variable_id_map());

	if (has_other_variables) {
		is_empty = this->is_empty();
		if (!is_empty) {
			is_bounded = false;
		}
	} else {
		//bool 	maximize (const Linear_Expression &expr, Coefficient &sup_n, Coefficient &sup_d, bool &maximum, Generator &point) const
		Coefficient sup_n, sup_d;
		bool is_maximum;
		Generator point = Generator::zero_dim_point();
		bool notempty_and_bounded_from_above =
				_mypoly.maximize(le, sup_n, sup_d, is_maximum, point);

		if (notempty_and_bounded_from_above) {
			is_empty = false;
			is_bounded = true;
			max_value = Rational(sup_n, sup_d * denominator);
			support_vec = convert_to_vdom_vector(point, get_index_to_variable_id_map());
		} else {
			if (this->is_empty()) {
				is_empty = true;
				is_bounded = true;
			} else {
				is_empty = false;
				is_bounded = false;
			}
		}
	}
}

void continuous_set_PPL_NNC::compute_support(const math::vdom_vector<double>& l,
		double& max_value, math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {

	Rational m;
	math::vdom_vector<Rational> lrat;
	math::vdom_vector<Rational> srat;
	lrat = l.convert_to<Rational> ();

	compute_support(lrat, m, srat, is_empty, is_bounded);

	support_vec = srat.convert_to<double> ();
	max_value = m.get_double();
}

continuous_set_PPL_NNC::const_ptr continuous_set_PPL_NNC::map_to_common_space_const(
		const continuous::continuous_set_const_ptr& p) {
	/* In this version the returned poly is a const versions q (only clones if necessary) */
	const_ptr mapped_p = boost::dynamic_pointer_cast<const/* automatic formatting fuses this with the next line */
	continuous_set_PPL_NNC, const continuous_set>(p);
	/* map both to a common space */
	if (get_index_to_variable_id_map() != mapped_p->get_index_to_variable_id_map()) {
		mapped_p = map_to_common_space(mapped_p);
	}
	return mapped_p;
}

//void continuous_set_PPL_NNC::get_common_space_versions(const continuous_set_const_ptr p,
//		const continuous_set_const_ptr q, ptr& mapped_p, ptr& mapped_q) {
//	/* In this version the returned poly
//	 * are non-const clones of p and q */
//	mapped_p = ptr(p->clone());
//	mapped_q = ptr(q->clone());
//	/* map both to a common space */
//	if (mapped_p->get_index_to_variable_id_map() != mapped_q->get_index_to_variable_id_map()) {
//		mapped_q = mapped_p->map_to_common_space(mapped_q);
//	}
//}

void continuous_set_PPL_NNC::get_common_space_versions_both_const(const continuous_set_const_ptr p,
		const continuous_set_const_ptr q, const_ptr& mapped_p, const_ptr& mapped_q) {
	/* In this version the returned poly are
	 * const versions of p and q (only clones if necessary) */
	mapped_p = boost::static_pointer_cast<const continuous_set_PPL_NNC, const continuous_set>(p);
	mapped_q = boost::static_pointer_cast<const continuous_set_PPL_NNC, const continuous_set>(q);
	/* map both to a common space */
	if (mapped_p->get_index_to_variable_id_map() != mapped_q->get_index_to_variable_id_map()) {
		ptr nonconst_p(mapped_p->clone());
		mapped_q = nonconst_p->map_to_common_space(mapped_q);
		mapped_p = nonconst_p;
	}
}

math::tribool continuous_set_PPL_NNC::is_disjoint_from(const continuous_set_const_ptr& p) const {
	const_ptr ppl_of_mapped_myself;
	const_ptr ppl_of_mapped_p;
	get_common_space_versions_both_const(p, get_const_ptr(), ppl_of_mapped_p, ppl_of_mapped_myself);

	return ppl_of_mapped_myself->_mypoly.is_disjoint_from(ppl_of_mapped_p->get_poly());
}

math::tribool continuous_set_PPL_NNC::contains(const continuous_set_const_ptr& p) const {
	const_ptr ppl_of_mapped_myself;
	const_ptr ppl_of_mapped_p;
	get_common_space_versions_both_const(p, get_const_ptr(), ppl_of_mapped_p, ppl_of_mapped_myself);
	//std::cout << ppl_of_mapped_p->get_poly().space_dimension() << "!" << ppl_of_mapped_p->get_index_to_variable_id_map() << std::endl;
	//std::cout << ppl_of_mapped_myself->get_poly().space_dimension() << "?" << ppl_of_mapped_myself->get_index_to_variable_id_map() << std::endl;
	return ppl_of_mapped_myself->_mypoly.contains(ppl_of_mapped_p->get_poly());
}

void continuous_set_PPL_NNC::embed_variables(const variable_id_set& id_set) {
	// add the ids to _index_to_variable_id_map_ptr
	index_to_variable_id_map_ptr p = get_index_to_variable_id_map()->get_map_with_ids_added(id_set);
	set_index_to_variable_id_map(p);
	_mypoly.add_space_dimensions_and_embed(p->dimensions() - get_dim());
}

void continuous_set_PPL_NNC::existentially_quantify_variables(const variable_id_set& id_set) {
	Variables_Set vs;
	index_type index;
	bool found;
	for (variable_id_set::const_iterator it = id_set.begin(); it != id_set.end(); ++it) {
		// Test if the variable is in the domain of *this
		index = get_index_to_variable_id_map()->check_for_index(*it, found);
		if (found) {
			vs.insert(Variable(index));
		}
	}
	_mypoly.remove_space_dimensions(vs);

	// fix _index_to_variable_id_map_ptr
	index_to_variable_id_map_ptr p = get_index_to_variable_id_map()->get_map_with_ids_removed(
			id_set);
	set_index_to_variable_id_map(p);
}

void continuous_set_PPL_NNC::simplify() {
	_mypoly.minimized_constraints();
}

void continuous_set_PPL_NNC::swap(ptr s_ppl) {
	// Swap the _index_to_variable_id_map_ptr
	index_to_variable_id_map_ptr p1 = get_index_to_variable_id_map();
	index_to_variable_id_map_ptr p2 = s_ppl->get_index_to_variable_id_map(); // Note: must use pointer to continuous_set_PPL_NNC because access to protected must pass via derived class
	set_index_to_variable_id_map(p2);
	s_ppl->set_index_to_variable_id_map(p1);

	// Swap the poly
	_mypoly.swap(s_ppl->_mypoly);
}

index_to_index_bimap continuous_set_PPL_NNC::map_to_common_iimap(
		const index_to_variable_id_map_ptr& p) {
	// map both to a common space
	dimension_t newdim;
	index_to_variable_id_map_ptr new_map = index_to_variable_id_map_ptr(
			new index_to_variable_id_map);
	index_to_index_bimap iimap;
	get_common_map(get_index_to_variable_id_map(), p, new_map, newdim, iimap);

	// Add missing dimensions to *this
	_mypoly.add_space_dimensions_and_embed(newdim - _mypoly.space_dimension());
	set_index_to_variable_id_map(new_map);

	//cout << "iimapped:" << _mypoly.space_dimension() << ":" << get_index_to_variable_id_map() << endl;

	return iimap;
}

continuous_set_PPL_NNC::ptr continuous_set_PPL_NNC::map_to_common_space(const const_ptr& p) {
	// map both to a common space
	index_to_index_bimap imap;
	imap = map_to_common_iimap(p->get_index_to_variable_id_map());

	// Create a mapped copy of p
	ptr ret_p(p->clone());
	ret_p->_mypoly.add_space_dimensions_and_embed(_mypoly.space_dimension()
			- ret_p->_mypoly.space_dimension());
	ret_p->_mypoly.map_space_dimensions(imap);
	ret_p->set_index_to_variable_id_map(get_index_to_variable_id_map());
	//cout << "Q:" << iq->_mypoly << flush << ":";
	//iq->print(cout);
	//cout << "commapped:" << ret_p->_mypoly.space_dimension() << ":" << ret_p->get_index_to_variable_id_map() << endl;

	return ret_p;
	//return continuous_set_ptr(iq);
}

void continuous_set_PPL_NNC::intersection_assign(const continuous_set_const_ptr p) {
	const_ptr ppl_of_mapped_p = map_to_common_space_const(p);
	_mypoly.intersection_assign(ppl_of_mapped_p->_mypoly);
}

void continuous_set_PPL_NNC::intersection_assign(const continuous_set_PPL_NNC& p) {
	// to do : work with references ?
	continuous_set_ptr ptr(p.clone()); // this unnecessary clone is just to get a pointer
	const_ptr ppl_of_mapped_p = map_to_common_space_const(ptr);
	_mypoly.intersection_assign(ppl_of_mapped_p->_mypoly);
}

void continuous_set_PPL_NNC::union_assign(const continuous_set_const_ptr p) {
	const_ptr ppl_of_mapped_p = map_to_common_space_const(p);
	_mypoly.poly_hull_assign(ppl_of_mapped_p->_mypoly);
}

void continuous_set_PPL_NNC::difference_assign(const continuous_set_const_ptr p) {
	const_ptr ppl_of_mapped_p = map_to_common_space_const(p);
	_mypoly.poly_difference_assign(ppl_of_mapped_p->_mypoly);
}

void continuous_set_PPL_NNC::cheap_difference_assign(const continuous_set_const_ptr p) {
	const_ptr ppl_of_mapped_p = map_to_common_space_const(p);
	if (ppl_of_mapped_p->_mypoly.contains(_mypoly)) {
		// Old version: keep the index_to_variable_id_map
		//		ip->_mypoly=NNC_Polyhedron(ip->get_index_to_variable_id_map()->dimensions(), Parma_Polyhedra_Library::EMPTY);

		// create an empty set
		ptr empty_set(create_empty());
		swap(empty_set);
	}
	// otherwise leave *this unchanged
}

/*! Adds the constraint \p c to \p *this.
 */
void continuous_set_PPL_NNC::add_constraint(const math::lin_constraint<Rational> &c, bool check_redundancy) {
	// Convert to linear expression and add the corresponding constraint.
	// Attention: The semantics of the sign is opposite to the PPL!
	Integer denominator;
	Linear_Expression l;

//	std::cerr << "adding " << c << std::endl;
//	std::cerr << this->get_index_to_variable_id_map();

	if (c.get_l().get_index_to_variable_id_map() != get_index_to_variable_id_map()) {
		// If there is a variable in c that is not yet defined in the polyhedron, remap the polyhedron
		map_to_common_iimap(c.get_l().get_index_to_variable_id_map());
	}
//	std::cerr << "remapped:" << this->get_index_to_variable_id_map();

	convert_to_Linear_Expression(c.get_l(), l, denominator, this->get_index_to_variable_id_map());
	if (c.get_sign() == LT) // lin_constraint<Rational>::sign
		add_constraint(l < Integer(0));
	else if (c.get_sign() == LE) // lin_constraint<Rational>::sign
		add_constraint(l <= Integer(0));
	else if (c.get_sign() == EQ) // lin_constraint<Rational>::sign
		add_constraint(l == Integer(0));
	else if (c.get_sign() == GE) // lin_constraint<Rational>::sign
		add_constraint(l >= Integer(0));
	else if (c.get_sign() == GT) // lin_constraint<Rational>::sign
		add_constraint(l > Integer(0));
}

void continuous_set_PPL_NNC::add_constraint(
		const Parma_Polyhedra_Library::Constraint &c, bool check_redundancy) {
	if (check_redundancy) {
		// deprecated in PPL 0.10
//		_mypoly.add_constraint_and_minimize(c);

		// only add if the constraint is non-redundant, i.e.,
		// the poly is not included in it
		if (!(_mypoly.relation_with(c)==Parma_Polyhedra_Library::Poly_Con_Relation::is_included()))
			_mypoly.add_constraint(c);
	} else {
		_mypoly.add_constraint(c);
	}
}

void continuous_set_PPL_NNC::remove_redundant_constraints() {
	_mypoly.minimized_constraints();
}

const Parma_Polyhedra_Library::Constraint_System& continuous_set_PPL_NNC::constraints() const {
	return _mypoly.constraints();
}

/** Compute constant bound time elapse. */
void continuous_set_PPL_NNC::assign_transformation(
		const constant_bound_time_elapse_transform& t) {
	// Bring the two sets onto the same variable space
	const_ptr time_elapse_set = map_to_common_space_const(t.get_set());

	// Apply the actual transform
	_mypoly.time_elapse_assign(time_elapse_set->get_poly());
}

void continuous_set_PPL_NNC::assign_transformation(
		const reset_function_transform& t) {
	// Get the id of the variable that's being assigned to.
	variable_id vid = t.get_variable_id();

	// Make sure the variable exists in *this
	variable_id_set vis;
	vis.insert(vid);
	embed_variables(vis);

	// get the linear expression corresponding to the function
	Rational_Linear_Expression e = convert_to_Rational_Linear_Expression(t.get_function_const(),
			get_index_to_variable_id_map());
	Linear_Expression l = e.get_LE();
	// Make sure it's the same dimension as my_poly
	if (get_dim() > 0)
		l += 0 * Variable(get_dim() - 1);

	//std::cout << e << std::endl;
	//std::cout << *this << std::endl;

	// Apply the actual transform
	_mypoly.affine_image(Variable(get_index(vid)), l, e.get_denominator());

	//std::cout << *this << std::endl;
}

/*	if (reset_value_transform* pt = dynamic_cast<reset_value_transform*>(t.get())) {
 // Get the id of the variable that's being assigned to.
 variable_id vid=pt->get_variable_id();

 // Make sure the variable exists in *this
 variable_id_set vis;
 vis.insert(vid);
 embed_variables(vis);

 // Apply the reset
 const Rational& r=pt->get_value<Rational>();
 _mypoly.affine_image(Variable(get_index(vid)),Linear_Expression(r.get_num()),r.get_den());
 } else */

/*	} else if (ode_time_elapse_transform* pt = dynamic_cast<ode_time_elapse_transform*>(t.get())) {
 // Construct a predicate tree that encodes the constraints
 //    x_i' == f_i
 if (pt->begin()!=pt->end())
 {
 tree::node::ptr p_func,p_var,p_comp,p_root;
 for (ode_time_elapse_transform::const_iterator it=pt->begin();it!=pt->end();++it) // otherwise continue with the second assignment
 {
 p_func=it->second; // get the function of the pair<variable_id,function_ptr>
 p_var=tree::node::ptr(new variable_node(variable::get_primed_id(it->first))); // get the variable id of the pair<variable_id,function_ptr>
 p_comp=tree::node::ptr(new comparison_node(EQ,p_var,p_func));
 if (it==pt->begin())
 {
 if (reset_function_transform* pt = dynamic_cast<reset_function_transform*>(t.get())) {

 } else if (constant_bound_time_elapse_transform* pt = dynamic_cast<constant_bound_time_elapse_transform*>(t.get())) {
 // Bring the two sets onto the same variable space
 const_ptr time_elapse_set=map_to_common_space_const(pt->get_set());

 // Apply the actual transform
 _mypoly.time_elapse_assign(time_elapse_set->get_poly());
 p_root=p_comp;
 } else
 {
 p_root=tree::node::ptr(new boolean_node(AND,p_root,p_comp)); // keep connecting the EQ nodes with AND nodes
 }
 }

 // Construct a set of constraints from the functions
 // This is then the set { x_1'==f_1 & x_2'==f_2 ... }.
 PPL_NNC_generator<Rational> gen;
 continuous_set_ptr cs=gen.convert(p_root);
 //cout << "Derivative constraints:" << cs << endl;
 // Intersect with the invariant
 cs=compute_or_assign_intersection(cs,pt->get_set());
 //cout << "with invariant:" << cs << endl;

 // Existentially quantify over x (unprimed)
 cs->existentially_quantify_variables(cs->get_primed_variables(0));
 //cout << "after quantify:" << cs << endl;
 // Go from x' back to x
 cs->decrease_primedness();
 //cout << "primedness_decreased:" << cs << endl;

 // Bring the two sets onto the same variable space
 const_ptr time_elapse_set=map_to_common_space_const(cs);

 // Apply the actual transform
 _mypoly.time_elapse_assign(time_elapse_set->get_poly());
 } */

void continuous_set_PPL_NNC::accept(dispatching::dispatcher<continuous_set_typelist>& d) const {
	d.dispatch(this);
}

math::lin_constraint_system<Rational>::const_ptr continuous_set_PPL_NNC::get_constraints() const {
	math::lin_constraint_system<Rational>::ptr p = math::lin_constraint_system<Rational>::ptr(
			new math::lin_constraint_system<Rational> ());
	const Parma_Polyhedra_Library::Constraint_System& cs = constraints();
	for (Parma_Polyhedra_Library::Constraint_System::const_iterator it = cs.begin(); it != cs.end(); ++it) {
		math::lin_constraint<Rational> c = convert_to_lin_constraint(*it,
				this->get_index_to_variable_id_map());
		p->insert(c);
	}
	//return boost::static_pointer_cast<const lin_constraint_system<Rational> >(p);
	return p;
}

}

