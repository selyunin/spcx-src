/*
 * polyhedron_collection.h
 *
 *  Created on: Oct 18, 2010
 *      Author: frehse
 */

#ifndef POLYHEDRON_COLLECTION_H_
#define POLYHEDRON_COLLECTION_H_

#include <list>
#include "polyhedron.h"

namespace continuous {

/**
 * A class for representing a collection of polyhedra.
 *
 * The semantics is in the sense of a set of linear
 * constraints in DNF. Consequently, the empty
 * collection presents the universe, while
 * an empty set of states is represented by
 * an empty polyhedron.
 *
 * However, the class needs to know what implementation
 * to use for polyhedra, e.g., when constraints are
 * added to a universe collection.
 * It is constructed with at least one element,
 * and keept at least one element in the collection
 * for this purpose.
 */
template<typename scalar_type>
class polyhedron_collection: public continuous_set {
public:
	typedef boost::shared_ptr<polyhedron_collection> ptr;
	typedef boost::shared_ptr<const polyhedron_collection> const_ptr;

	typedef typename polyhedron<scalar_type>::ptr element_type;
	typedef std::list<element_type> container_type;
	typedef typename container_type::iterator iterator;
	typedef typename container_type::const_iterator const_iterator;

	/** Intermediate constructor -- not in a consistent state. */
	polyhedron_collection();

	/** Initialize with adopting an element. */
	explicit polyhedron_collection(typename polyhedron<scalar_type>::ptr element);

	/** Destructor */
	virtual ~polyhedron_collection();

	/** Create a deep copy. */
	virtual polyhedron_collection* clone() const;

	virtual polyhedron_collection* create_universe() const;
	virtual polyhedron_collection* create_empty() const;

	/** Return a representative element that can be used to
	 * create new elements via virtual constructors.
	 *
	 * The default element can be used to create
	 * polyhedra of the same type as (at least one) in
	 * the collection.
	 */
	virtual typename polyhedron<scalar_type>::const_ptr default_element() const;

	virtual iterator begin();
	virtual iterator end();
	virtual const_iterator begin() const;
	virtual const_iterator end() const;

	virtual unsigned int size() const;

	virtual int get_memory() const;
	virtual const variable_id_set& get_variable_ids() const;

	virtual unsigned int get_dim() const;

	virtual math::tribool is_empty() const;

	virtual void embed_variables(const variable_id_set& id_set);
	virtual void existentially_quantify_variables(const variable_id_set& id_set);

	/** Returns whether there is an element that contains p */
	virtual math::tribool
	element_wise_contains(typename polyhedron<scalar_type>::const_ptr p) const;

	/** Returns whether all elements of polys are contained in some element */
	virtual math::tribool
	element_wise_contains(const polyhedron_collection<scalar_type>& polys) const;

	/** Delete the elements that are contained in polyhedron p */
	virtual void delete_if_contained_in(
			typename polyhedron<scalar_type>::const_ptr p);

	/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
	 */
	virtual void reassign_primedness(unsigned int d, unsigned int p = 0);

	/** Increase the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, increase all. */
	virtual void increase_primedness(unsigned int d = 0);

	/** Decrease the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, decrease all. */
	virtual void decrease_primedness(unsigned int d = 0);

	/** Insert p into the container.
	 *
	 * With redundancy_check=true,
	 * the element is only inserted if it is not contained already; also,
	 * all elements contained in p are removed. */
	virtual void insert(typename polyhedron<scalar_type>::ptr p,
			bool redundancy_check = false);

	/** Insert the elements of a polyhedron_collection into the container, without redundancy check.
	 *
	 * The elements of polys are inserted without cloning (adopted).
	 * Use with care and clone before if in doubt.
	 */
	virtual void insert(const polyhedron_collection<scalar_type>& polys);

	/*! Adds the constraint \p c to \p *this.
	 *
	 * Does not remove any empty sets in the collection.
	 */
	virtual void add_constraint(const math::lin_constraint<scalar_type> &c, bool check_redundancy = false);

	/** Add the constraints in the set con_set.
	 *
	 * Removes empty sets in the collection until the collection
	 * is of size 1.
	 */
	virtual void add_constraints(
			const math::lin_constraint_system<scalar_type>& con_set, bool check_redundancy = false);

	/** Remove all redundant constraints.
	 * */
	virtual void remove_redundant_constraints();

	/** Returns true if compute_support returns a support vector
	 * and false otherwise.
	 *
	 * The statement must hold for the current state of the object,
	 * and remain until the object is modified.
	 */
	virtual bool computes_support_vector() const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<Rational>& l,
			Rational& max_value, math::vdom_vector<Rational>& support_vec,
			bool& is_empty, bool& is_bounded) const;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<double>& l,
			double& max_value, math::vdom_vector<double>& support_vec,
			bool& is_empty, bool& is_bounded) const;

	/** Returns the set in predicate form.*/
	virtual continuous_set_predicate::ptr get_predicate() const;

	/** Accept a visitor. */
	virtual void
	accept(dispatching::dispatcher<continuous_set_typelist>& d) const;

	class output_format {
	public:
		output_format() { // default values
			preamble = "(";
			epilogue = ")";
			element_separator = " | ";
		}
		;
		static output_format matlab() {
			output_format f;
			f.preamble = "[";
			f.epilogue = "]";
			f.element_separator = ",";
			return f;
		}
		;
		static output_format space_separated() {
			output_format f;
			f.preamble = "";
			f.epilogue = "";
			f.element_separator = " ";
			return f;
		}
		;
		std::string preamble; // written before matrix
		std::string epilogue; // written after matrix
		std::string element_separator; // written between elements
	};
	static const output_format& get_output_format() {
		return my_output_format;
	}
	;
	static void set_output_format(output_format new_format) {
		my_output_format = new_format;
	}
	;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const;

	/** Remove empty elements in *this.
	 *
	 * Leaves one empty element if there are no others
	 * so that the container has a default element. */
	virtual void remove_empty();

	/** Swap contents */
	virtual void swap(polyhedron_collection<scalar_type>& polys);

protected:
	/** Insert element at end of list without emptiness checking. */
	virtual void push_back(typename polyhedron<scalar_type>::ptr p);

	/** Remove empty elements in *this in the range [beg,end).
	 *
	 * Leaves one empty element if there are no others
	 * so that the container has a default element. */
	virtual void remove_empty(const iterator& ibeg,const iterator& iend);

private:
	container_type my_container;
	static output_format my_output_format;
};

template<typename scalar_type> typename polyhedron_collection<scalar_type>::output_format
		polyhedron_collection<scalar_type>::my_output_format;

}

#include "polyhedron_collection.hpp"

#endif /* POLYHEDRON_COLLECTION_H_ */
