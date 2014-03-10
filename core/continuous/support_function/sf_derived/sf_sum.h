/*
 * sf_sum.h
 *
 *  Created on: Apr 2, 2010
 *      Author: notroot
 */

#ifndef SF_SUM_H_
#define SF_SUM_H_

#include "core/continuous/support_function/sf_base/sf_binary.h"


namespace continuous {
namespace support_function {

/** An implicit representation of the Minkowski sum of two continuous sets by
 * means of their support function.
 */

template<typename scalar_type> class sf_sum: public sf_binary<scalar_type> {
public:
	typedef typename sf_binary<scalar_type>::affine_map affine_map;

	/** Construct a support function set representation of the minkowski sum of
	 * two continuous sets. */
	sf_sum(const support_function_provider::const_ptr& s1,
			const support_function_provider::const_ptr& s2) :
		sf_binary<scalar_type> (s1, s2) {
	}
	;

	/** Construct a support function set representation of the minkowski sum of
	 * two continuous sets, transformed by an affine map M. */
	sf_sum(const support_function_provider::const_ptr& s1,
			const support_function_provider::const_ptr& s2, const affine_map& M) :
		sf_binary<scalar_type> (s1, s2, M) {
	}
	;

	virtual ~sf_sum() {
	}
	;

	/** Shallow copy */
	virtual sf_sum<scalar_type>* clone() const {
		sf_sum<scalar_type>* p;
		if (this->get_map()) {
			p = new sf_sum<scalar_type> (this->my_set1, this->my_set2,
					*(this->get_map()));
		} else {
			p = new sf_sum<scalar_type> (this->my_set1, this->my_set2);
		}
		return p;
	}
	;

	/** Return whether the set is empty */
	virtual math::tribool is_empty() const {
		return (this->my_set1->is_empty()) || (this->my_set2->is_empty());
	}
	;

	//-------------------------------------------------------
	// support_function_provider functions
	//-------------------------------------------------------
	/** Returns true if compute_support returns a support vector
	 * and false otherwise.
	 */
	virtual bool computes_support_vector() const {
		return this->my_set1->computes_support_vector() && this->my_set2->computes_support_vector();
	}
	;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<Rational>& l,
			Rational& max_value, math::vdom_vector<Rational>& support_vec,
			bool& is_empty, bool& is_bounded) const {
		Rational max_value1;
		math::vdom_vector<Rational> support_vec1;
		bool is_empty1;
		bool is_bounded1;
		// 	compute_support_mapped<Rational> (my_set, my_map, l, max_value, support_vec, is_empty, is_bounded);
		sf_set<scalar_type>::template compute_support_mapped<Rational>(
				*this->my_set1, this->get_map(), l, max_value1, support_vec1,
				is_empty1, is_bounded1);
		if (is_bounded1) {
			Rational max_value2;
			math::vdom_vector<Rational> support_vec2;
			bool is_empty2;
			bool is_bounded2;
			sf_set<scalar_type>::template compute_support_mapped<Rational>(
					*this->my_set2, this->get_map(), l, max_value2, support_vec2,
					is_empty2, is_bounded2);
			if (is_bounded2) {
				max_value = max_value1 + max_value2;
				if (support_vec1.size() > 0 && support_vec2.size() > 0) {
					support_vec = support_vec1 + support_vec2;
				}
				is_empty = is_empty1 && is_empty2;
				is_bounded = true;
			} else {
				// set1 is unbounded
				is_empty = is_empty1 && is_empty2;
				is_bounded = false;
			}
		} else {
			// set1 is unbounded
			is_empty = is_empty1;
			is_bounded = false;
		}
	}
	;

	/** Computes the support function for a given vector l, i.e.,
	 * the supremum of l*x over all x in *this.
	 * If *this is empty, then max_value=0.
	 * If l has a variable that is not constrained in *this, then the result is
	 * considered unbounded, so is_bounded=false.
	 */
	virtual void compute_support(const math::vdom_vector<double>& l,
			double& max_value, math::vdom_vector<double>& support_vec,
			bool& is_empty, bool& is_bounded) const {
		double max_value1;
		math::vdom_vector<double> support_vec1;
		bool is_empty1;
		bool is_bounded1;
		// 	compute_support_mapped<double> (my_set, my_map, l, max_value, support_vec, is_empty, is_bounded);
		sf_set<scalar_type>::template compute_support_mapped<double>(
				*this->my_set1, this->get_map(), l, max_value1, support_vec1,
				is_empty1, is_bounded1);
		if (is_bounded1) {
			double max_value2;
			math::vdom_vector<double> support_vec2;
			bool is_empty2;
			bool is_bounded2;
			sf_set<scalar_type>::template compute_support_mapped<double>(
					*this->my_set2, this->get_map(), l, max_value2, support_vec2,
					is_empty2, is_bounded2);
			if (is_bounded2) {
				max_value = max_value1 + max_value2;
				if (support_vec1.size() > 0 && support_vec2.size() > 0) {
					support_vec = support_vec1 + support_vec2;
				}
				is_empty = is_empty1 && is_empty2;
				is_bounded = true;
			} else {
				// set1 is unbounded
				is_empty = is_empty1 && is_empty2;
				is_bounded = false;
			}
		} else {
			// set1 is unbounded
			is_empty = is_empty1;
			is_bounded = false;
		}
	}
	;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const {
		os << "{";
		if (this->get_map()) {
			os << this->get_map() << " * ";
		}
		os << "(" << this->my_set1 << "+" << this->my_set2 << ")";
		os << "}";
	}
};

}
}

#endif /* SF_SUM_H_ */
