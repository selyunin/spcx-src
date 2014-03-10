#ifndef CONTINUOUS_SET_OPERATORS_H_
#define CONTINUOUS_SET_OPERATORS_H_

#include "boost/shared_ptr.hpp"
//#include "continuous_set.h"
#include "math/tribool.h"

/** This file declares the (binary) operators on continuous_set.
 * The actual implementation is given in continuous_set_operator_implementations. */

/** Forward declarations */
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;
class continuous_set_transform;
typedef boost::shared_ptr<const continuous_set_transform>
		continuous_set_transform_const_ptr;
}

namespace continuous {

/** Returns true iff p1 contains p2. */
math::tribool containment_test(const continuous_set* const p1,
		const continuous_set* const p2);

/** Returns true iff p1 contains p2. */
math::tribool containment_test(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2);

// Operators that possibly take as arguments two continuous_sets of different derived classes
continuous_set_ptr compute_intersection(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2);

continuous_set_ptr compute_or_assign_intersection(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2);

continuous_set_ptr compute_or_assign_union(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2);

continuous_set_ptr compute_or_assign_difference(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2);

/** The cheap difference is empty if p1 is contained in p2 and p1 otherwise. */
continuous_set_ptr compute_or_assign_cheap_difference(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2);

continuous_set_ptr compute_cheap_difference(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2);

// --------------------------------------------
/** \name Operations specific to using continuous_sets as relations
 *  \{ */
// --------------------------------------------

/** Compute the concatenation of the two relations *this and *ps. Formally, the concatenation of two relations
 * \f$R\f$ and \f$S\f$ is \f$R \circ S = \{ (x,z) \mid \exists y : (x,y) \in R \wedge (y,z) \in S \}.\f$ */
continuous_set_ptr compute_or_assign_relation_concatenation(
		continuous_set_ptr p1, const continuous_set_const_ptr& p2);

/** \} */

// --------------------------------------------
/** \name Applying transformations
 *  \{ */
// --------------------------------------------

/** Returns the result of applying the continuous transformation t
 * (a map from a continuous set to a continuous set). */
continuous_set_ptr compute_transformation(continuous_set_const_ptr p,
		const continuous_set_transform_const_ptr& t);

/** Returns the result of applying the continuous transformation t
 * (a map from a continuous set to a continuous set).
 * If possible, the result is assigned to *this, otherwise a new set is created.
 * In both cases a ptr to the result is returned. */
continuous_set_ptr compute_or_assign_transformation(continuous_set_ptr p,
		const continuous_set_transform_const_ptr& t);

/** \} */

}

//#include "continuous_set_operator_implementations/compute_transformation.h"

#endif /*CONTINUOUS_SET_OPERATORS_H_*/
