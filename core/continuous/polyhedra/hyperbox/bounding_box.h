/*
 * bounding_box.h
 *
 *  Created on: Nov 25, 2009
 *      Author: frehse
 */

// In order to break circular dependency with polyhedron_output.hpp,
// we apply a two-phase lock on the header
#if !defined(BOUNDING_BOX_H_) || !defined(BOUNDING_BOX_H__BODY)
//#ifndef BOUNDING_BOX_H_
#define BOUNDING_BOX_H_

#include "core/continuous/polyhedra/hyperbox/hyperbox.h"

#ifndef BOUNDING_BOX_H__BODY
#define BOUNDING_BOX_H__BODY

namespace continuous {

/** A bounding box is a hyperbox that tightly overapproximated a
 * continuous set.
 *
 * The bounds are computed only when needed (component-wise).
 * However, an affine transform requires all bounds.
 */
template<typename scalar_type>
class bounding_box: public hyperbox<scalar_type> {
public:
	typedef typename hyperbox<scalar_type>::value_type value_type;
	typedef typename hyperbox<scalar_type>::size_type size_type;
	typedef typename hyperbox<scalar_type>::point_type point_type;

	/** Construct a bounding box for the continuous set s.
	 *
	 * The set s has to be a support_function_provider.
	 */
	explicit bounding_box(const support_function_provider::const_ptr& s);
	virtual ~bounding_box();

	virtual const point_type& get_l() const;
	virtual const point_type& get_u() const;

protected:
	/** Compute the supremum on variable with index i (pos=true) or
	 * the infimum (neg=true).
	 */
	virtual value_type compute_sup(size_type j, bool pos) const;

	/** All access to elements of l and u is via get and set methods,
	 * use on-the-fly computation. */
	virtual const value_type& get_l(size_type i) const;
	virtual const value_type& get_u(size_type i) const;

private:
	void update_l(size_type i) const;
	void update_u(size_type i) const;
	void update_l() const;
	void update_u() const;

	support_function_provider::const_ptr my_root;
	bool my_l_is_uptodate;
	bool my_u_is_uptodate;
};

template<typename scalar_type>
hyperbox<scalar_type> compute_bounding_box(const support_function_provider& s);

template<typename scalar_type>
hyperbox<scalar_type> compute_bounding_box(const support_function_provider& s, const math::affine_map<scalar_type>& t);

}

#include "bounding_box.hpp"

#endif /* BOUNDING_BOX_H__BODY */
#endif /* BOUNDING_BOX_H_ */
