#ifndef INTERSECTION_TRANSFORM_H_
#define INTERSECTION_TRANSFORM_H_

#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"

namespace continuous {

/** Intersects this with the set p.  */
class intersection_transform : public continuous_set_transform {
public:
	typedef boost::shared_ptr<intersection_transform> ptr;
	typedef boost::shared_ptr<const intersection_transform> const_ptr;

	intersection_transform(const continuous_set_ptr& p);
	virtual ~intersection_transform();

	const continuous_set_ptr& get_set() const;

	virtual void get_used_and_modif_variables(variable_id_set& used_vars,
			variable_id_set& modif_vars) const;

	virtual void print(std::ostream& os) const;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const;
protected:
	continuous_set_ptr my_set;
};

/** Pre- and post-intersection differ in parallel composition:
 * one type can not be composed with the other.
 * Pre-intersection is supposed to be carried out before the other operations,
 * while post-intersection is carried out after. */
class pre_intersection_transform : public intersection_transform {
public:
	pre_intersection_transform(const continuous_set_ptr& p);
	virtual ~pre_intersection_transform();
	virtual continuous_set_const_ptr get_relation(continuous_set_const_ptr cset) const;
};

class post_intersection_transform : public intersection_transform {
public:
	post_intersection_transform(const continuous_set_ptr& p);
	virtual ~post_intersection_transform();
	virtual continuous_set_const_ptr get_relation(continuous_set_const_ptr cset) const;
};

}

#endif /*INTERSECTION_TRANSFORM_H_*/
