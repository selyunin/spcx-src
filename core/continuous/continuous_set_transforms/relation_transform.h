#ifndef RELATION_TRANSFORM_H_
#define RELATION_TRANSFORM_H_

#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"

namespace continuous {

/** Represents the image compuation of a relation R,
 * i.e., \f$ \exists x' : (x,x') \in R\f$. The relation is represented as
 * a \p continous_set, in which the unprimed variables represent x and the
 * (singly-) primed variables represent x'.  */
class relation_transform : public continuous_set_transform {
public:
	typedef boost::shared_ptr<relation_transform> ptr;
	typedef boost::shared_ptr<const relation_transform> const_ptr;

	relation_transform(const relation_ptr& p);
	virtual ~relation_transform();

	virtual relation_const_ptr get_relation(continuous_set_const_ptr cset) const;
	virtual relation_ptr get_relation();
	virtual relation_const_ptr get_relation_const() const;

	virtual void get_used_and_modif_variables(variable_id_set& used_vars,
			variable_id_set& modif_vars) const;

	virtual void print(std::ostream& os) const;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const;

	relation_ptr relation;
};

}

#endif /*RELATION_TRANSFORM_H_*/
