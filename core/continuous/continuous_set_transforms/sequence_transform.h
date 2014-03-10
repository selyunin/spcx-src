#ifndef SEQUENCE_TRANSFORM_H_
#define SEQUENCE_TRANSFORM_H_

#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"

namespace continuous {

/** A set of transformations to be applied in sequence. */
class sequence_transform : public continuous_set_transform {
public:
	typedef boost::shared_ptr<sequence_transform> ptr;
	typedef boost::shared_ptr<const sequence_transform> const_ptr;

	typedef std::list<continuous_set_transform::const_ptr> transform_list_type;
	typedef std::list<continuous_set_transform::const_ptr>::const_iterator const_iterator;

	/** Initialize to an empty sequence. */
	sequence_transform();

	/** Initialize with the singleton sequence (t1), flattening it if t1 is itself a sequence. */
	sequence_transform(const continuous_set_transform::const_ptr& t1);

	/** Initialize with the sequence (t1,t2), joining them if one of them is a sequence. */
	sequence_transform(const continuous_set_transform::const_ptr& t1,
			const continuous_set_transform::const_ptr& t2);

	/** Initialize with the sequence (t1,t2,t3), joining them if one of them is a sequence. */
	sequence_transform(const continuous_set_transform::const_ptr& t1,
			const continuous_set_transform::const_ptr& t2, const continuous_set_transform::const_ptr& t3);

	/** Initialize with the sequence (*it1,...,*it2). */
	sequence_transform(const const_iterator& it1, const const_iterator& it2);

	virtual ~sequence_transform();

	const_iterator begin() const;

	const_iterator end() const;

	/** Returns the first element of the sequence. List must not be empty. */
	continuous_set_transform::const_ptr front() const;

	/** Returns the last element of the sequence. List must not be empty. */
	continuous_set_transform::const_ptr back() const;

	unsigned int size() const;

	void push_front(const continuous_set_transform::const_ptr& t);

	void push_back(const continuous_set_transform::const_ptr& t);

	virtual continuous_set_const_ptr get_relation(continuous_set_const_ptr cset) const;

	virtual void get_used_and_modif_variables(variable_id_set& used_vars,
			variable_id_set& modif_vars) const;

	virtual void print(std::ostream& os) const;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const;

private:
	transform_list_type transform_list;

	/** \todo compute set_difference  more effiently without tmp_set. */
	virtual void get_used_and_modif_variables(variable_id_set& used_vars,
			variable_id_set& modif_vars, const const_iterator& bit, const const_iterator& eit) const;

};

}

#endif /*SEQUENCE_TRANSFORM_H_*/
