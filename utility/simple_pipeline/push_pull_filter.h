/*
 * push_pull_pull_filter.h
 *
 *  Created on: Sep 28, 2009
 *      Author: frehse
 */

#ifndef PUSH_PULL_FILTER_H_
#define PUSH_PULL_FILTER_H_

#include "boost/shared_ptr.hpp"
#include <vector>

namespace simple_pipeline {

/** \brief A filter class for pipelines with push_pull semantics.
 *
 * In push_pull semantics, objects can be pushed through the
 * filter (the filter is called with an input object and propagates
 * it to the successors)
 * or pulled through the filter (the filter is called to produce
 * an output objects and requests any necessary input from its
 * predecessors).
 *
 * This class is flexible in the sense that one input object may
 * produce many output objects. However, it is assumed that
 * only one input object is necessary for an output to be produced.
 *
 * The objects that are passed through the filter are adopted by it,
 * i.e., it must delete any objects that are not propagated.
 * This may lead to problems if a filter has several successors:
 * one of the successors may delete the object while the others
 * are trying to access it.
 * The solution used in this class is to use smart pointers.
 * There is no need to actively delete or clone objects.
 *
 * @note A data source is a push_pull filter without any
 * predecessors. A data sink is a push_pull filter without any
 * successors.
 */

class push_pull_filter;
typedef boost::shared_ptr<push_pull_filter> push_pull_filter_ptr;

template<typename object_type> class push_pull_filter {
public:
	typedef boost::shared_ptr<object_type> object_type_ptr;

	push_pull_filter() {
	}
	;

	virtual ~push_pull_filter() {
	}
	;

	/** Accept *o as the input to the filter and prepare
	 * to return the filtered objects with get_output().
	 */
	virtual bool submit_input(object_type_ptr o) = 0;

	/** Return a filtered object. The input object is the last
	 * object submitted by submit_input().
	 */
	virtual object_type_ptr get_output() = 0;

	/** Return true if all filtered objects have been
	 * returned, i.e., if a new input object is needs to
	 * be submitted before any output objects can be
	 * produced.
	 */
	virtual bool input_exhausted() = 0;

	/** Return a predecessor from which to request the
	 * input object.
	 */
	virtual const push_pull_filter_ptr& get_predecessor() = 0;

	/** Return the next successor to propagate output.
	 * Every output object is pushed to successors given
	 * by this function until it returns 0.
	 */
	virtual const push_pull_filter_ptr& get_successor() = 0;

	/** Filter the object *o and propagate the resulting
	 * objects to successors until get_successor() returns 0.
	 * Returns true if filtering and all pushes were successful.
	 */
	bool push(object_type_ptr) {
		bool success = submit_input(o);
		while (success && !input_exhausted()) {
			const push_pull_filter_ptr& p = get_successor();
			while (success && p) {
				success = p->push(get_output());
				p = get_successor();
			}
		}
		return success;
	}
	;

	/** Returns a filtered object in o.
	 * The input to the filter is an object from the filter given by
	 * get_predecessor().
	 */
	bool pull(object_type_ptr) {
		success = true;
		if (input_exhausted()) {
			const push_pull_filter_ptr& p = get_predecessor();
			if (p) {
				success = p->pull(o) && submit_input(o);
			}
		}
		o = get_output();
		return success;
	}
	;

	/** Add a new predecessor. */
	void add_predecessor(push_pull_filter_ptr p) {
		register_predecessor(p);
		p->register_successor(get_ptr());
	}
	;
	/** Add a new successor. */
	void add_successor(push_pull_filter_ptr p) {
		register_successor(p);
		p->register_predecessor(get_ptr());
	}
	;
protected:
	/** Add a new predecessor. */
	void register_predecessor(push_pull_filter_ptr p) {
		my_predecessors.push_back(p);
	}
	;
	/** Add a new successor. */
	void register_successor(push_pull_filter_ptr p) {
		my_successors.push_pull_back(p);
	}
	;
private:
	merge_type my_merge_type;
	branch_type my_branch_type;
	std::vector<push_pull_filter_ptr> my_predecessors;
	std::vector<push_pull_filter_ptr> my_successors;
};

}

#endif /* PUSH_PULL_FILTER_H_ */
