/*
 * push_filter.h
 *
 *  Created on: Sep 26, 2009
 *      Author: frehse
 */

#ifndef PUSH_FILTER_H_
#define PUSH_FILTER_H_

#include "boost/shared_ptr.hpp"
#include <vector>

namespace simple_pipeline {

/** \brief A filter class for pipelines with push semantics.
 *
 * In push semantics, a filter actively propagates objects to
 * the next one.
 *
 * @note A data sink is a push filter without any successors.
 */

class push_filter;
typedef boost::shared_ptr<push_filter> push_filter_ptr;

template<typename object_type> class push_filter {
public:
	/** An AND_BRANCH is only successful if all successors have processed
	 * the object. An OR_BRANCH tries each of the successors until one
	 * has processed it, and then stops.
	 */
	enum branch_type {
		AND_BRANCH, OR_BRANCH
	};

	push_filter(branch_type bt = AND_BRANCH) :
		my_branch_type(bt) {
	}
	;

	virtual ~push_filter() {
	}
	;

	/** Change the branch type to bt. */
	void set_branch_type(branch_type bt) {
		my_branch_type = bt;
	}
	;

	/** The process method does its magic on the
	 * object, and calls propagate() for all objects
	 * that might need to be passed on to the successors
	 * of the filter.
	 */
	virtual bool process(object_type* o) = 0;

	bool push(object_type* o) {
		return process(o);
	}
	;

	/** Call all successor filters with the object.*/
	virtual bool propagate(object_type* o) const {
		for (std::vector<push_filter_ptr>::iterator it = my_successors.begin(); it
				= my_successors.end(); ++it) {
			bool success = (*it)->push(o);
			if (my_branch_type == AND_BRANCH && !success)
				return false;
			if (my_branch_type == OR_BRANCH && success)
				return true;
		}
		if (my_branch_type == OR_BRANCH) // none of the successors worked
			return false;
		else {
			// all successors worked
			return true;
		}
	}
	;

	/** Add a new successor. */
	void add_successor(push_filter_ptr p) {
		my_successors.push_back(p);
	}
	;
private:
	branch_type my_branch_type;
	std::vector<push_filter_ptr> my_successors;
};

}

#endif /* PUSH_FILTER_H_ */
