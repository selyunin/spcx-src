/*
 * pull_filter.h
 *
 *  Created on: Sep 26, 2009
 *      Author: frehse
 */

#ifndef PULL_FILTER_H_
#define PULL_FILTER_H_

#include "boost/shared_ptr.hpp"
#include <vector>

namespace simple_pipeline {

/** \brief A filter class for pipelines with pull semantics.
 *
 * In pull semantics, a filter actively requests objects from
 * its predecessors.
 *
 * The filter is buffered: The objects from predecessors are
 * stored in a list and then processed after each predecessor
 * returns from the pull.
 *
 * @note A data source is a pull filter without any predecessors.
 */

class pull_filter;
typedef boost::shared_ptr<pull_filter> pull_filter_ptr;

template<typename object_type> class pull_filter {
public:
	typedef std::list<object_type*> container_type;

	/** An AND_MERGE is only successful if all predecessors have processed
	 * the object. An OR_MERGE tries each of the predecessors until one
	 * has processed it, and then stops.
	 */
	enum merge_type {
		AND_MERGE, OR_MERGE
	};

	pull_filter(merge_type bt = AND_MERGE) :
		my_merge_type(bt) {
	}
	;

	virtual ~pull_filter() {
	}
	;

	/** Change the merge type to bt. */
	void set_merge_type(merge_type bt) {
		my_merge_type = bt;
	}
	;

	/** The process method does its magic on the
	 * object, and calls prepagate() for all objects
	 * that might need to be passed on to the predecessors
	 * of the filter.
	 */
	virtual container_type process(object_type* o) = 0;

	/** Call all predecessor filters with the object.*/
	virtual container_type pull() {
		container_type res;
		for (std::vector<pull_filter_ptr>::iterator it = my_predecessors.begin(); it
				= my_predecessors.end(); ++it) {
			container_type pull_res = (*it)->pull();
			bool success=!pull_res.empty();
			// Process the pulled objects
			for (container_type::const_iterator pit=pull_res.begin();pit!=pull_res.end();++pit) {
				container_type proc_pull_res=process(*pit);
				// add the processed objects to the output
				res.splice(res.end(),proc_pull_res);
			}
			if (my_merge_type == AND_MERGE && !success)
				return container_type();
			if (my_merge_type == OR_MERGE && success)
				return res;
		}
		return res;
	}
	;

	/** Add a new predecessor. */
	void add_predecessor(pull_filter_ptr p) {
		my_predecessors.push_back(p);
	}
	;
private:
	merge_type my_merge_type;
	std::vector<pull_filter_ptr> my_predecessors;
};

}



#endif /* PULL_FILTER_H_ */
