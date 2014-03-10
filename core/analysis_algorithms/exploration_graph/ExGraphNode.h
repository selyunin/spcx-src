/*
 * ExGraphNode.h
 *
 *  Created on: Jan 8, 2012
 *      Author: goyal
 */

#ifndef EXGRAPHNODE_H_
#define EXGRAPHNODE_H_

#include "boost/shared_ptr.hpp"
#include <stdexcept>
#include <set>

namespace exploration_graph {

class ExGraphTransition;

/**
 * ExGraphNode class represents a node of an ExGraph.
 *
 * Each ExGraphNode defined by a vector of ExGraphTransition
 * which denotes the outgoing transitions from this ExGraphNode.
 */
class ExGraphNode: public boost::enable_shared_from_this<ExGraphNode> {

public:

	typedef boost::shared_ptr<ExGraphNode> ExNode_ptr;

	/** Node status. */
	typedef enum {
		IS_VISITED, IS_NOT_VISITED
	} status;

	/** Return a shared_ptr to *this. */
	ExNode_ptr get_ptr() {
		return boost::enable_shared_from_this<ExGraphNode>::shared_from_this();
	};

	/** ExGraphNode destructor */
	~ExGraphNode();

	/** ExGraphNode constructor */
	ExGraphNode();

	/** Return the vector of outgoing transitions of the ExGraphNode */
	std::vector<ExGraphTransition> getTransitions() const;

	/** Add a new ExGraph transition. */
	void addTransition(const ExGraphNode& succNode, const trans_type& tType,
			const trans_data_type& tData);

	/** Return the status of an ExGraph node. */
	status getStatus() const;

	/** Set the status of an ExGraph node. */
	void setStatus(const status st);

private:

	std::vector<ExGraphTransition> myTransitions; /*!< list of outgoing transitions of an ExGraph-node */

	status myStatus; /*!< ExGraph-node status */

	friend class ExGraph;
};
}

#include "ExGraphTransition.h"
#include "ExGraphNode.hpp"
#endif /* EXGRAPHNODE_H_ */
