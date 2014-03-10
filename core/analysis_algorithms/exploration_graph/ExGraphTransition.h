/*
 * ExGraphTransition.h
 *
 *  Created on: Jan 8, 2012
 *      Author: goyal
 */

#ifndef EXGRAPHTRANSITION_H_
#define EXGRAPHTRANSITION_H_

#include "ExGraphTransition_Info.h"

/** forward declarations */
namespace exploration_graph {

class ExGraphNode;
typedef boost::shared_ptr<ExGraphNode> ExNode_ptr;

class ExGraphTransition_Info;
typedef ExGraphTransition_Info::trans_type trans_type;
typedef ExGraphTransition_Info::trans_data_type trans_data_type;
typedef ExGraphTransition_Info::transId_type transId_type;
typedef ExGraphTransition_Info::interval_type interval_type;
typedef ExGraphTransition_Info::trans_bad_type trans_bad_type;
typedef ExGraphTransition_Info::ExGraphTransition_Info_ptr
		ExGraphTransition_Info_ptr;
}

/**
 * ExGraphTransition class represents an ExGraph transition.

 * An ExGraphTransition is defined by a successor node and additional information
 * such as transition type and transition data. */

namespace exploration_graph {

class ExGraphTransition {

public:

	/** destructor */
	~ExGraphTransition();

	/** constructors */
	ExGraphTransition();

	ExGraphTransition(const ExGraphNode& succNode, const trans_type& tType,
			const trans_data_type& tData);

	/** Return reference to the successor node. */
	const ExGraphNode& getSuccessorNode() const {
		return *successorNode;
	}

	/** Return transition information. */
	const ExGraphTransition_Info_ptr getInfo() const {
		return myInfo;
	}

private:

	ExNode_ptr successorNode; /*!< pointer to successor node */
	ExGraphTransition_Info_ptr myInfo; /*!< transition information */
};
}

#include "ExGraphNode.h"
#include "ExGraphTransition.hpp"
#endif /* EXGRAPHTRANSITION_H_ */
