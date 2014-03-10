/*
 * ExGraphTransition_Info.h
 *
 *  Created on: Jan 19, 2012
 *      Author: goyal
 */

#ifndef EXGRAPHTRANSITION_INFO_H_
#define EXGRAPHTRANSITION_INFO_H_

#include "math/numeric/interval.h"
#include <boost/variant.hpp>

namespace exploration_graph {

/** The class represents the information contained in
 * an ExGraph transition.
 *
 * The information consists of two fields: (1) Transition type
 * (2) Transition data. */
class ExGraphTransition_Info {

public:

	typedef unsigned int transId_type; /*!< Discrete transition id type */

	typedef scalar_with_infinity<double> interval_data_type; /*!< Time Elapse interval data type */

	typedef math::numeric::interval<interval_data_type> interval_type; /*!< Time Elapse interval type */

	/** Types of ExGraph transitions */
	typedef enum {
		DISCRETE, TIME_ELAPSE, CONTAINMENT
	} trans_type;

	/** Signify whether a transition corresponds to a bad path or not.
	 *
	 * By default, a transition is not bad. */
	typedef enum {
		IS_BAD, IS_NOT_BAD
	} trans_bad_type;

	/** ExGraph transition data type: It can either be an id or a time interval. */
	typedef boost::variant<transId_type, interval_type> trans_data_type;

	/** Transition information. */
	typedef boost::shared_ptr<ExGraphTransition_Info>
			ExGraphTransition_Info_ptr;

	/** Constructors */
	ExGraphTransition_Info() {
		myType = CONTAINMENT;
	}

	ExGraphTransition_Info(trans_type tType, trans_data_type tData, trans_bad_type is_bad = IS_NOT_BAD) {
		myType = tType;
		myData = tData;
		myBad = is_bad;
	}

	/** Destructor */
	~ExGraphTransition_Info(void) {

	}

	/** Return true if two transition informations are same, else false. */
	bool IsEqual(ExGraphTransition_Info_ptr transInfo) {

		trans_type type1, type2;
		type1 = this->getType();
		type2 = transInfo->getType();

		trans_data_type data1, data2;
		data1 = this->getData();
		data2 = transInfo->getData();
		if((type1 == DISCRETE) && (type2 == DISCRETE)) {
			transId_type *tId1 = boost::get<transId_type>(&data1);
			transId_type *tId2 = boost::get<transId_type>(&data2);
			if(*tId1 == *tId2)
				return true;
			else
				return false;
		}
		else if((type1 == TIME_ELAPSE) && (type2 == TIME_ELAPSE)) {
			interval_type *interval1 = boost::get<interval_type>(&data1);
			interval_type *interval2 = boost::get<interval_type>(&data2);
			if(((*interval1).lower() == (*interval2).lower()) && ((*interval1).upper() == (*interval2).upper()))
				return true;
			else
				return false;
		}
		else if((type1 == CONTAINMENT) && (type2 == CONTAINMENT))
			return true;
		else
			return false;
	}

	/** Return ExGraph transition type. */
	trans_type getType() {
		return myType;
	}

	/** Return an ExGraph transition data. */
	trans_data_type getData() {
		return myData;
	}

	/** Return the bad type of a transition. */
	trans_bad_type getBadType() {
		return myBad;
	}

	/** set the bad type. */
	void setBadType(trans_bad_type bad) {
		myBad = bad;
	}

private:
	trans_type myType;  /*!< ExGraph transition type */
	trans_data_type myData; /*!< ExGraph transition data */
	trans_bad_type myBad; /*!< ExGraph transition bad type */
};

}

#endif /* EXGRAPHTRANSITION_INFO_H_ */
