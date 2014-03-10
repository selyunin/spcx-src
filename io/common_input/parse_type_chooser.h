/*
 * parse_type_chooser.h
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#ifndef PARSE_TYPE_CHOOSER_H_
#define PARSE_TYPE_CHOOSER_H_

#include "global/global_types.h"

/** This class stores globally which types the parser should instantiate. */
class parse_type_chooser {
public:
	typedef global_types::coefficient_type type;

	static type get_bool();
	static type get_number();
	static void set_bool(type t);
	static void set_number(type t);
private:
	static type my_bool;
	static type my_number;
};
#endif /* PARSE_TYPE_CHOOSER_H_ */
