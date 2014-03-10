/*
 * spacetime_plif.cpp
 *
 *  Created on: Mar 6, 2013
 *      Author: frehse
 */

#include "spacetime_plif.h"

#include "math/numeric/comp.h"

namespace spacetime {

spacetime_plif::error_type spacetime_plif::my_numeric_error =
		spacetime_plif::error_type(1e-08, 1e-12);

}

