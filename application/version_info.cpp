/*
 * version_info.cpp
 *
 *  Created on: Oct 19, 2010
 *      Author: frehse
 */

#include "version_info.h"

namespace sspaceex {

std::string get_version() {
	return "0.9.8b";
}

std::string get_compilation_time() {
	return std::string(__DATE__) + ", " + std::string(__TIME__);
}

}
