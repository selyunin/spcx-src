#include "parse_exception.h"

namespace parser {

parse_exception::parse_exception(const std::string& msg) :
	basic_exception(msg) {
}

parse_exception::parse_exception(const std::string& msg,
		const parse_exception& cause) :
	basic_exception(msg, cause) {
}

parse_exception::parse_exception(const std::string& msg,
		const std::exception& cause) :
	basic_exception(msg, cause) {
}

}
