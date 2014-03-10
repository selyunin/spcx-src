#ifndef CONTINUOUS_SET_PPL_NNC_UTILITY_H_
#define CONTINUOUS_SET_PPL_NNC_UTILITY_H_

#include <string>
#include "boost/shared_ptr.hpp"

/** Forward declaration */
namespace ppl_polyhedron {
class continuous_set_PPL_NNC;
typedef boost::shared_ptr<continuous_set_PPL_NNC> continuous_set_PPL_NNC_ptr;
typedef boost::shared_ptr<const continuous_set_PPL_NNC> continuous_set_PPL_NNC_const_ptr;
}

namespace ppl_polyhedron {

/** Returns a pointer to a continuous_set_PPL_NNC constructed from the expr
 * by calling parse_predicate. */
continuous_set_PPL_NNC_ptr parse_PPL_NNC(std::string expr);

}

#endif /*CONTINUOUS_SET_PPL_NNC_UTILITY_H_*/
