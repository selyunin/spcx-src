#ifndef CONTINUOUS_SET_DECLARATIONS_H_
#define CONTINUOUS_SET_DECLARATIONS_H_

#include "utility/typelist.h"
#include "global/global_types.h"

/** Forward declaration of all classes that derive from continuous_set and
 * that wish to accept a dispatcher. */

class Rational;
namespace continuous {
class support_function_provider;
template <typename scalar_type> class polyhedron;
template <typename scalar_type> class constr_polyhedron;
namespace support_function {
template <typename scalar_type> class sfm_cont_set;
}
template <typename scalar_type> class spacetime_flowpipe;
class predicate_continuous_set;

template<typename scalar_type> class continuous_set_simulation;
}
namespace ppl_polyhedron {
class continuous_set_PPL_NNC;
}

namespace continuous {
typedef typelist<support_function_provider,
typelist<ppl_polyhedron::continuous_set_PPL_NNC,
//typelist<constr_polyhedron,
typelist<constr_polyhedron<global_types::rational_type>,
typelist<constr_polyhedron<global_types::float_type>,
//typelist<sfm_cont_set<global_types::rational_type>,
typelist<support_function::sfm_cont_set<global_types::float_type>,
typelist<spacetime_flowpipe<global_types::float_type>,
typelist<predicate_continuous_set,
typelist<continuous_set_simulation<global_types::float_type>,
null_typelist>
> > > > > > > continuous_set_typelist;
/*
 typedef typelist<ppl_polyhedron::continuous_set_PPL_NNC,
 typelist<polyhedron<double>,
 typelist<polyhedron<Rational>,
 typelist<constr_polyhedron<double>,
 typelist<constr_polyhedron<Rational>,
 null_typelist>
 > > > > continuous_set_typelist;
 */
}

#endif /*CONTINUOUS_SET_DECLARATIONS_H_*/
