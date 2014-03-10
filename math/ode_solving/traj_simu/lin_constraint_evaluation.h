#ifndef LIN_CONSTRAINT_EVALUATION_H
#define LIN_CONSTRAINT_EVALUATION_H

#include "math/vdom/lin_constraint.h"
#include "math/vdom/vdom_vector.h"

namespace math{

template<typename scalar_type>
scalar_type lin_constraint_evaluation(const lin_constraint<scalar_type> & c, const vdom_vector<scalar_type> & s);

template<typename scalar_type>
scalar_type lin_constraint_evaluation(const vdom_vector<scalar_type> & normal_v, scalar_type can_coeff, const vdom_vector<scalar_type> & s);

}

#include "math/ode_solving/traj_simu/lin_constraint_evaluation.hpp"

#endif
