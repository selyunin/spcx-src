#include "lin_constraint_evaluation.h"
#include "math/numeric/approx_comparator.h"

namespace math{

template<typename scalar_type>
scalar_type lin_constraint_evaluation(const lin_constraint<scalar_type> & c, const vdom_vector<scalar_type> & s){
	return lin_constraint_evaluation<scalar_type>(c.get_normal(),c.get_canonic_inh_coeff(),s);
}

template<typename scalar_type>
scalar_type lin_constraint_evaluation(const vdom_vector<scalar_type> & normal_v, scalar_type can_coeff, const vdom_vector<scalar_type> & s){
	scalar_type temp =  scalar_product(normal_v,s)+can_coeff;
	/*
	if ( math::numeric::approx_comparator<scalar_type>::is_maybe_zero(temp) )
		return scalar_type(0);
	else return  scalar_product(normal_v,s)+can_coeff;
	*/
	return  scalar_product(normal_v,s)+can_coeff;
}

}
