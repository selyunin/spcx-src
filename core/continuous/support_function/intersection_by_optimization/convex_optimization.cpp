#include "convex_optimization.h"

namespace continuous {
namespace support_function{


template<>
Rational theta_functor<Rational>::compute(const Rational& theta) const{
	bool is_empty, is_bounded;
	math::vdom_vector<Rational> sv;

	unsigned int n = my_v.size();
	math::matrix<Rational> pi_transpose(n,2); // pi is the n-dim space to 2-dim space linear transformation matrix.

	// Initialize pi_transpose
	for(unsigned int i=0;i<n;i++){
		pi_transpose(i,0) = my_c[i];
		pi_transpose(i,1) = my_v[i];
	}

	math::vector<Rational> direction(2);
	direction[0] = Rational(std::cos(theta.get_double()));
	direction[1] = Rational(std::sin(theta.get_double()));

	// to be continued from here
	math::vector<Rational> trans_direction = pi_transpose.multiply_vector(direction);
	math::vdom_vector<Rational> my_vdom_v(trans_direction, my_v.get_index_to_variable_id_map());

	Rational sopt;
	my_S->compute_support(my_vdom_v, sopt, sv, is_empty, is_bounded);

	return (sopt - my_d * Rational(std::cos(theta.get_double())) )/Rational(std::sin(theta.get_double()));
};
}
}
