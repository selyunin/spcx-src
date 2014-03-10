
/*
 * golden_section_search.cpp
 *
 *  Created on: Feb 21, 2010
 *  Author: ray
 */

//#include "golden_section_search.h"

namespace continuous {
namespace support_function{

template <class scalar_type>
scalar_type get_min(scalar_type t1, scalar_type t2, scalar_type t3)
{
	math::numeric::approx_comparator<scalar_type> my_comp;
	if((my_comp.is_definitely_strictly_smaller(t1,t2) || my_comp.is_maybe_equal(t1,t2))&&
			(my_comp.is_definitely_strictly_smaller(t1,t3) || my_comp.is_maybe_equal(t1,t3)))
		return t1;
	else if((my_comp.is_definitely_strictly_smaller(t2,t1) || my_comp.is_maybe_equal(t2,t1)) &&
			(my_comp.is_definitely_strictly_smaller(t2,t3) || my_comp.is_maybe_equal(t2,t3)))
		return t2;
	else
		return t3;
}

template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type, functor>::min_interval golden_section_search<scalar_type,functor>::section_search(
					const golden_section_search<scalar_type, functor>::interval& search_interval) {

	const scalar_type R(0.6180339887); // Golden Ratio
	const scalar_type C = 1-R;

	// Keeping count of number of samples.
	unsigned int no_samples = 0;
	math::numeric::approx_comparator<scalar_type> my_comp;
//	std::cout << "GSS: special_constant:" << special_constant << std::endl;

	typename convex_opt<scalar_type, functor>::sample_type lambda1,lambda2, lambda3, lambda4; // three pivots required for golden section search
	if( search_interval.upper().is_infinity() || search_interval.lower().is_infinity() ){
		throw std::runtime_error("golden_section_search: Cannot run golden section search algorithm on unbounded search interval");
	}

	if(search_interval.upper() == search_interval.lower()){
		interval my_x_interval(search_interval.upper(),search_interval.upper());
		lambda1 = sample(search_interval.upper().get_val());
		interval my_y_interval = interval(scalar_with_infinity<scalar_type>(lambda1.f_lambda));
		typename golden_section_search<scalar_type, functor>::min_interval min_structure = {my_x_interval,my_y_interval};
		return min_structure;
	}

	scalar_type tolerance = convex_opt<scalar_type, functor>::get_interval_tolerance();
	assert(tolerance > 0.0); // This to be replaced with bound = machine's double precision
							 // given that tolerance is equal or less than the machines double precision.
//	std::cout << "GS: tolerance:" << tolerance << std::endl;
	scalar_with_infinity<scalar_type> upper(scalar_with_infinity<scalar_type>::pos_infty());
	scalar_with_infinity<scalar_type> lower(scalar_with_infinity<scalar_type>::neg_infty());


	lambda1.lambda = search_interval.lower().get_val();
	lambda4.lambda = search_interval.upper().get_val();
	scalar_type d = lambda4.lambda - lambda1.lambda;
	scalar_type x = d * R;

	lambda3.lambda = lambda1.lambda + x;
	lambda2.lambda = lambda4.lambda - x;

//	std::cout << "GS:lambda1:" << lambda1.lambda << std::endl;
//	std::cout << "GS:lambda2:" << lambda2.lambda << std::endl;
//	std::cout << "GS:lambda3:" << lambda3.lambda << std::endl;

	lambda2 = sample(lambda2.lambda);
	lambda3 = sample(lambda3.lambda);
	add_sample(lambda2.lambda,lambda2.f_lambda);
	add_sample(lambda3.lambda, lambda3.f_lambda);
	no_samples+=2;

//	std::cout << "GSS: L1 = " << std::setprecision(6) << lambda1.lambda << std::endl;
//	std::cout << "GSS: L2 = " << std::setprecision(6) << lambda2.lambda << std::endl;
//	std::cout << "GSS: L3 = " << std::setprecision(6) << lambda3.lambda << std::endl;
//	std::cout << "GSS: L4 = " << std::setprecision(6) << lambda4.lambda << std::endl;

	unsigned int cnt = 0;
	while(my_comp.is_definitely_strictly_larger(d , tolerance)){ // We keep track of 4 pivot points at any time.
		if(my_comp.is_definitely_strictly_smaller(lambda2.f_lambda, lambda3.f_lambda)) {
				lambda4 = lambda3;
				lambda3 = lambda2;
				lambda2.lambda = R*lambda2.lambda + C*lambda1.lambda;
				no_samples++;
				lambda2 = sample(lambda2.lambda);
				add_sample(lambda2.lambda, lambda2.f_lambda );
		}
		else{
				lambda1 = lambda2;
				lambda2 = lambda3;
				lambda3.lambda = R*lambda3.lambda+C*lambda4.lambda;
				lambda3 = sample(lambda3.lambda);
				no_samples++;
				add_sample(lambda3.lambda, lambda3.f_lambda);
		}

//		std::cout << "GSS: L1 = " << std::setprecision(6) << lambda1.lambda << std::endl;
//		std::cout << "GSS: L2 = " << std::setprecision(6) << lambda2.lambda << std::endl;
//		std::cout << "GSS: L3 = " << std::setprecision(6) << lambda3.lambda << std::endl;
//		std::cout << "GSS: L4 = " << std::setprecision(6) << lambda4.lambda << std::endl;

		d = lambda4.lambda - lambda1.lambda;
//		std::cout << "d:" << d <<  std::endl;
//		cnt++;
		//cpmparison purpose
		//if(no_samples == 6)
		//	break;
		//-----

	}
	if(my_comp.is_definitely_strictly_smaller(lambda2.f_lambda, lambda3.f_lambda))
		upper = scalar_with_infinity<scalar_type>(lambda2.f_lambda);
	else
		upper = scalar_with_infinity<scalar_type>(lambda3.f_lambda);

	interval my_x_interval(scalar_with_infinity<scalar_type>(lambda1.lambda),scalar_with_infinity<scalar_type>(lambda4.lambda));
	interval my_y_interval;
	my_y_interval = interval(upper,upper);
	typename golden_section_search<scalar_type, functor>::min_interval min_structure = {my_x_interval,my_y_interval};
	//LOGGER(DEBUG1,"GOLDEN_SECTION:","no_samples = "+to_string(no_samples));
	//std::cout<< "GOLDEN_SECTION: no_samples = " << no_samples << std::endl;
	return min_structure;
};

template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type, functor>::min_interval golden_section_search<scalar_type,functor>::section_search() {
	throw std::runtime_error("golden_section_search: Cannot run golden section search algorithm on unbounded search interval");
}


}//end of support_function namespace
}//end of continuous namespace

