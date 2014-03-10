
#ifndef _CONVEX_OPT_HPP
#define _CONVEX_OPT_HPP

//#define GOLD_SAMPLE
//#define OPT_INIT_SAMPLE__

#include "convex_optimization.h"

namespace continuous {
namespace support_function{


template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type,functor>::sample_type convex_opt<scalar_type,functor>::sample(const scalar_type& lopt) const{

	LOGGER(DEBUG7,"convex_opt:sample:","function sampled at lambda = "+to_string(lopt));
	LOGGERSWOC(DEBUG5,"convex_opt:sample","convex optimization function sampled");
	scalar_type sopt = my_functor.compute(lopt);
	typename convex_opt<scalar_type,functor>::sample_type new_s = {lopt, sopt};
	return new_s;
}

template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type,functor>::sample_type convex_opt<scalar_type,functor >::line_intersection(
		convex_opt<scalar_type,functor>::sample_type s1, convex_opt<scalar_type,functor>::sample_type s2,
		convex_opt<scalar_type,functor>::sample_type p1, convex_opt<scalar_type,functor>::sample_type p2) {

//	LOGGER(DEBUG2,"convex_optimization::line_intersection","4 points are: " + to_string(s1)+", "+to_string(s2)+"and "+to_string(p1)+", "+to_string(p2));

//	std::cout << "convex_opt: s1.lambda:" << s1.lambda << ", s1.f_lambda:" << s1.f_lambda << std::endl;
//	std::cout << "convex_opt: s2.lambda:" << s2.lambda << ", s2.f_lambda:" << s2.f_lambda << std::endl;
//	std::cout << "convex_opt: p1.lambda:" << p1.lambda << ", p1.f_lambda:" << p1.f_lambda << std::endl;
//	std::cout << "convex_opt: p2.lambda:" << p2.lambda << ", p2.f_lambda:" << p2.f_lambda << std::endl;

	math::numeric::approx_comparator<scalar_type> my_comp;
	sample_type intersection;

	// Check if one of the lines is perpendicular to the x axis.

	if(my_comp.is_maybe_equal(s2.lambda, s1.lambda) && !(my_comp.is_maybe_equal(p2.lambda, p1.lambda))){
		if(my_comp.is_maybe_equal(p2.f_lambda,p1.f_lambda)){
			intersection.lambda = s1.lambda;
			intersection.f_lambda = p2.f_lambda;
		}
		else{
			scalar_type m2 = (p2.f_lambda - p1.f_lambda)/(p2.lambda - p1.lambda);
			scalar_type c2 = p1.f_lambda - m2*p1.lambda;
			intersection.lambda = s1.lambda;
			intersection.f_lambda = m2*intersection.lambda + c2;
		}
	}
	else if(my_comp.is_maybe_equal(p2.lambda, p1.lambda) && !(my_comp.is_maybe_equal(s2.lambda, s1.lambda))){
		if(my_comp.is_maybe_equal(s1.f_lambda, s2.f_lambda)){
			intersection.lambda = p1.lambda;
			intersection.f_lambda = s2.f_lambda;
		}
		else{
			scalar_type m1 = (s2.f_lambda - s1.f_lambda)/(s2.lambda - s1.lambda);
			scalar_type c1 = s1.f_lambda - m1*s1.lambda;
			intersection.lambda = p1.lambda;
			intersection.f_lambda = m1*intersection.lambda + c1;
		}
	}
	else if(my_comp.is_maybe_equal(p2.lambda, p1.lambda) && my_comp.is_maybe_equal(s2.lambda, s1.lambda)){
		throw std::runtime_error("line intersection:Cannot find intersection of parallel lines\n");
	}
	// Also check if the lines are perpendicular to y axis.

	else if(my_comp.is_maybe_equal(s2.f_lambda, s1.f_lambda) && !(my_comp.is_maybe_equal(p2.f_lambda, p1.f_lambda))){
		if(my_comp.is_maybe_equal(p2.lambda,p1.lambda)){
			intersection.lambda = p1.lambda;
			intersection.f_lambda = s1.f_lambda;
		}
		else{
			scalar_type m2 = (p2.f_lambda - p1.f_lambda)/(p2.lambda - p1.lambda);
			scalar_type c2 = p1.f_lambda - m2*p1.lambda;
			intersection.f_lambda = s1.f_lambda;
			intersection.lambda = (intersection.f_lambda - c2)/m2;
		}

	}
	else if(my_comp.is_maybe_equal(p2.f_lambda, p1.f_lambda) && !(my_comp.is_maybe_equal(s2.f_lambda, s1.f_lambda))){
		if(my_comp.is_maybe_equal(s1.lambda, s2.lambda)){
			intersection.lambda = s1.lambda;
			intersection.f_lambda = p1.f_lambda;
		}
		else{
			scalar_type m1 = (s2.f_lambda - s1.f_lambda)/(s2.lambda - s1.lambda);
			scalar_type c1 = s1.f_lambda - m1*s1.lambda;
			intersection.f_lambda = p1.f_lambda;
			intersection.lambda = (intersection.f_lambda - c1)/m1;
		}

	}
	else if(my_comp.is_maybe_equal(p2.f_lambda, p1.f_lambda) && (my_comp.is_maybe_equal(s2.f_lambda, s1.f_lambda))){
		throw std::runtime_error("line intersection:Cannot find intersection of parallel lines\n");
	}
	else{

		scalar_type m1 = (s2.f_lambda - s1.f_lambda)/(s2.lambda - s1.lambda);
		scalar_type m2 = (p2.f_lambda - p1.f_lambda)/(p2.lambda - p1.lambda);
		/*DEBUG*/
/*
		std::cout << "m1:" << m1 << std::endl ;
		std::cout << "m2" << m2 << std::endl ;
*/

		scalar_type c1 = s1.f_lambda - m1*s1.lambda;
	//	std::cout << "c1 = " << c1 << std::endl;
		scalar_type c2 = p1.f_lambda - m2*p1.lambda;
	//	std::cout << "c2 = " << c2 << std::endl;


		//std::cout << "c2-c1:" << c2-c1 <<std::endl;

		if(my_comp.is_maybe_equal(c1,c2)){ // not sure if this is what we should do.
	//		std::cout << "Inside c1 maybe equal c2" << std::endl;
			intersection.lambda = scalar_type(0);
			intersection.f_lambda = c1;
			return intersection;
		}
		else
			if(my_comp.is_maybe_equal(m1,m2)){
	//			std::cout << "m1:" << m1 << ", m2:" << m2 << std::endl;
				throw std::runtime_error("line intersection: Cannot find intersection of parallel lines\n");
			}
			else
				intersection.lambda = (c2 - c1)/(m1-m2);
	//	std::cout << "intersect.x " << intersection.lambda << std::endl;
		intersection.f_lambda = m1*intersection.lambda + c1;
	//	std::cout << "intersect.y " << intersection.f_lambda << std::endl;
	}

	return intersection;
}

template<class scalar_type, template<typename > class functor >
	template<typename float_type>
typename convex_opt<scalar_type,functor>::sample_type convex_opt<scalar_type,functor >::line_intersection_stable(
		sample_type s1, sample_type s2,
		sample_type p1, sample_type p2,
		bool& convergent) {

	math::numeric::approx_comparator<scalar_type> my_comp;
	sample_type intersection;

	float_type x1, x2, x3, x4, y1, y2, y3, y4;
	x1 = convert_element<float_type, scalar_type>(s1.lambda);
	x2 = convert_element<float_type, scalar_type>(s2.lambda);
	x3 = convert_element<float_type, scalar_type>(p1.lambda);
	x4 = convert_element<float_type, scalar_type>(p2.lambda);

	y1 = convert_element<float_type, scalar_type>(s1.f_lambda);
	y2 = convert_element<float_type, scalar_type>(s2.f_lambda);
	y3 = convert_element<float_type, scalar_type>(p1.f_lambda);
	y4 = convert_element<float_type, scalar_type>(p2.f_lambda);

	// Check if one of the lines is perpendicular to the x axis.

	if(my_comp.is_maybe_equal(s2.lambda, s1.lambda) && !(my_comp.is_maybe_equal(p2.lambda, p1.lambda))){
		if(my_comp.is_maybe_equal(p2.f_lambda,p1.f_lambda)){
			intersection.lambda = s1.lambda;
			intersection.f_lambda = p2.f_lambda;
		}
		else{
			float_type m2 = (y4 - y3)/(x4 - x3);
			float_type c2 = y3 - m2*x3;
			intersection.lambda = s1.lambda;
			intersection.f_lambda = convert_element<scalar_type, float_type>(m2*x1 + c2);
		}
	}
	else if(my_comp.is_maybe_equal(p2.lambda, p1.lambda) && !(my_comp.is_maybe_equal(s2.lambda, s1.lambda))){
		if(my_comp.is_maybe_equal(s1.f_lambda, s2.f_lambda)){
			intersection.lambda = p1.lambda;
			intersection.f_lambda = s2.f_lambda;
		}
		else{
			float_type m1 = (y2 - y1)/(x2 - x1);
			float_type c1 = y1 - m1*x1;
			intersection.lambda = p1.lambda;
			intersection.f_lambda = convert_element<scalar_type, float_type>(m1*x3 + c1);
		}
	}
	else if(my_comp.is_maybe_equal(p2.lambda, p1.lambda) && my_comp.is_maybe_equal(s2.lambda, s1.lambda)){
		if(my_comp.is_maybe_equal(p1.lambda, s1.lambda)){ // overlapping lines.
			convergent = true;
			return p1;
		}
		else
			throw std::runtime_error("line intersection:Cannot find intersection of parallel lines\n");
	}
	// Also check if the lines are perpendicular to y axis.

	else if(my_comp.is_maybe_equal(s2.f_lambda, s1.f_lambda) && !(my_comp.is_maybe_equal(p2.f_lambda, p1.f_lambda))){
		if(my_comp.is_maybe_equal(p2.lambda,p1.lambda)){
			intersection.lambda = p1.lambda;
			intersection.f_lambda = s1.f_lambda;
		}
		else{
			float_type m2 = (y4 - y3)/(x4 - x3);
			float_type c2 = y3 - m2*x3;
			intersection.f_lambda = s1.f_lambda;
			intersection.lambda = convert_element<scalar_type, float_type>((y1 - c2)/m2);
		}

	}
	else if(my_comp.is_maybe_equal(p2.f_lambda, p1.f_lambda) && !(my_comp.is_maybe_equal(s2.f_lambda, s1.f_lambda))){
		if(my_comp.is_maybe_equal(s1.lambda, s2.lambda)){
			intersection.lambda = s1.lambda;
			intersection.f_lambda = p1.f_lambda;
		}
		else{
			float_type m1 = (y2 - y1)/(x2 - x1);
			float_type c1 = y1 - m1*x1;
			intersection.f_lambda = p1.f_lambda;
			intersection.lambda = convert_element<scalar_type, float_type>((y3 - c1)/m1);
		}

	}
	else if (my_comp.is_maybe_equal(p2.f_lambda, p1.f_lambda) && (my_comp.is_maybe_equal(s2.f_lambda, s1.f_lambda))){
		if(my_comp.is_maybe_equal(p1.f_lambda,s1.f_lambda)){
			convergent = true; // Overlapping lines
			return p1;
		}
		throw std::runtime_error("line intersection:Cannot find intersection of parallel lines\n");
	}
	else{
		float_type x12,x34,y12,y34, a, b, c;

		x12 = x1 - x2;
		x34 = x3 - x4;
		y12 = y1 - y2;
		y34 = y3 - y4;

		c = (x12*y34) - (y12*x34);
		a = x1*y2 - y1*x2;
		b = x3*y4 - y3*x4;


		float_type x = (a*x34 - b*x12)/c;
		float_type y = (a*y34 - b*y12)/c;
		intersection.lambda  = convert_element<scalar_type, float_type>(x);
		intersection.f_lambda = convert_element<scalar_type, float_type>(y);
		if(my_comp.is_maybe_equal(intersection.lambda, s2.lambda)){
			// heuristics: intersection point is one of the line segment end point
			//convergent = true;
			return intersection;
		}
		if(my_comp.is_maybe_equal(intersection.lambda, p1.lambda)){
			//convergent = true;
			return intersection;
		}
	}

	return intersection;
}

template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type,functor>::sample_type convex_opt<scalar_type,functor>::sample_line(
		sample_type s1, sample_type s2, scalar_type lambda) {

	scalar_type m1 = (s2.f_lambda - s1.f_lambda)/(s2.lambda - s1.lambda);
	scalar_type c1 = s1.f_lambda - m1*s1.lambda;

	sample_type line_sample;
	line_sample.lambda = lambda;
	line_sample.f_lambda = m1*lambda + c1;

	return line_sample;
}

template<class scalar_type, template<typename > class functor>
bool convex_opt<scalar_type,functor>::is_colinear(sample_type s1, sample_type s2, sample_type s3){

//	std::cout << "colinear: s1.lambda:" << s1.lambda << ", s1.f_lambda:" << s1.f_lambda << std::endl;
//	std::cout << "colinear: s2.lambda:" << s2.lambda << ", s2.f_lambda:" << s2.f_lambda << std::endl;
//	std::cout << "colinear: s3.lambda:" << s3.lambda << ", s3.f_lambda:" << s3.f_lambda << std::endl;

	math::numeric::approx_comparator<scalar_type> my_comp;
	if(my_comp.is_maybe_equal(s1.lambda, s2.lambda) && my_comp.is_maybe_equal(s2.lambda, s3.lambda)) // 3 points on a vertical line
		return true;
	if(my_comp.is_maybe_equal(s1.f_lambda,s2.f_lambda) && my_comp.is_maybe_equal(s2.f_lambda,s3.f_lambda)){ // 3 points on a horizontal line
		return true;
	}
	scalar_type m1 = (s2.f_lambda - s1.f_lambda)/(s2.lambda - s1.lambda);
	scalar_type c1 = s1.f_lambda - m1*s1.lambda;

	if(my_comp.is_maybe_equal((m1*s3.lambda+c1), s3.f_lambda))
		return true;
	else
		return false;
}

template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type,functor>::sample_type convex_opt<scalar_type,functor>::sample_interval(const interval& my_interval, bool side) const {

	// mid point sample
	//scalar_type lambda = (my_interval.upper().get_val() + my_interval.lower().get_val())/2;
	scalar_type lambda = sample_request(my_interval,side);

	LOGGER(DEBUG7,"convex_opt:sample_interval:","function sampling interval: "+to_string(my_interval));
	LOGGER(DEBUG7,"convex_opt:sample_interval:","function sampled at lambda = "+to_string(lambda));
	scalar_type sopt = my_functor.compute(lambda);
	typename convex_opt<scalar_type,functor>::sample_type new_s = {lambda, sopt};
	return new_s;
}

template<class scalar_type, template<typename > class functor>
scalar_type convex_opt<scalar_type,functor>::sample_request(const interval& my_interval, bool side) const {

#ifdef GOLD_SAMPLE
	scalar_type d = my_interval.upper().get_val() - my_interval.lower().get_val();
	const scalar_type R(0.6180339887);
	scalar_type lambda;

	if(side){
		//left
		lambda = my_interval.upper().get_val() - d*R;
	}
	else{ //right
		lambda = my_interval.lower().get_val() + d*R;
	}
#else
//	mid point sample
	scalar_type lambda = (my_interval.upper().get_val() + my_interval.lower().get_val())/scalar_type(2.0);
#endif
	return lambda;
}

/**
 * Routine to bracket the minima with 3 points. We also return a 4th point for our
 * novel minima finding approach. This routine uses the constant factor 1.618034 to
 * magnify the downhill decend.
 *
 * @param t1
 * @param t2
 * @param t3
 * @param t4
 */
template<class scalar_type, template<typename > class functor>
void convex_opt<scalar_type,functor>::minbrak_simple(
		convex_opt<scalar_type, functor>::sample_type& t1,
		convex_opt<scalar_type, functor>::sample_type& t2,
		convex_opt<scalar_type, functor>::sample_type& t3,
		convex_opt<scalar_type, functor>::sample_type& t4,
		const unsigned int bound = 0){

	const double GOLD=1.618034;
	math::numeric::approx_comparator<scalar_type> comp;
 	t1 = sample(scalar_type(0)); // guessed initial starting points as 0 and 1.

 	if(my_dynamics_map.is_empty()){
 		t2 = sample(scalar_type(1));
 	}
 	else{
#ifdef OPT_INIT_SAMPLE__
 		std::list<scalar_type> init_lambdas  = get_init_sample();
 		t2 = sample(init_lambdas.front()); // intelligently choosing initial directions of sampling.
#else
 		t2 = sample(scalar_type(1));
#endif
 	}

	add_sample(t1.lambda, t1.f_lambda); // sample added to plotter
	add_sample(t2.lambda, t2.f_lambda); // sample added to plotter

	if(bound!=0 && this->get_size() >= bound)
		return;

	if(comp.is_maybe_equal(t2.f_lambda, t1.f_lambda)){
		t4 = t2;
		t2 = sample((t4.lambda - t1.lambda)/scalar_type(3));
		add_sample(t2.lambda,t2.f_lambda);
		if(bound!=0 && this->get_size() >= bound)
			return;
		t3 = sample(scalar_type(2)*(t4.lambda - t1.lambda)/scalar_type(3));
		add_sample(t3.lambda,t3.f_lambda);
/*
		std::cout << "minbrak_simple: t1 = t4 case reached\n";
		std::cout << "t1.f_lambda = " << t1.f_lambda << std::endl;
		std::cout << "t2.f_lambda = " << t2.f_lambda << std::endl;
		std::cout << "t3.f_lambda = " << t3.f_lambda << std::endl;
		std::cout << "t4.f_lambda = " << t4.f_lambda << std::endl;
*/
		return;
	}

	if(comp.is_definitely_strictly_larger(t2.f_lambda, t1.f_lambda)){
		swap(t1,t2);
	} // switch t1, t2 so that we can go downhill in the direction from a to b.
	t3  = sample(t2.lambda + (t2.lambda - t1.lambda)* scalar_type(GOLD)); // First guess for t3
	add_sample(t3.lambda, t3.f_lambda); // sample added to plotter

	if(bound!=0 && this->get_size() >= bound)
		return;

	while(comp.is_definitely_strictly_smaller(t3.f_lambda, t2.f_lambda)){
		shift2(t1,t2,t3);
		t3 = sample(t2.lambda + (t2.lambda - t1.lambda)*scalar_type(GOLD));
		add_sample(t3.lambda, t3.f_lambda); // sample added to plotter
		if(bound!=0 && this->get_size() >= bound)
			return;
	}
	if(comp.is_definitely_strictly_larger(t1.lambda, t2.lambda)){ // Rename the pivots to increasing order
		swap(t1,t3);
	}
	t4 = t3;
	t3 = sample(t2.lambda + (t4.lambda - t2.lambda)/scalar_type(2));
	add_sample(t3.lambda,t3.f_lambda);
}
;
/*
 * Brackets the function minima with parabolic extrapolation.
 * Adapted from Numerical Recipes in C++.
 */
template<class scalar_type, template<typename > class functor>
void convex_opt<scalar_type,functor>::minbrak_para_ext(
		convex_opt<scalar_type, functor>::sample_type& t1,
		convex_opt<scalar_type, functor>::sample_type& t2,
		convex_opt<scalar_type, functor>::sample_type& t3,
		convex_opt<scalar_type, functor>::sample_type& t4,
		const unsigned int bound = 0){

	const double GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
	scalar_type ulimit, u, r, q;

	/* GOLD is the default ratio by which successive intervals are magnified; GLIMIT is the maximum magnification allowed for a parabolic fit step.
	 */
	math::numeric::approx_comparator<scalar_type> comp;
	typename convex_opt<scalar_type, functor>::sample_type t, temp;
	t1 = sample(scalar_type(0)); // guessed initial starting points as 0 and 1.
	if(my_dynamics_map.is_empty())
		t2 = sample(scalar_type(1));
	else{
#ifdef 	OPT_INIT_SAMPLE__
		std::list<scalar_type> init_lambdas  = get_init_sample();
		t2 = sample(init_lambdas.front()); // Intelligently choosing initial directions of sampling.
#else
		t2 = sample(scalar_type(1));
#endif
	}

	add_sample(t1.lambda, t1.f_lambda); // sample added to plotter
	add_sample(t2.lambda, t2.f_lambda); // sample added to plotter
	if(bound!=0 && this->get_size() >= bound)
			return;

	if(comp.is_definitely_strictly_larger(t2.f_lambda, t1.f_lambda)){
		swap(t1,t2);
	} // switch t1, t2 so that we can go downhill in the direction from a to b.
	t3  = sample(t2.lambda + (t2.lambda - t1.lambda)* scalar_type(GOLD)); // First guess for t3
	add_sample(t3.lambda, t3.f_lambda); // sample added to plotter

	if(bound!=0 && this->get_size() >= bound)
			return;
	while(comp.is_definitely_strictly_smaller(t3.f_lambda,t2.f_lambda)){
		r = (t2.lambda - t1.lambda)*(t2.f_lambda - t3.f_lambda);
		q = (t2.lambda - t3.lambda)*(t2.f_lambda - t1.f_lambda);

		u = t2.lambda - ( ((t2.lambda - t3.lambda)*q - (t2.lambda - t1.lambda)*r)/
				(scalar_type(2.0) * SIGN(MAX(abs(q-r),scalar_type(TINY)),q-r)) );
		ulimit = t2.lambda + scalar_type(GLIMIT)*(t3.lambda - t2.lambda);// We wont go farther than this
		if((t2.lambda - u)*(u - t3.lambda) > scalar_type(0.0)) { //parabolic u is between b and c
			t = sample(u);
			add_sample(t.lambda, t.f_lambda);
			if(bound!=0 && this->get_size() >= bound)
					return;
			if(comp.is_definitely_strictly_smaller(t.f_lambda, t3.f_lambda)){
				t1 = t2;
				t2 = t;
				// Rename the pivots to increasing order
				if(t1.lambda > t2.lambda){
					swap(t1,t3);
				}
				// add another pivot point for the novel approach of finding minima
				t4 = sample(t3.lambda + scalar_type(5));
				add_sample(t4.lambda, t4.f_lambda); // sample added to plotter
				return;
			} else if(comp.is_definitely_strictly_larger(t.f_lambda, t2.f_lambda)) {
				t3 = t;
				// Rename the pivots to increasing order
				if(t1.lambda > t2.lambda){
					swap(t1,t3);
				}
				// add another pivot point for the novel approach of finding minima
				t4 = sample(t3.lambda + scalar_type(5));
				add_sample(t4.lambda, t4.f_lambda); // sample added to plotter
				return;
			}
			t = sample(t3.lambda + (t3.lambda - t2.lambda)*scalar_type(GOLD)); // parabolic fit was of no use
			add_sample(t.lambda, t.f_lambda);

			if(bound!=0 && this->get_size() >= bound)
					return;
		}
		else if( (t3.lambda - u) * (u - ulimit) > scalar_type(0.0)){
			t = sample(u);
			add_sample(t.lambda, t.f_lambda);
			if(bound!=0 && this->get_size() >= bound)
					return;
			temp = sample(t3.lambda+(t3.lambda-t2.lambda)* scalar_type(GOLD));
			add_sample(temp.lambda, temp.f_lambda);
			if(bound!=0 && this->get_size() >= bound)
					return;
			if(comp.is_definitely_strictly_smaller(t.f_lambda, t3.f_lambda)){
				shift3(t2, t3, t, temp);
			}
		}
		else if((u - ulimit)*(ulimit - t3.lambda) >=scalar_type(0.0)){ // Limit parabolic u to maximum allowed value
			u = ulimit;
			t = sample(u);
			add_sample(t.lambda, t.f_lambda);
			if(bound!=0 && this->get_size() >= bound)
					return;
		}
		else { // Reject parabolic u, use default magnification
			t = sample(t3.lambda + (t3.lambda - t2.lambda)*scalar_type(GOLD));
			add_sample(t.lambda, t.f_lambda);

			if(bound!=0 && this->get_size() >= bound)
					return;
		}
		shift3(t1,t2,t3,t);
	}
	// Rename the pivots to increasing order
	if(t1.lambda > t2.lambda){
		swap(t1,t3);
	}
	// add another pivot point for the novel approach of finding minima
	t4 = sample(t3.lambda + scalar_type(5));
	add_sample(t4.lambda, t4.f_lambda); // sample added to plotter
}
template<class scalar_type, template<typename > class functor>
std::list<scalar_type> convex_opt<scalar_type,functor>::get_init_sample() const {

	assert(!my_dynamics_map.is_empty());
	//assert(my_U);

	std::list<scalar_type> i_pivots;
	typename support_function_provider::const_ptr S = my_functor.get_S();

	// construct an sf provider object for AX+b+U

	typename support_function_provider::const_ptr S_map =
			typename support_function_provider::const_ptr(new support_function::sf_unary<scalar_type>(S,my_dynamics_map));
	typename support_function_provider::const_ptr S_map_sum_U = S_map;
	if (my_U) {
		S_map_sum_U = typename support_function_provider::const_ptr(new support_function::sf_sum<scalar_type>(S_map,my_U));
	}
	//
	bool is_empty,is_bounded;
	scalar_type max_val;
	math::vdom_vector<scalar_type> sup_vec;
	math::vdom_vector<scalar_type> l = my_functor.get_direction() - my_functor.get_cons_normal(); // l-n
	S_map_sum_U->compute_support(l,max_val,sup_vec,is_empty,is_bounded);
	//
	scalar_type lambda = sup_vec.scalar_product(my_functor.get_direction()) / sup_vec.scalar_product(my_functor.get_cons_normal());

	i_pivots.push_back(lambda);
	return i_pivots;
}

} // end of namespace support_function
} // end of namespace continuous

#endif /* _convex_opt_hpp*/
