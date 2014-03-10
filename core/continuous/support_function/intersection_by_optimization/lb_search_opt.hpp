#ifndef _LB_SEARCH_OPT_HPP
#define _LB_SEARCH_OPT_HPP

namespace continuous{
	namespace support_function{
/**
 * Checks if the passed lambda is out of the bracketing interval, in which case,
 * the function return true saying the the passed lambda is redundant.
 * Otherwise, this returns false.
 *
 * @param lambda
 * @return
 */
template<class scalar_type, template<typename > class functor >
bool lb_search_opt<scalar_type, functor>::sample_redundant(scalar_type lambda){
	math::numeric::approx_comparator<scalar_type> my_comp;
	if(prob_state != four_point){
		if(my_comp.is_definitely_strictly_smaller(lambda, s1.lambda) || my_comp.is_definitely_strictly_larger(lambda, s4.lambda)){
//			std::cout << "debug position 1:" << std::endl;
			return true;
		}
	}
	else{
		if(my_comp.is_definitely_strictly_smaller(lambda, s1.lambda) || my_comp.is_definitely_strictly_larger(lambda, s5.lambda)){
//			std::cout << "debug position 2:" << std::endl;
//			std::cout << "sample.x=" << sample.lambda << "sample.y=" << sample.f_lambda <<  std::endl;
			return true;
		}
	}
	// check if the passed sample is one of the sample points.
	if(my_comp.is_maybe_equal(lambda, s1.lambda) ||
	   my_comp.is_maybe_equal(lambda, s2.lambda) ||
	   my_comp.is_maybe_equal(lambda, s3.lambda) ||
	   my_comp.is_maybe_equal(lambda, s4.lambda))
		return true;

	if(prob_state == four_point && my_comp.is_maybe_equal(lambda, s5.lambda)){
//		std::cout << "debug position 4:" << std::endl;
		return true;
 	}
	return false;
}
/**
 * Chooses a sampling point out of the many choices. The choice could be based on a
 * number of criteria. Currently, the first sampling point in the list is choosen
 * arbitrarily.
 *
 * @param lambdas
 * @return
 */
template<class scalar_type, template<typename > class functor>
scalar_type lb_search_opt<scalar_type,functor>::choose_lambda(std::list<scalar_type>& lambdas){
	if(!lambdas.empty()){
		return lambdas.front();
	}// the front choice is returned. Blind decision, can be improved
	else
		throw std::runtime_error("lb_search_opt: lambdas choice list empty\n");
}
;
/**
 * Maps parameter lambda to a direction, given by, dir = lambda*c-v
 *
 * @param v direction
 * @param c Normal to the guard
 * @param lambda parameter
 * @return direction given by (lambda*c-v)
 */
template<typename scalar_type, template<typename > class functor>
typename math::vector<scalar_type> lb_search_opt<scalar_type,functor>::map_to_direction(const math::vdom_vector<scalar_type>& v,const math::vdom_vector<scalar_type>& c, scalar_type lambda){
	math::vector<scalar_type> my_v, my_c;
	my_v = v.get_vector();
	my_c = c.get_vector();
	LOGGER(DEBUG7,"lb_search_opt:map to direction","lambda:"+to_string(lambda)+",l:"+to_string(my_v)+",n:"+to_string(my_c)+"maps to:"+to_string(my_v - lambda*my_c));
	return my_v - lambda*my_c;
};

/**
 * Given a function sample, this function updates the bracketing interval to
 * close down on the function minimum.
 *
 * @param sample
 * @return
 */
template<class scalar_type, template<typename > class functor>
typename math::numeric::interval<scalar_type> lb_search_opt<scalar_type,functor>::update_bounds(lb_search_opt<scalar_type,functor>::sample_type sample){
	math::numeric::approx_comparator<scalar_type> my_comp;
	typename lb_search_opt<scalar_type, functor>::sample_type intersect1,intersect2,intersect,s;
	scalar_type candidate;

	// Check if the passed sample lies between s1 and s2.
	if(my_comp.is_definitely_strictly_larger(sample.lambda, s1.lambda) &&
	   my_comp.is_definitely_strictly_smaller(sample.lambda, s2.lambda)){
//		std::cout << "new sample between s1,s2" <<std::endl;
		if(prob_state == one){
			s1 = sample;
			s.lambda = s4.lambda;
			s.f_lambda = scalar_type(0);

			// by convexity, sample.f_lambda has to be > s2.f_lambda.

			// Update lower bound
			bool conv_pts = false;
			intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s,s4,conv_pts);
			if(conv_pts){// heuristic to speed up computation and to bypass numerical error
				// minima reached heuristics
				prob_state = stop;
				min_bounds.set_lower(s3.f_lambda);
				min_bounds.set_upper(s3.f_lambda);
				return min_bounds;
			}
			if(my_comp.is_definitely_strictly_smaller(intersect1.f_lambda, intersect2.f_lambda)){
				intersect = intersect1;
			}
			else
				intersect = intersect2;

			if(min_bounds.lower().is_finite()){
				get_maximum(min_bounds.lower().get_val(), intersect.f_lambda);
			}
			else{
				min_bounds.set_lower(intersect.f_lambda);
			}
			return min_bounds;
		}
		else if(prob_state == two){
			if(my_comp.is_definitely_strictly_larger(sample.f_lambda, s2.f_lambda)){
				// we go to 5 point state
				s5 = s4;
				s4 = s3;
				s3 = s2;
				s2 = sample;
				prob_state = four_point;
				// Change the bounds
				bool conv_pts = false;
				intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5,conv_pts);

				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(my_comp.is_definitely_strictly_smaller(intersect1.f_lambda,intersect2.f_lambda) )
					candidate = intersect1.f_lambda;
				else if(my_comp.is_definitely_strictly_larger(intersect1.f_lambda, intersect2.f_lambda))
					candidate = intersect2.f_lambda;
				else if(my_comp.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}
				else // choose one of the triangles arbitrarily.
					candidate = intersect1.f_lambda;

				//update lower bound
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(candidate, min_bounds.lower().get_val()))
						min_bounds.set_lower(candidate);
				}
				else
					min_bounds.set_lower(candidate);

				//update upper bound
				if(min_bounds.upper().is_finite()){
					if(my_comp.is_definitely_strictly_smaller(s3.f_lambda, min_bounds.upper().get_val()))
						min_bounds.set_upper(s3.f_lambda);
				}
				else
					min_bounds.set_upper(s3.f_lambda);
				return min_bounds;
			}
			else if(my_comp.is_definitely_strictly_smaller(sample.f_lambda, s2.f_lambda)){
				s4 = s3;
				s3 = s2;
				s2 = sample;
				//Update the upper bound
				candidate = s2.f_lambda;
				if(min_bounds.upper().is_finite()){
					if(my_comp.is_definitely_strictly_smaller(candidate, min_bounds.upper().get_val()))
						min_bounds.set_upper(candidate);
				}
				else
					min_bounds.set_upper(candidate);

				//Update the lower bound

				typename lb_search_opt<scalar_type, functor>::sample_type min_1,min_2,s;
				s.lambda = s1.lambda;
				s.f_lambda = scalar_type(0);
				bool conv_pts = false;
				min_1 = this->template line_intersection_stable<precise_float>(s2,s3,s,s1,conv_pts);
				min_2 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s2.f_lambda);
					min_bounds.set_upper(s2.f_lambda);
					return min_bounds;
				}

				if(min_bounds.lower().is_finite()){
					min_bounds.set_lower(get_maximum(min_bounds.lower().get_val(),get_minimum(min_1.f_lambda,min_2.f_lambda)));
				}
				else{
					min_bounds.set_lower(get_minimum(min_1.f_lambda,min_2.f_lambda));
				}

				return min_bounds;

			}
			else if(my_comp.is_maybe_equal(sample.f_lambda, s2.f_lambda)){
//			else if(sample.f_lambda == s2.f_lambda){
				s4 = s3;
				s3 = s2;
				s2 = sample;
				prob_state = three;

				bool conv_pts = false;
				intersect = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);

				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(intersect.f_lambda == s2.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s2.f_lambda);
					min_bounds.set_upper(s2.f_lambda);
					return min_bounds;
				}
				// Change the lower bound
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(intersect.f_lambda, min_bounds.lower().get_val()))
						min_bounds.set_lower(intersect.f_lambda);
				}
				else
					min_bounds.set_lower(intersect.f_lambda);
				//Change the upper bound
				if(min_bounds.upper().is_finite()){
					if(my_comp.is_definitely_strictly_smaller(s2.f_lambda, min_bounds.upper().get_val()))
						min_bounds.set_upper(s2.f_lambda);
				}
				else
					min_bounds.set_upper(s2.f_lambda);

				return min_bounds;
			}
			else{
				throw std::runtime_error("lb_search:STATE 2: update_bounds: should not reach here by convexity\n");
			}
		}
		else if(prob_state == three){
			if(my_comp.is_maybe_equal(sample.f_lambda, s2.f_lambda)){
//			if(sample.f_lambda == s2.f_lambda){
				prob_state = stop; // minimum reached state
				min_bounds.set_lower(sample.f_lambda);
				min_bounds.set_upper(sample.f_lambda);
				return min_bounds;
			}
			else if(my_comp.is_definitely_strictly_larger(sample.f_lambda, s2.f_lambda)){
				s1 = sample;
				// Change the bounds
				typename lb_search_opt<scalar_type, functor>::sample_type intersect;
				bool conv_pts = false;
				intersect = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);

				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

/*
				if(intersect.f_lambda == s2.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s2.f_lambda);
					min_bounds.set_upper(s2.f_lambda);
					return min_bounds;
				}
*/
				//change the lower bound;
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(intersect.f_lambda, min_bounds.lower().get_val()))
						min_bounds.set_lower(intersect.f_lambda);
				}
				else
					min_bounds.set_lower(intersect.f_lambda);
				//upper bound will not change
				return min_bounds;
			}
			else{
				throw std::runtime_error("update_interval: STATE 3: by convexity, this state should not be reached!");
			}
		}
		else if(prob_state == four_point){
			if(my_comp.is_definitely_strictly_larger(sample.f_lambda, s2.f_lambda)){
				s1 = sample;
				// Change the bounds here
				bool conv_pts = false;
				intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5,conv_pts);

				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(my_comp.is_definitely_strictly_smaller(intersect1.f_lambda,intersect2.f_lambda) )
					candidate = intersect1.f_lambda;
				else if(my_comp.is_definitely_strictly_larger(intersect1.f_lambda, intersect2.f_lambda))
					candidate = intersect2.f_lambda;
				else if(my_comp.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
	//			else if(intersect1.f_lambda == s3.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				else // choose one of the triangles arbitrarily.
					candidate = intersect1.f_lambda;

				//update lower bound
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(candidate, min_bounds.lower().get_val()))
						min_bounds.set_lower(candidate);
				}
				else
					min_bounds.set_lower(candidate);
				//upper bound doesn't change
				return min_bounds;
			}
			else{
				print_pivots();
				this->plot_graph();
				this->plot_function(-1,0.5,0.001);
				std::cout << "new_sample.x:" << sample.lambda << ", new_sample.fx:" << sample.f_lambda << std::endl;
				throw std::runtime_error("update_bounds: State 4: should not reach here due to convexity");
			}
		}
	}
	// Check if the passed sample lies between s2 and s3.
	else if(my_comp.is_definitely_strictly_larger(sample.lambda, s2.lambda) &&
			my_comp.is_definitely_strictly_smaller(sample.lambda, s3.lambda)){
//		std::cout << "new sample between s2,s3" <<std::endl;
		//debug:
//		std::cout << "prob_state:" << get_problem_state();
//		print_pivots();
		if(prob_state == one || prob_state == four_point){
			if(my_comp.is_definitely_strictly_larger(sample.f_lambda, s3.f_lambda)){
				s1 = s2;
				s2 = sample;
				// Change the bounds
				bool conv_pts = false;
				if(prob_state == four_point){
					intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
					intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5,conv_pts);
					if(conv_pts){// heuristic to speed up computation and to bypass numerical error
						// minima reached heuristics
						prob_state = stop;
						min_bounds.set_lower(s3.f_lambda);
						min_bounds.set_upper(s3.f_lambda);
						return min_bounds;
					}

					if(my_comp.is_definitely_strictly_smaller(intersect1.f_lambda,intersect2.f_lambda) )
						candidate = intersect1.f_lambda;
					else if(my_comp.is_definitely_strictly_larger(intersect1.f_lambda, intersect2.f_lambda))
						candidate = intersect2.f_lambda;
					else if(my_comp.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
//					else if(intersect1.f_lambda == s3.f_lambda){
						prob_state = stop;
						min_bounds.set_lower(s3.f_lambda);
						min_bounds.set_upper(s3.f_lambda);
						return min_bounds;
					}
					else // choose one of the triangles arbitrarily.
						candidate = intersect1.f_lambda;
					// update lower bounds.
					if(min_bounds.lower().is_finite()){
						if(my_comp.is_definitely_strictly_larger(candidate, min_bounds.lower().get_val()))
							min_bounds.set_lower(candidate);
					}
					else
						min_bounds.set_lower(candidate);

					return min_bounds;
				}
				else{ // prob_state == one
					typename lb_search_opt<scalar_type, functor>::sample_type s;
					bool conv_pts = false;
					s.lambda = s4.lambda;
					s.f_lambda = scalar_type(0);
					intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
					intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s,s4,conv_pts);

					if(conv_pts){// heuristic to speed up computation and to bypass computation error
						// minima reached heuristics
						prob_state = stop;
						min_bounds.set_lower(s3.f_lambda);
						min_bounds.set_upper(s3.f_lambda);
						return min_bounds;
					}

					if(min_bounds.lower().is_finite()){
						min_bounds.set_lower(get_maximum(min_bounds.lower().get_val(),get_minimum(intersect1.f_lambda,intersect2.f_lambda)));
					}
					else{
						min_bounds.set_lower(get_minimum(intersect1.f_lambda,intersect2.f_lambda));
					}

				}
			}
			else if(my_comp.is_maybe_equal(sample.f_lambda, s3.f_lambda)){
//			else if(sample.f_lambda == s3.f_lambda){
				s1 = s2;
				s2 = sample;
				prob_state = three;
				// Change the bounds
				bool conv_pts = false;
				intersect = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(intersect.f_lambda == s2.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}
				//update lower bound
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(intersect.f_lambda, min_bounds.lower().get_val()))
						min_bounds.set_lower(intersect.f_lambda);
				}
				else
					min_bounds.set_lower(intersect.f_lambda);
				return min_bounds;
			}
			else if(my_comp.is_definitely_strictly_smaller(sample.f_lambda, s3.f_lambda)){
				s5 = s4;
				s4 = s3;
				s3 = sample;
				prob_state = four_point;
				// Change the bounds
				//change the lower bound
				bool conv_pts = false;
				intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5,conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(my_comp.is_definitely_strictly_smaller(intersect1.f_lambda,intersect2.f_lambda) )
					candidate = intersect1.f_lambda;
				else if(my_comp.is_definitely_strictly_larger(intersect1.f_lambda, intersect2.f_lambda))
					candidate = intersect2.f_lambda;
				else if(my_comp.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
//				else if(intersect1.f_lambda == s3.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}
				else // choose one of the triangles arbitrarily.
					candidate = intersect1.f_lambda;
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(candidate, min_bounds.lower().get_val()))
						min_bounds.set_lower(candidate);
				}
				else
					min_bounds.set_lower(candidate);
				//change the upper bound
				if(min_bounds.upper().is_finite()){
					if(my_comp.is_definitely_strictly_smaller(s3.f_lambda, min_bounds.upper().get_val()))
						min_bounds.set_upper(s3.f_lambda);
				}
				else
					min_bounds.set_upper(s3.f_lambda);
				return min_bounds;
			}
			else{
				print_pivots();
				this->plot_graph();
				this->plot_function(-1,0.5,0.001);
				std::cout << "new_sample.x:" << sample.lambda << ", new_sample.fx:" << sample.f_lambda << std::endl;
				throw std::runtime_error("update_interval: this state should not be reached!");
			}
		}
		else if(prob_state == two){
			if(my_comp.is_definitely_strictly_larger(sample.f_lambda, s2.f_lambda)){
				s4 = s3;
				s3 = sample;
				// update the lower bound

				typename lb_search_opt<scalar_type, functor>::sample_type s;
				s.lambda = s1.lambda;
				s.f_lambda = scalar_type(0);
				bool conv_pts = false;
				intersect1 = this->template line_intersection_stable<precise_float>(s2,s3,s,s1,conv_pts);
				intersect2 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass computation error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s2.f_lambda);
					min_bounds.set_upper(s2.f_lambda);
					return min_bounds;
				}

				if(min_bounds.lower().is_finite()){
					min_bounds.set_lower(get_maximum(min_bounds.lower().get_val(),get_minimum(intersect1.f_lambda,intersect2.f_lambda)));
				}
				else{
					min_bounds.set_lower(get_minimum(intersect1.f_lambda,intersect2.f_lambda));
				}
				return min_bounds;
			}
			else if(my_comp.is_maybe_equal(sample.f_lambda, s2.f_lambda)){
//			else if(sample.f_lambda == s2.f_lambda){
				s4 = s3;
				s3 = sample;
				prob_state = three;
				// Change the bounds
				// only the lower bound may change
				bool conv_pts = false;
				intersect = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s2.f_lambda);
					min_bounds.set_upper(s2.f_lambda);
					return min_bounds;
				}

				if(intersect.f_lambda == s2.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}
				//update the lower bound
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(intersect.f_lambda, min_bounds.lower().get_val()))
						min_bounds.set_lower(intersect.f_lambda);
				}
				else
					min_bounds.set_lower(intersect.f_lambda);
				return min_bounds;
			}
			else if(my_comp.is_definitely_strictly_smaller(sample.f_lambda, s2.f_lambda)){
				s5 = s4;
				s4 = s3;
				s3 = sample;
				prob_state = four_point;
				// Change the bounds
				//change the lower bound
				bool conv_pts = false;
				intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4, conv_pts);
				intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5, conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(my_comp.is_definitely_strictly_smaller(intersect1.f_lambda,intersect2.f_lambda) )
					candidate = intersect1.f_lambda;
				else if(my_comp.is_definitely_strictly_larger(intersect1.f_lambda, intersect2.f_lambda))
					candidate = intersect2.f_lambda;
				else if(my_comp.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
//				else if(intersect1.f_lambda == s3.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}
				else // choose one of the triangles arbitrarily.
					candidate = intersect1.f_lambda;
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(candidate, min_bounds.lower().get_val()))
						min_bounds.set_lower(candidate);
				}
				else
					min_bounds.set_lower(candidate);

				//update the lower bound
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(candidate, min_bounds.lower().get_val()))
						min_bounds.set_lower(candidate);
				}
				else
					min_bounds.set_lower(candidate);
				//change the upper bound
				if(min_bounds.upper().is_finite()){
					if(my_comp.is_definitely_strictly_smaller(s3.f_lambda, min_bounds.upper().get_val()))
						min_bounds.set_upper(s3.f_lambda);
				}
				else
					min_bounds.set_upper(s3.f_lambda);

				return min_bounds;
			}
			else{
				print_pivots();
				this->plot_graph();
				this->plot_function(-1,0.5,0.001);
				std::cout << "new_sample.x:" << sample.lambda << ", new_sample.fx:" << sample.f_lambda << std::endl;
				throw std::runtime_error("update_interval: this state should not be reached!");
			}
		}
		else if(prob_state == three){
			if(my_comp.is_definitely_strictly_smaller(sample.f_lambda, s2.f_lambda)){
				s5 = s4;
				s4 = s3;
				s3 = sample;
				prob_state = four_point;
				// Change the bounds
				//change the lower bound
				bool conv_pts = false;
				intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5, conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(my_comp.is_definitely_strictly_smaller(intersect1.f_lambda,intersect2.f_lambda) )
					candidate = intersect1.f_lambda;
				else if(my_comp.is_definitely_strictly_larger(intersect1.f_lambda, intersect2.f_lambda))
					candidate = intersect2.f_lambda;
				else if(my_comp.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
//				else if(intersect1.f_lambda == s3.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}
				else // choose one of the triangles arbitrarily.
					candidate = intersect1.f_lambda;

				//update lower bound
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(candidate, min_bounds.lower().get_val()))
						min_bounds.set_lower(candidate);
				}
				else
					min_bounds.set_lower(candidate);

				//change the upper bound
				if(min_bounds.upper().is_finite()){
					if(my_comp.is_definitely_strictly_smaller(s3.f_lambda, min_bounds.upper().get_val()))
						min_bounds.set_upper(s3.f_lambda);
				}
				else
					min_bounds.set_upper(s3.f_lambda);

				return min_bounds;
			}
			else if(my_comp.is_maybe_equal(sample.f_lambda, s2.f_lambda)){
//			else if(sample.f_lambda == s2.f_lambda){
				prob_state = stop; // minimum reached state
				min_bounds.set_lower(sample.f_lambda);
				min_bounds.set_upper(sample.f_lambda);
				return min_bounds;
			}
			else{
				print_pivots();
				this->plot_graph();
				this->plot_function(-1,0.5,0.001);
				std::cout << "new_sample.x:" << sample.lambda << ", new_sample.fx:" << sample.f_lambda << std::endl;
				throw std::runtime_error("update_interval: by convexity, this state should not be reached!");
			}
		}
		else{}
	}
	// Check if the passed sample lies between s3 and s4.
	else if(my_comp.is_definitely_strictly_larger(sample.lambda, s3.lambda) &&
			my_comp.is_definitely_strictly_smaller(sample.lambda, s4.lambda)){
//		std::cout << "new sample between s3,s4" <<std::endl;
//		std::cout << std::setprecision(20) << "sample.lambda:" << sample.lambda << ", s3.lambda:" << s3.lambda << std::endl;
		if(my_comp.is_maybe_equal(s3.lambda,sample.lambda)){
//			std::cout << "s3-s4, maybe test pass" << std::endl;
		}
		if(prob_state == one || prob_state == four_point){
			if(my_comp.is_definitely_strictly_larger(sample.f_lambda, s3.f_lambda)){
				s5 = s4;
				s4 = sample;
				prob_state = four_point;
				// Change the bounds
				//Upper bound doesn't change

				//Lower bound may change
				bool conv_pts = false;
				intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5,conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(my_comp.is_definitely_strictly_smaller(intersect1.f_lambda,intersect2.f_lambda) )
					candidate = intersect1.f_lambda;
				else if(my_comp.is_definitely_strictly_larger(intersect1.f_lambda, intersect2.f_lambda))
					candidate = intersect2.f_lambda;
				else if(my_comp.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
//				else if(intersect1.f_lambda == s3.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}
				else // choose one of the triangles arbitrarily.
					candidate = intersect1.f_lambda;
				//update lower bound
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(candidate, min_bounds.lower().get_val()))
						min_bounds.set_lower(candidate);
				}
				else
					min_bounds.set_lower(candidate);

				return min_bounds;
			}
			else if(my_comp.is_maybe_equal(sample.f_lambda, s3.f_lambda)){
//			else if(sample.f_lambda == s3.f_lambda){
				s1 = s2 ;
				s2 = s3;
				s3 = sample;
				prob_state = three;
				// Change the bounds
				// Upper bound doesn't change

				//Lower bound may change
				bool conv_pts = false;
				intersect = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(intersect.f_lambda == s2.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}
				// update lower bound
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(intersect.f_lambda, min_bounds.lower().get_val()))
						min_bounds.set_lower(intersect.f_lambda);
				}
				else
					min_bounds.set_lower(intersect.f_lambda);

				return min_bounds;

			}
			else if(my_comp.is_definitely_strictly_smaller(sample.f_lambda, s3.f_lambda)){
				s1 = s2;
				s2 = s3;
				s3 = sample;
				// Change the bounds

				//change the upper bound
				if(min_bounds.upper().is_finite()){
					if(my_comp.is_definitely_strictly_smaller(s3.f_lambda, min_bounds.upper().get_val()))
						min_bounds.set_upper(s3.f_lambda);
				}
				else
					min_bounds.set_upper(s3.f_lambda);
				//change the lower bound if the state is in 4
				if(prob_state == four_point){
					bool conv_pts = false;
					intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4, conv_pts);
					intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5, conv_pts);
					if(conv_pts){// heuristic to speed up computation and to bypass numerical error
						// minima reached heuristics
						prob_state = stop;
						min_bounds.set_lower(s3.f_lambda);
						min_bounds.set_upper(s3.f_lambda);
						return min_bounds;
					}

					if(my_comp.is_definitely_strictly_smaller(intersect1.f_lambda,intersect2.f_lambda) )
						candidate = intersect1.f_lambda;
					else if(my_comp.is_definitely_strictly_larger(intersect1.f_lambda, intersect2.f_lambda))
						candidate = intersect2.f_lambda;
					else if(my_comp.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
//					else if(intersect1.f_lambda == s3.f_lambda){
						prob_state = stop;
						min_bounds.set_lower(s3.f_lambda);
						min_bounds.set_upper(s3.f_lambda);
						return min_bounds;
					}
					else // choose one of the triangles arbitrarily.
						candidate = intersect1.f_lambda;

					if(min_bounds.lower().is_finite()){
						if(my_comp.is_definitely_strictly_larger(candidate, min_bounds.lower().get_val()))
							min_bounds.set_lower(candidate);
					}
					else
						min_bounds.set_lower(candidate);
				}
				else{ // prob_state == 1
					typename lb_search_opt<scalar_type, functor>::sample_type s;
					bool conv_pts = false;
					s.lambda = s4.lambda;
					s.f_lambda = scalar_type(0);
					intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
					intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s,s4,conv_pts);
					if(conv_pts){// heuristic to speed up computation and to bypass numerical error
						// minima reached heuristics
						prob_state = stop;
						min_bounds.set_lower(s3.f_lambda);
						min_bounds.set_upper(s3.f_lambda);
						return min_bounds;
					}

					if(min_bounds.lower().is_finite()){
						min_bounds.set_lower(get_maximum(min_bounds.lower().get_val(),get_minimum(intersect1.f_lambda,intersect2.f_lambda)));
					}
					else{
						min_bounds.set_lower(get_minimum(intersect1.f_lambda,intersect2.f_lambda));
					}
				}
				return min_bounds;
			}
			else{}
		}
		else if(prob_state == two){
			s4 = sample;
			// Update lower bound here
			typename lb_search_opt<scalar_type, functor>::sample_type s;
			bool conv_pts = false;
			s.lambda = s1.lambda;
			s.f_lambda = scalar_type(0);
			intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s,s1,conv_pts);

			if(conv_pts){// heuristic to speed up computation and to bypass computation error
				// minima reached heuristics
				prob_state = stop;
				min_bounds.set_lower(s2.f_lambda);
				min_bounds.set_upper(s2.f_lambda);
				return min_bounds;
			}

			if(min_bounds.lower().is_finite()){
				min_bounds.set_lower(get_maximum(min_bounds.lower().get_val(),get_minimum(intersect1.f_lambda,intersect2.f_lambda)));
			}
			else{
				min_bounds.set_lower(get_minimum(intersect1.f_lambda,intersect2.f_lambda));
			}

		}
		else if(prob_state == three){
			if(my_comp.is_definitely_strictly_larger(sample.f_lambda, s3.f_lambda)){
				s4 = sample;
				// Update the lower bound
				typename lb_search_opt<scalar_type, functor>::sample_type intersect;
				bool conv_pts = false;
				intersect = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(intersect.f_lambda == s2.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s2.f_lambda);
					min_bounds.set_upper(s2.f_lambda);
					return min_bounds;
				}
				//change the lower bound;
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(intersect.f_lambda, min_bounds.lower().get_val()))
						min_bounds.set_lower(intersect.f_lambda);
				}
				else
					min_bounds.set_lower(intersect.f_lambda);

			}
			else if(my_comp.is_maybe_equal(sample.f_lambda, s3.f_lambda)){
//			else if(sample.f_lambda == s3.f_lambda){
				prob_state = stop; // minimum reached state
				min_bounds.set_lower(sample.f_lambda);
				min_bounds.set_upper(sample.f_lambda);
				return min_bounds;
			}
			else{
				print_pivots();
				this->plot_graph();
				this->plot_function(-1,0.5,0.001);
				std::cout << "new_sample.x:" << sample.lambda << ", new_sample.fx:" << sample.f_lambda << std::endl;
				throw std::runtime_error("update_interval: by convexity, this state should not be reached!");
			}
		}
		else{}
	}
	// Check if the passed sample lies between s4 and s5.
	else if(my_comp.is_definitely_strictly_larger(sample.lambda, s4.lambda) &&
			my_comp.is_definitely_strictly_smaller(sample.lambda, s5.lambda)){
//		std::cout << "new sample between s4,s5" <<std::endl;
		if(prob_state!=four_point)
			return min_bounds;
		else{
			if(my_comp.is_definitely_strictly_larger(sample.f_lambda, s4.f_lambda)){
				s5 = sample;
				// Change the bounds
				//Upper bound doesn't change

				//lower bound may change
				bool conv_pts = false;
				intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
				intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5,conv_pts);
				if(conv_pts){// heuristic to speed up computation and to bypass numerical error
					// minima reached heuristics
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}

				if(my_comp.is_definitely_strictly_smaller(intersect1.f_lambda,intersect2.f_lambda) )
					candidate = intersect1.f_lambda;
				else if(my_comp.is_definitely_strictly_larger(intersect1.f_lambda, intersect2.f_lambda))
					candidate = intersect2.f_lambda;
				else if(my_comp.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
//				else if(intersect1.f_lambda == s3.f_lambda){
					prob_state = stop;
					min_bounds.set_lower(s3.f_lambda);
					min_bounds.set_upper(s3.f_lambda);
					return min_bounds;
				}
				else // choose one of the triangles arbitrarily.
					candidate = intersect1.f_lambda;
				//update the lower bound
				if(min_bounds.lower().is_finite()){
					if(my_comp.is_definitely_strictly_larger(candidate, min_bounds.lower().get_val()))
						min_bounds.set_lower(candidate);
				}
				else
					min_bounds.set_lower(candidate);
				return min_bounds;
			}
			else{
				print_pivots();
				this->plot_graph();
				this->plot_function(-1,0.5,0.001);
				std::cout << "new_sample.x:" << sample.lambda << ", new_sample.fx:" << sample.f_lambda << std::endl;
				throw std::runtime_error("update_interval: by convexity, this state should not be reached!");
			}
		}
	}
	else{
		// throw std::runtime_error("new sample out of bounds\n");
	}
	return min_bounds;
}
;

/**
 * This method chooses the next sampling point in the function domain.
 *
 * @return The sampling point.
 */
template<class scalar_type, template<typename > class functor>
scalar_type lb_search_opt<scalar_type,functor>::next_sample(){

	math::numeric::approx_comparator<scalar_type> my_comparator;
	scalar_type request;
	typename lb_search_opt<scalar_type, functor>::sample_type min_1,min_2,s;
	bool conv_pts = false;
	typename lb_search_opt<scalar_type, functor>::interval sample_interval;

//	std::cout << "next_sample:" << std::endl;
//	print_pivots();
	if(prob_state == one){
		s.lambda = s4.lambda;
		s.f_lambda = scalar_type(0);
		min_1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
		min_2 = this->template line_intersection_stable<precise_float>(s2,s3,s,s4,conv_pts);

		if(conv_pts){// heuristic to speed up computation and to bypass computation error
			// minima reached heuristics
			prob_state = stop;
			min_bounds = typename lb_search_opt<scalar_type, functor>::interval(
								scalar_with_infinity<scalar_type>(s3.f_lambda),
								scalar_with_infinity<scalar_type>(s3.f_lambda));
//			std::cout << "STATE 1, requested sample:" << s3.lambda;
			return s3.lambda; // So the the problem goes redundant
		}
		/*
			Choose the subinterval which has a lower minima
		*/
		bool left = false, right = false;
		if(my_comparator.is_definitely_strictly_smaller(min_1.f_lambda, min_2.f_lambda)){
			sample_interval = typename lb_search_opt<scalar_type, functor>::interval(
					scalar_with_infinity<scalar_type>(s2.lambda),
					scalar_with_infinity<scalar_type>(s3.lambda));
			left = true;
		}
		else
		{
			sample_interval = typename lb_search_opt<scalar_type, functor>::interval(
									scalar_with_infinity<scalar_type>(s3.lambda),
									scalar_with_infinity<scalar_type>(s4.lambda));
			right = true;
		}

		request = sample_request(sample_interval,left);
//		std::cout << "STATE 1, requested sample:" << request;
//		std::cout << "next_sample: state 1" << std::endl;
	}
	else if(prob_state == two){
		s.lambda = s1.lambda;
		s.f_lambda = scalar_type(0);

		min_1 = this->template line_intersection_stable<precise_float>(s2,s3,s,s1,conv_pts);
		min_2 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);

		if(conv_pts){// heuristic to speed up computation and to bypass computation error
			// minima reached heuristics
			prob_state = stop;
			min_bounds = typename lb_search_opt<scalar_type, functor>::interval(
								scalar_with_infinity<scalar_type>(s2.f_lambda),
								scalar_with_infinity<scalar_type>(s2.f_lambda));
//			std::cout << "STATE 2, requested sample:" << s2.lambda;
			return s2.lambda; // So the the problem goes redundant
		}

		/*
			Choose the subinterval which has a lower minima
		*/
		bool left = false, right = false;
		if(my_comparator.is_definitely_strictly_smaller(min_1.f_lambda, min_2.f_lambda)){
			sample_interval = typename lb_search_opt<scalar_type, functor>::interval(
					scalar_with_infinity<scalar_type>(s1.lambda),
					scalar_with_infinity<scalar_type>(s2.lambda));
			left = true;
		}
		else
		{
			sample_interval = typename lb_search_opt<scalar_type, functor>::interval(
									scalar_with_infinity<scalar_type>(s2.lambda),
									scalar_with_infinity<scalar_type>(s3.lambda));
			right = true;
		}

		request = sample_request(sample_interval,left);
//		std::cout << "STATE 2, requested sample:" << request;
	}
	else if(prob_state == three){
		bool conv_pts = false;
		typename lb_search_opt<scalar_type, functor>::sample_type intersect = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
		if(conv_pts){// heuristic to speed up computation and to bypass computation error
			// minima reached heuristics
			prob_state = stop;
			min_bounds = typename lb_search_opt<scalar_type, functor>::interval(
								scalar_with_infinity<scalar_type>(s3.f_lambda),
								scalar_with_infinity<scalar_type>(s3.f_lambda));
//			std::cout << "STATE 3, requested sample:" << s3.lambda << std::endl;
			return s3.lambda; // So the the problem goes redundant

		}
		if(my_comparator.is_maybe_equal(intersect.lambda,s3.lambda) ||
				my_comparator.is_maybe_equal(intersect.lambda,s2.lambda)){
			prob_state = stop;
			min_bounds = typename lb_search_opt<scalar_type, functor>::interval(
								scalar_with_infinity<scalar_type>(intersect.f_lambda),
								scalar_with_infinity<scalar_type>(intersect.f_lambda));
//			std::cout << "STATE 3, requested sample:" << intersect.lambda << std::endl;
			return intersect.lambda;
		}
		request = intersect.lambda;
//		std::cout << "STATE 3, requested sample:" << request;
//		std::cout << "next_sample: state 3" << std::endl;
	}
	else if(prob_state == four_point){
		typename lb_search_opt<scalar_type, functor>::sample_type intersect1,intersect2;
		bool conv_pts = false;
		intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
//		std::cout << "next_sample:intersect1.f:" << std::setprecision(15) << intersect1.f_lambda << std::endl;
//		std::cout << "next_sample:intersect1.x:" << std::setprecision(15) << intersect1.lambda << std::endl;
		intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5,conv_pts);
//		std::cout << "next_sample:intersect2.f:" << std::setprecision(15) << intersect2.f_lambda << std::endl;
//		std::cout << "next_sample:intersect2.x:" << std::setprecision(15) << intersect2.lambda << std::endl;

		if(abs(intersect1.lambda - s3.lambda) < 10e-14 || abs(intersect2.lambda - s3.lambda) < 10e-14)
			conv_pts = true;
		if(conv_pts){// heuristic to speed up computation and to bypass computation error
					// minima reached heuristics
			//std::cout << "next_sample: conv pts\n";
			prob_state = stop;
			min_bounds = typename lb_search_opt<scalar_type, functor>::interval(
								scalar_with_infinity<scalar_type>(s3.f_lambda),
								scalar_with_infinity<scalar_type>(s3.f_lambda));
//			std::cout << "STATE 3, requested sample:" << s3.lambda << std::endl;
			return s3.lambda; // So the the problem goes redundant

		}

		if(my_comparator.is_definitely_strictly_smaller(intersect1.f_lambda,intersect2.f_lambda) ){
			request = intersect1.lambda;
			//std::cout << "Left triangle choosen" << std::endl;

		}
		else if (my_comparator.is_definitely_strictly_larger(intersect1.f_lambda , intersect2.f_lambda)){
			request = intersect2.lambda;
			//std::cout << "Right triangle choosen" << std::endl;
		}
		else if(my_comparator.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
//			std::cout << "Inside next_sample: int1 == s3 case\n";
			prob_state = stop;
			min_bounds = typename lb_search_opt<scalar_type, functor>::interval(
								scalar_with_infinity<scalar_type>(s3.f_lambda),
								scalar_with_infinity<scalar_type>(s3.f_lambda));
//			std::cout << "STATE 4, requested sample:" << s3.lambda;
			return s3.lambda; // So the the problem goes redundant

		}
		else // intersect1.f_lambda == intersect2.f_lambda
			request = intersect1.lambda; // choose the left subinterval arbitrarily for further search.
//		std::cout << "STATE 4, requested sample:" << request;
		//debug
//		std::cout << "s3.f_lambda" << s3.f_lambda  << std::endl;
//		std::cout << "next_sample: state 4" << std::endl;
	}
	else if(prob_state == stop){
//		std::cout << "next_sample: stop state" << std::endl;
		throw std::runtime_error("next_sample: No sample request in end state!\n");
	}
	else{
//		std::cout << "next_sample: error state" << std::endl;
		throw std::runtime_error("next_sample: Problem state undefined.\n");
	}
//	std::cout << "next_sample: request:" << request <<std::endl;
	std::cout.flush();
	return request;
}
;
/**
 * Computes the initial 4 sampling points within the passed my_interval.
 *
 * @param my_interval
 * @param t1
 * @param t2
 * @param t3
 * @param t4
 * @return
 */

template<class scalar_type, template<typename > class functor>
bool lb_search_opt<scalar_type,functor>::get_pivots(
		const lb_search_opt<scalar_type, functor>::interval& my_interval,
		lb_search_opt<scalar_type, functor>::sample_type& t1,
		lb_search_opt<scalar_type, functor>::sample_type& t2,
		lb_search_opt<scalar_type, functor>::sample_type& t3,
		lb_search_opt<scalar_type, functor>::sample_type& t4) {


	sample_type s1,s2,s;
	math::numeric::approx_comparator<scalar_type> my_comparator;
	double GOLD = 1.618034;

	if(my_interval.lower().is_infinity() && my_interval.upper().is_infinity())
	{
		LOGGER(DEBUG7,"lb_search_opt:get_pivots[-infty, +infty]","Initiating pivots search");
		// sample at 2 arbitrary points, say 0 and 1.
//		std::cout << "pivots: inf bound on both sides" << std::endl;
		s1 = sample(scalar_type(0));
		s2 = sample(scalar_type(1));
		add_sample(s1.lambda, s1.f_lambda);
		add_sample(s2.lambda, s2.f_lambda);



		if(my_comparator.is_maybe_equal(s1.f_lambda,s2.f_lambda)){
			t1 = s1;
			t4 = s2;
			t2 = sample(scalar_type(1/3));
			add_sample(t2.lambda, t2.f_lambda);
			if(my_comparator.is_maybe_equal(t2.f_lambda,s2.f_lambda)){
				t1 = t2 = t3 = t4 = s2;
				LOGGER(DEBUG7,"lb_search_opt:get_pivots[-infty,+infty]","Pivots search terminated successfully with the minima found");
				return true;

			}
			t3 = sample(scalar_type(2/3));
			add_sample(t3.lambda, t3.f_lambda);

			return false;
		}
		else if(my_comparator.is_definitely_strictly_smaller(s1.f_lambda, s2.f_lambda)){
//			std::cout << "s1.f_lambda" << s1.f_lambda << std::endl;
//			std::cout << "s2.f_lambda" << s2.f_lambda << std::endl;

			t4 = s2;
			t3 = s1;

			// go left of s1.
			//t2 = sample(scalar_type(-10));
			// Downhill decend with an increase with a constant factor

			t2 = sample(t3.lambda - (t4.lambda - t3.lambda)*scalar_type(GOLD));
			add_sample(t2.lambda, t2.f_lambda);
//			std::cout << "t2.lambda:" << t2.lambda << std::endl;
//			std::cout << "t2.f_lambda:" << t2.f_lambda << std::endl;

			s = sample(t2.lambda - (t3.lambda - t2.lambda)*scalar_type(GOLD));
			add_sample(s.lambda, s.f_lambda);

			while(my_comparator.is_definitely_strictly_smaller(s.f_lambda,t2.f_lambda)){
				s1 = s;
				s = sample(s.lambda - (t2.lambda - s.lambda)*scalar_type(GOLD));
				add_sample(s.lambda, s.f_lambda);
				t2 = s1;

				LOGGER(DEBUG7,"lb_search_opt:get_pivots[-infty,+infty]","Searching a pivot point towards -infty, currently at:"+to_string(s.lambda));

			};
			t1 = s;
		}
		else{ // s1.f_lambda > s2.f_lambda
			t1 = s1;
			t2 = s2;

			// go right of s2
			//t3 = sample(scalar_type(20));
			t3 = sample(t2.lambda + (t2.lambda-t1.lambda)*scalar_type(GOLD)); // Downhill decend
			add_sample(t3.lambda, t3.f_lambda);

			// now keep going right until sample is more that t3 sample.
			s = sample(t3.lambda + (t3.lambda - t2.lambda)*scalar_type(GOLD));
			add_sample(s.lambda, s.f_lambda);

			while(my_comparator.is_definitely_strictly_smaller(s.f_lambda,t3.f_lambda)){
				s1 = s;
				s = sample(s.lambda + (s.lambda - t3.lambda)*scalar_type(GOLD));
				add_sample(s.lambda, s.f_lambda);
				t3 = s1;
				LOGGER(DEBUG7,"lb_search_opt:get_pivots[-infty,+infty]","Searching a pivot point towards +infty, currently at:"+to_string(s.lambda));

			};
			t4 = s;
		}
	}
	else if(my_interval.lower().is_infinity() && !my_interval.upper().is_infinity()){
		LOGGER(DEBUG7,"lb_search_opt:get_pivots[-infty,R]","Initiating pivots search");
		// we have finite bound on the right side of the search domain
//		std::cout << "finite bound on the right side of the search domain" << std::endl;
//		std::cout << "upper limit:" << my_interval.upper().get_val() << std::endl;

		s1 = sample(my_interval.upper().get_val());
		s2 = sample(my_interval.upper().get_val()-scalar_type(1)); // Initial points arbitrarily chosen

		add_sample(s1.lambda, s1.f_lambda);
		add_sample(s2.lambda, s2.f_lambda);

		if(my_comparator.is_definitely_strictly_larger(s1.f_lambda, s2.f_lambda)){
			t4 = s1;
			t3 = s2;

			// go to further left
			//t2 = sample(my_interval.upper().get_val()-scalar_type(2*MOVE_PARAM));
			t2 = sample(t3.lambda - (t4.lambda - t3.lambda)*scalar_type(GOLD));
			add_sample(t2.lambda, t2.f_lambda);

			s = sample(t2.lambda - (t3.lambda - t2.lambda)*scalar_type(GOLD));
			add_sample(s.lambda, s.f_lambda);

			while(my_comparator.is_definitely_strictly_smaller(s.f_lambda,t2.f_lambda)){
				s1 = s;
				s = sample(s.lambda - (t2.lambda - s.lambda)*scalar_type(GOLD));
				add_sample(s.lambda, s.f_lambda);
				t2 = s1;

				LOGGER(DEBUG7,"lb_search_opt:get_pivots[-infty,R]","Searching a pivot point towards -infty, currently at:"+to_string(s.lambda));
			};
			t1 = s;
		}
		else if(my_comparator.is_definitely_strictly_smaller(s1.f_lambda, s2.f_lambda )){
			// Notice that we are sure of the min lying between t1 and t4 here.

//			std::cout << "case 2" << std::endl;
			t1 = s2;
			t4 = s1;
			t2 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(2));
			add_sample(t2.lambda, t2.f_lambda);

			// Check for collinearity
			if(is_colinear(t1,t2,t4)){
				t1 = t2 = t3 = t4 = t4;
				LOGGER(DEBUG7,"lb_search_opt:get_pivots[-infty,R]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
				return true;
			}

			// move towards t4 until sample is less than t4 sample.

			s = sample(t2.lambda + (t4.lambda - t2.lambda)/scalar_type(2));
			add_sample(s.lambda, s.f_lambda);

			sample_type s_prev = t2;
			while(math::maybe(my_comparator.is_LE(t4.f_lambda, s.f_lambda))){ // This while loop will get stuck if the min is close to t4.
				if(is_colinear(t4,s,s_prev)){
					t1 = t2 = t3 = t4 = t4;
					LOGGER(DEBUG7,"lb_search_opt:get_pivots[-infty,R]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
					return true;
				}
				s_prev = s;
				s = sample(s.lambda + (t4.lambda - s.lambda)/scalar_type(2));
				add_sample(s.lambda, s.f_lambda);
//				std::cout << "right bounded case: while loop" << std::endl;
//				std::cout << "t4.x" << t1.lambda << "t4.f" << t1.f_lambda << std::endl;
//				std::cout << "s.x:" << s.lambda << ", s.f:" << s.f_lambda << std::endl;
//				std::cout << "get_pivots: while 2\n";
				LOGGER(DEBUG7,"lb_search_opt:get_pivots[-infty,R]","Searching a pivot point towards k, currently at:"+to_string(s.lambda));
			}
			t3 = s;

		}
		else{//s1.f_lambda == s2.f_lambda
			// Notice that we are sure of the min lying between t1 and t4 here.
//			std::cout << "case 3" << std::endl;
			t1 = s2;
			t4 = s1;
			t2 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(3));
			add_sample(t2.lambda, t2.f_lambda);
			t3 = sample(t1.lambda + scalar_type(2)*(t4.lambda - t1.lambda)/scalar_type(3));
			add_sample(t3.lambda, t3.f_lambda);


			if(my_comparator.is_maybe_equal(t2.f_lambda, t4.f_lambda) || my_comparator.is_maybe_equal(t3.f_lambda, t4.f_lambda)){ // meaning the function plot between t1 to t4 is a straight line
				t1 = t2 = t3 = t4;
				LOGGER(DEBUG7,"lb_search_opt:get_pivots[-infty,R]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
				return true;// true return means min is found already
			}
		}
	}
	else if(!my_interval.lower().is_infinity() && my_interval.upper().is_infinity()){
		// we have finite bound on the left side of the search domain

		LOGGER(DEBUG7,"lb_search_opt:get_pivots[L,+infty]","Initiating pivots search");

		s1 = sample(my_interval.lower().get_val());
		add_sample(s1.lambda, s1.f_lambda);
		s2 = sample(my_interval.lower().get_val()+scalar_type(1)); // lets call this value(1) as pivot movement parameter!
		add_sample(s2.lambda, s2.f_lambda);


		if(my_comparator.is_definitely_strictly_larger(s1.f_lambda, s2.f_lambda )){
			t1 = s1;
			t2 = s2;

			// go to further right
			//t3 = sample(my_interval.upper().get_val()+scalar_type(2*10));
			t3 = sample(t2.lambda + (t2.lambda - t1.lambda)*scalar_type(GOLD));
			add_sample(t3.lambda, t3.f_lambda);
			s = sample(t3.lambda + (t3.lambda - t2.lambda)*scalar_type(GOLD));
			add_sample(s.lambda, s.f_lambda);

			while(my_comparator.is_definitely_strictly_smaller(s.f_lambda,t3.f_lambda)){
				s1 = s;
				s = sample(s.lambda + (s.lambda - t3.lambda)*scalar_type(GOLD));
				add_sample(s.lambda, s.f_lambda);
				t3 = s1;

				LOGGER(DEBUG7,"lb_search_opt:get_pivots[L,+infty]","Searching a pivot point towards +infty, currently at:"+to_string(s.lambda));

			};
			t4 = s;
		}
		else if(my_comparator.is_definitely_strictly_smaller(s1.f_lambda, s2.f_lambda)){
			// Notice that we are sure of the min lying between t1 and t4 here.
			t1 = s1;
			t4 = s2;
			t3 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(2));
			add_sample(t3.lambda, t3.f_lambda);
//			std::cout << "t3.x" << t3.lambda << "t3.f" << t3.f_lambda << std::endl;
			// move left towards t1 until sample is less than t1 sample.
			// check for colinearity
			if(is_colinear(t1,t3,t4)){
				LOGGER(DEBUG7,"lb_search_opt:get_pivots[L,+infty]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
				t1 = t2 = t3 = t4 = t1;
				return true;
			}
			else{
				s = sample(t1.lambda + (t3.lambda - t1.lambda)/scalar_type(2));
				add_sample(s.lambda, s.f_lambda);
	//			std::cout << "s.x" << s.lambda << "s.f" << s.f_lambda <<std::endl;
				sample_type s_prev = t3;
				//while(my_comparator.is_LE(t1.f_lambda, s.f_lambda)){// This while loop will get stuck if the min is at or close to t1
				while(math::maybe(my_comparator.is_LE(t1.f_lambda, s.f_lambda))){// This while loop will get stuck if the min is at or close to t1
					if(is_colinear(t1,s,s_prev)){
						//my_comparator.is_definitely_strictly_smaller(t1.f_lambda,s3.f_lambda)
	//					std::cout << "t1.f = " << t1.f_lambda << ", s.f =" << s.f_lambda << ", s_prev.f=" << s_prev.f_lambda << std::endl;
						//if(my_comparator.is_LE(t1.f_lambda, s.f_lambda) && my_comparator.is_LE(s.f_lambda, s_prev.f_lambda)){
						//todo: is_LE seems not working
						// min found
	//						std::cout << "left bounded case: min found" << std::endl;
						t1 = t2 = t3 = t4 = t1;
						LOGGER(DEBUG7,"lb_search_opt:get_pivots[k,+infty]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
						return true;
					}
					// Another check to see if the min is at t1
					if(my_comparator.is_maybe_equal(t1.f_lambda,s.f_lambda)){
						t1 = t2 = t3 = t4 = t1;
	//					std::cout << "min reached" << std::endl;
						LOGGER(DEBUG7,"lb_search_opt:get_pivots[k,+infty]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
						return true;
					}

					s_prev = s;
					s = sample(t1.lambda +(s.lambda - t1.lambda)/scalar_type(2));
					add_sample(s.lambda, s.f_lambda);
	//				std::cout << "left bounded case: while loop" << std::endl;
	//				std::cout << "t1.x" << t1.lambda << "t1.f" << t1.f_lambda << std::endl;
	//				std::cout << "s.x:" << s.lambda << ", s.f:" << s.f_lambda << std::endl;
				}
				t2 = s;
			}
		}
		else{
			// Notice that we are sure of the min lying between t1 and t4 here.
			t1 = s1;
			t4 = s2;
			t2 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(3));
			add_sample(t2.lambda, t2.f_lambda);
			t3 = sample(t1.lambda + scalar_type(2)*(t4.lambda - t1.lambda)/scalar_type(3));
			add_sample(t3.lambda, t3.f_lambda);

			if(my_comparator.is_maybe_equal(t2.f_lambda, t4.f_lambda) || my_comparator.is_maybe_equal(t3.f_lambda, t4.f_lambda)){
				// meaning the function plot between t1 to t4 is a straight line
				t1 = t2 = t3 = t4;
				LOGGER(DEBUG7,"lb_search_opt:get_pivots[k,+infty]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
				return true;// true return means min is found already
			}
		}
	}
	else{ //we have finite bounds on the search space for minimum
		LOGGER(DEBUG7,"lb_search_opt:get_pivots[L,R]","Initiating pivots search");
		sample_type s;
		t1 = sample(my_interval.lower().get_val());
		add_sample(t1.lambda, t1.f_lambda);
		t4 = sample(my_interval.upper().get_val());
		add_sample(t4.lambda, t4.f_lambda);

		sample_type s_prev;
		if(my_comparator.is_definitely_strictly_larger(t1.f_lambda, t4.f_lambda)){
			t2 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(2));
			add_sample(t2.lambda, t2.f_lambda);
			if(is_colinear(t1,t2,t4)){
				t1 = t2 = t3 = t4 = t4;
				LOGGER(DEBUG7,"lb_search_opt:get_pivots[k1,k2]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
				return true;
			}
			// move towards t4 until sample is less than t4 sample.

			s = sample(t2.lambda + (t4.lambda - t2.lambda)/scalar_type(2));
			add_sample(s.lambda, s.f_lambda);

			s_prev = t2;
			while(my_comparator.is_LE(t4.f_lambda, s.f_lambda)){ // This while loop will get stuck if the min is at or close to t4.
				s_prev = s;
				s = sample(s.lambda + (t4.lambda - s.lambda)/scalar_type(2));
				add_sample(s.lambda, s.f_lambda);

				if(is_colinear(t4,s,s_prev)){
					t1 = t2 = t3 = t4 = t4;
					LOGGER(DEBUG7,"lb_search_opt:get_pivots[k1,k2]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
					return true;
				}

			}
			t3 = s;
		}
		else if(my_comparator.is_definitely_strictly_smaller(t1.f_lambda, t4.f_lambda)){
			t3 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(2));
			add_sample(t3.lambda, t3.f_lambda);
			if(is_colinear(t1,t3,t4)){
				t1 = t2 = t3 = t4 = t1;
				LOGGER(DEBUG7,"lb_search_opt:get_pivots[k1,k2]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
				return true;
			}
			// move left towards t1 until sample is less than t1 sample.
			s = sample(t1.lambda + (t3.lambda - t1.lambda)/scalar_type(2));
			add_sample(s.lambda, s.f_lambda);

			s_prev = t3;
			while(math::maybe(my_comparator.is_LE(t1.f_lambda, s.f_lambda))){// This while loop will get stuck if the min is at or close to t1
				s_prev = s;
				s = sample(t1.lambda +(s.lambda - t1.lambda)/scalar_type(2));
				add_sample(s.lambda, s.f_lambda);

				if(is_colinear(t1,s,s_prev)){
					LOGGER(DEBUG7,"lb_search_opt:get_pivots[k1,k2]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
					t4 = t3 = t2 = t1;
					return true;
				}

//				std::cout << "get_pivots: while 2\n";
			}
			t2 = s;
		}
		else{//(t1.f_lambda == t4.f_lambda)
			t2 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(3));
			add_sample(t2.lambda, t2.f_lambda);
			t3 = sample(t1.lambda + scalar_type(2)*(t4.lambda - t1.lambda)/scalar_type(3));
			add_sample(t3.lambda, t3.f_lambda);

			if(my_comparator.is_maybe_equal(t2.f_lambda,t4.f_lambda) || my_comparator.is_maybe_equal(t3.f_lambda, t4.f_lambda)){ // meaning the function plot between t1 to t4 is a straight line
				t1 = t2 = t3 = t4;
				LOGGER(DEBUG7,"lb_search_opt:get_pivots[k1,k2]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
				return true;// true return means min is found already

			}
		}
	}
	return false;
}

template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type, functor>::min_interval lb_search_opt<scalar_type,functor>::section_search
		(lb_search_opt<scalar_type, functor>::sample_type& s1,
		 lb_search_opt<scalar_type, functor>::sample_type& s2,
		 lb_search_opt<scalar_type, functor>::sample_type& s3,
		 lb_search_opt<scalar_type, functor>::sample_type& s4){

	/* CURRENT IMPLEMENTATION HAS BUGS */
	/**
	 * A 4 pivot search algorithm, specialized to merge to exact minimum
	 * for pwa functions, on finite number of iterations.
	 *
	 * Assumption is: f(s1.lambda > f(s2.lambda) and f(s3.lambda) < f(s4.lambda))
	 */

	unsigned int STATE = 0;
	scalar_type minimum, min_intv;

	scalar_with_infinity<scalar_type> lower_bound(scalar_with_infinity<scalar_type>::neg_infty());
	scalar_with_infinity<scalar_type> upper_bound(scalar_with_infinity<scalar_type>::pos_infty());
	bool minimum_found = false,upper_bound_initialised = false, lower_bound_initialised = false;
	typename lb_search_opt<scalar_type, functor>::sample_type new_sample, s5, intersect; // sample for the 5 point case

	math::numeric::approx_comparator<scalar_type> my_comparator;
	while(1){
		if(STATE == 0 && my_comparator.is_definitely_strictly_larger(s2.f_lambda,s3.f_lambda)){ // case 1
			//some renaming
//			std::cout << "STATE 0: 1\n";
			STATE = 1;
		}
		else if(STATE == 0 && my_comparator.is_definitely_strictly_smaller(s2.f_lambda,s3.f_lambda)){ // case 2
//			std::cout << "STATE 0: 2\n";
			STATE = 2;
		}
		else if( STATE == 0 && my_comparator.is_maybe_equal(s2.f_lambda,s3.f_lambda)){
//			std::cout << "STATE 0: 3\n";
			STATE = 3;
		}
		else if(STATE == 1){
			typename lb_search_opt<scalar_type, functor>::interval my_interval(
					scalar_with_infinity<scalar_type>(s3.lambda),
					scalar_with_infinity<scalar_type>(s4.lambda));

			/* code for maintaining the bounds on min */

			if(!upper_bound_initialised){
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(s2.f_lambda,s3.f_lambda));
				upper_bound_initialised = true;
			}
			else{
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(),get_minimum(s2.f_lambda,s3.f_lambda)));
			}
			/*----*/
			bool left = false;
			new_sample = sample_interval(my_interval,left);

			if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda,s3.f_lambda)){
				s1 = s2;
				s2 = s3;
				s3 = new_sample;
				STATE = 1;
//				std::cout << "STATE 1: 1\n";
			}
			else if(my_comparator.is_maybe_equal(new_sample.f_lambda,s3.f_lambda)){
				//renaming
				s1 = s2;
				s2 = s3;
				s3 = new_sample;
				STATE = 3;
//				std::cout << "STATE 1: 3\n";
			}
			else{ // new_sample.f_lambda > s3.f_lambda
				//renaming
				s5 = s4;
				s4 = new_sample;
				STATE = 4; // The 5 point state
//				std::cout << "STATE 1: 4\n";
			}
		}
		else if(STATE == 2){
			typename lb_search_opt<scalar_type, functor>::interval my_interval(
					scalar_with_infinity<scalar_type>(s1.lambda),
					scalar_with_infinity<scalar_type>(s2.lambda));
			/* code for maintaining the bounds on min*/

			if(!upper_bound_initialised)
				//upper_bound = scalar_with_infinity<scalar_type>(lb_search_opt<scalar_type, functor>::get_minimum(s2.f_lambda,s3.f_lambda));
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(s2.f_lambda,s3.f_lambda));
			else{
//				upper_bound = scalar_with_infinity<scalar_type>(
//						lb_search_opt<scalar_type, functor>::get_minimum(upper_bound.get_val(),
//						lb_search_opt<scalar_type, functor>::get_minimum(s2.f_lambda,s3.f_lambda)));
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(),get_minimum(s2.f_lambda,s3.f_lambda)));
			}
			/*----*/
			bool left = true;
			new_sample = sample_interval(my_interval,true);

			if(my_comparator.is_definitely_strictly_larger(new_sample.f_lambda, s2.f_lambda)){
				STATE = 4; // The 5 point state.
				//renaming
				s5 = s4;
				s4 = s3;
				s3 = s2;
				s2 = new_sample;
//				std::cout << "STATE 2: 4\n";
			}
			else if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda , s2.f_lambda)){
				//renaming
				s4 = s3;
				s3 = s2;
				s2 = new_sample;
				STATE = 2;
//				std::cout << "STATE 2: 2\n";
			}
			else{ // new_sample.f_lambda == s2.f_lambda
				s4 = s3;
				s3 = s2;
				s2 = new_sample;
				STATE = 3;
//				std::cout << "STATE 2: 3\n";
			}

		}
		else if(STATE == 3){
			/* code for maintaining the bounds on min*/

			if(!upper_bound_initialised)
				upper_bound = scalar_with_infinity<scalar_type>(s2.f_lambda);
			else{
//				upper_bound = scalar_with_infinity<scalar_type>(
//				lb_search_opt<scalar_type, functor>::get_minimum(upper_bound.get_val(), s2.f_lambda));
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(), s2.f_lambda));
			}
			/*----*/

			bool conv_pts = false;
			intersect = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			new_sample = sample(intersect.lambda);

			if(my_comparator.is_maybe_equal(new_sample.f_lambda,s2.f_lambda)){
				STATE = 6; //stop state
				minimum_found = true;
				min_intv = s2.lambda;
				minimum = s2.f_lambda;
//				std::cout << "STATE 3: 6\n";
//				std::cout << "minimum at state 3:" << minimum << std::endl;
			}
			else if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda , s2.f_lambda)){
				//renaming
				s5 = s4;
				s4 = s3;
				s3 = new_sample;
				STATE = 4; // The 5 point state.
//				std::cout << "STATE 3: 4\n";
			}
			else{
				std::cout << "Something is wrong. Should not reach here logically\n";
			}

		}
/*
		else if(STATE == 4 && true){ // fishy case
			lb_search_opt<scalar_type, functor>::sample_type new_sample1, new_sample2, intersect1, intersect2;
			intersect1 = line_intersect(s1,s2,s3,s4);
			intersect2 = line_intersect(s2,s3,s4,s5);

//			new_sample1 = sample(intersect1.lambda);
//			new_sample2 = sample(intersect2.lambda);

			new_sample = intersect1.f_lambda < intersect2.f_lambda? sample(intersect1.lambda): sample(intersect2.lambda);

			if(new_sampe.f_lambda < s3.f_lambda){
				//renaming
				s5 = s4;
				s4 = s3;
				s3 = new_sample;
				STATE = 4;
			}
			else if(new_sample.f_lambda > s3.f_lambda){
				//renaming
				s1 = s2;
				s2 = new_sample;
				STATE = 4;
			}
			else if(new_sample.f_lambda == s3.f_lambda){

			}

		}
*/
		else if(STATE == 4){ // fishy case
			typename lb_search_opt<scalar_type, functor>::sample_type intersect1, intersect2;

			/* code for maintaining the bounds on min*/

			if(!upper_bound_initialised)
				//upper_bound = scalar_with_infinity<scalar_type>(lb_search_opt<scalar_type, functor>::get_minimum(s3.f_lambda,s4.f_lambda));
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(s3.f_lambda,s4.f_lambda));
			else{
				//upper_bound = scalar_with_infinity<scalar_type>(lb_search_opt<scalar_type, functor>::get_minimum(upper_bound.get_val(),
				//												lb_search_opt<scalar_type, functor>::get_minimum(s3.f_lambda,s4.f_lambda)));
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(),
																				get_minimum(s3.f_lambda,s4.f_lambda)));
			}

			/*----*/

			bool left_triangle = false, right_triangle = false;
			/*DEBUG */

			bool conv_pts = false;
			intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5,conv_pts);

			/*DEBUG*/
//			std::cout << "intersect1.x:" << intersect.lambda <<", y:" << intersect1.f_lambda << std::endl;
//			std::cout << "intersect2.x:" << intersect.lambda <<", y:" << intersect2.f_lambda << std::endl;

			if (my_comparator.is_definitely_strictly_smaller(intersect1.f_lambda , intersect2.f_lambda)) {
				new_sample = sample(intersect1.lambda);

				left_triangle = true;
				/* Code for maintaining the minimum */
				if(!lower_bound_initialised)
					lower_bound = scalar_with_infinity<scalar_type>(intersect1.f_lambda);
				else{
					if(my_comparator.is_definitely_strictly_larger(lower_bound.get_val(), intersect1.f_lambda))
							lower_bound = scalar_with_infinity<scalar_type>(intersect1.f_lambda);
				}
			}
			else if (my_comparator.is_definitely_strictly_larger(intersect1.f_lambda , intersect2.f_lambda)){
				new_sample = sample(intersect2.lambda);

				right_triangle = true;
				/* Code for maintaining the minimum */
				if(!lower_bound_initialised)
					lower_bound = scalar_with_infinity<scalar_type>(intersect2.f_lambda);
				else{
					if(my_comparator.is_definitely_strictly_larger(lower_bound.get_val(), intersect2.f_lambda))
							lower_bound = scalar_with_infinity<scalar_type>(intersect2.f_lambda);
				}
			}
			else{ // intersect1.f_lambda == intersect2.f_lambda
//				std::cout << "STATE 4: 2 triangles are at same height";
				if(my_comparator.is_maybe_equal(intersect1.f_lambda, s3.f_lambda)){
//				if(intersect1.f_lambda == s3.f_lambda){

					// I claim that this means the minimum is indeed s3.f_lambda.

					STATE = 6; //stop state
					minimum_found = true;
					min_intv = s3.lambda;
					minimum = s3.f_lambda;
//					std::cout << "STATE 4: 6\n";
					break;
				}
				//choose a triangle randomly.
				else{
					new_sample = sample(intersect1.lambda); // left triangle choosen arbitrarily.

					left_triangle = true;
					/* Code for maintaining the minimum */
					if(!lower_bound_initialised)
						lower_bound = scalar_with_infinity<scalar_type>(intersect1.f_lambda);
					else{
						if(my_comparator.is_definitely_strictly_larger(lower_bound.get_val(), intersect1.f_lambda))
								lower_bound = scalar_with_infinity<scalar_type>(intersect1.f_lambda);
					}
				}
			}

			if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda , s3.f_lambda)){
				//renaming
				if(left_triangle){
					s5 = s4;
					s4 = s3;
					s3 = new_sample;
					STATE = 4;
//					std::cout << "LTR:sample < s3.f_lambda: STATE 4: 4\n";
				}
				else{ // right_triangle
				//renaming
					s1 = s2;
					s2 = s3;
					s3 = new_sample;
					STATE = 4;
//					std::cout << "RTR:sample < s3.f_lambda: STATE 4: 4\n";
				}
			}
			else if(my_comparator.is_definitely_strictly_larger(new_sample.f_lambda , s3.f_lambda)){
				//renaming
				if(left_triangle){
					s1 = s2;
					s2 = new_sample;
					STATE = 4;
//					std::cout << "LTR:sample > s3.f_lambda: STATE 4: 4\n";
				}
				else{//right triangle
					s5 = s4;
					s4 = new_sample;
					STATE = 4;
//					std::cout << "RTR:sample > s3.f_lambda: STATE 4: 4\n";
//					std::cout << "new_sample.x:" << new_sample.lambda <<" and s3.x:" << s3.lambda << std::endl;
//					std::cout << "new_sample.y:" << new_sample.f_lambda <<" and s3.y:" << s3.f_lambda << std::endl;
//					break;
				}
			}
			else{ //new_sample.f_lambda == s3.f_lambda) // this case needs thought

				// I claim that this means the minimum is indeed s3.f_lambda.

				STATE = 6; //stop state
				minimum_found = true;
				min_intv = s3.lambda;
				minimum = s3.f_lambda;
//				std::cout << "STATE 4: 6\n";
//				std::cout << "minimum at state 3:" << minimum << std::endl;

				/*
				if(left_triangle){
					//renaming
					s1 = s2;
					s2 = new_sample;
					STATE = 3;
					std::cout << "STATE 4: 3\n";
				}
				else{
					//renaming
					s1 = s2;
					s2 = s3;
					s3 = new_sample;
					STATE = 3;
					std::cout << "STATE 4: 3\n";
				} */
			}

		}

		else if(STATE == 6){ //stop state.
//			std::cout << "STATE 6\n";
			break;

		}
		else{}
	}//end of while
	if(minimum_found){
		typename math::numeric::interval<scalar_type> my_y_interval = math::numeric::interval<scalar_type>(scalar_with_infinity<scalar_type>(minimum),scalar_with_infinity<scalar_type>(minimum));
		typename math::numeric::interval<scalar_type> my_x_interval = math::numeric::interval<scalar_type>(
				scalar_with_infinity<scalar_type>(min_intv),scalar_with_infinity<scalar_type>(min_intv));

		typename convex_opt<scalar_type, functor>::min_interval min_structure;
		min_structure.lambda_interval = my_x_interval;
		min_structure.f_lambda_interval = my_y_interval;
//		std::cout << "NUMBER OF ITERATIONS:" << no_samples << std::endl;
		return min_structure;
	}
	else{
		typename lb_search_opt<scalar_type, functor>::interval my_x_interval(
						scalar_with_infinity<scalar_type>(s1.lambda),scalar_with_infinity<scalar_type>(s4.lambda));

		typename lb_search_opt<scalar_type, functor>::interval my_y_interval(lower_bound,upper_bound);

		typename convex_opt<scalar_type, functor>::min_interval min_structure = {my_x_interval, my_y_interval};
//		std::cout << "NUMBER OF ITERATIONS:" << no_samples << std::endl;
		return min_structure;
	}

}
;
template<class scalar_type, template<typename > class functor >
typename convex_opt<scalar_type, functor>::min_interval lb_search_opt<scalar_type,functor>::section_search(){
	typename math::numeric::interval<scalar_type> my_interval(scalar_with_infinity<scalar_type>::neg_infty(),
																  scalar_with_infinity<scalar_type>::pos_infty());
	return section_search(my_interval);
}
;
template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type, functor>::min_interval lb_search_opt<scalar_type,functor>::section_search(
		const math::numeric::interval<scalar_type>& my_interval) {
	sample_type t1,t2,t3,t4;
	typename convex_opt<scalar_type, functor>::min_interval min_structure;
	bool min_found_flag = false;

	if(this->get_minbrak_type() == "gold_desc"){
		//min_found_flag = get_pivots(my_interval, t1,t2,t3,t4);
		minbrak_simple(t1,t2,t3,t4);
	}
	else if(this->get_minbrak_type() == "parab_desc"){
		minbrak_para_ext(t1,t2,t3,t4);
	}
	else{
		throw std::runtime_error("lb_search_opt: minima bracketing method not known");
	}
	if(min_found_flag){ // min is already found;
//		std::cout << "min found after pivot selection\n";
		typename lb_search_opt<scalar_type,functor>::interval my_x_interval(scalar_with_infinity<scalar_type>(t1.lambda),
																 	        scalar_with_infinity<scalar_type>(t2.lambda));
		typename lb_search_opt<scalar_type,functor>::interval my_y_interval(scalar_with_infinity<scalar_type>(t1.f_lambda),
																		 	scalar_with_infinity<scalar_type>(t2.f_lambda));

		// GF: The following creates a compilation error
		//min_structure = typename convex_opt<scalar_type, functor>::min_interval {my_x_interval,my_y_interval};
		min_structure.lambda_interval = my_x_interval;
		min_structure.f_lambda_interval = my_y_interval;
		return  min_structure;
	}
	min_structure = section_search(t1,t2,t3,t4);
	return min_structure;
}
;

template<class scalar_type, template<typename >class functor>
typename lb_search_opt<scalar_type, functor>::interval lb_search_opt<scalar_type, functor>::section_search_opt(){
	math::numeric::approx_comparator<scalar_type> my_comp;
	scalar_type lambda,f_lambda;
	typename lb_search_opt<scalar_type,sup_functor>::sample_type new_sample;
	scalar_type eps = this->get_interval_tolerance();

	while(!min_bounds.is_finite() || (my_comp.is_definitely_strictly_larger(min_bounds.upper().get_val() - min_bounds.lower().get_val(), eps))){

		lambda = next_sample();
//		std::cout << "Requested next sample:" << lambda << std::endl;
		new_sample = sample(lambda);
		add_sample(new_sample.lambda, new_sample.f_lambda);
		min_bounds = update_bounds(new_sample);
//		std::cout << "min bounds after updation: " << min_bounds << std::endl;
//		std::cout << "Prb currently in state: " << get_problem_state() << std::endl;
		//Comparison purpose
		//if(this->get_size()==6)
		//	break;
	}
	return min_bounds;
}

template<class scalar_type, template<typename >class functor>
bool lb_search_opt<scalar_type, functor>::is_min_found(){
	math::numeric::approx_comparator<scalar_type> comp;
	//check s1,s2,s3
	if(comp.is_maybe_equal(s1.f_lambda,s2.f_lambda) && comp.is_maybe_equal(s2.f_lambda,s3.f_lambda)){
		min_bounds.set_lower(s1.f_lambda);
		min_bounds.set_upper(s1.f_lambda);
		return true;
	}
	//check s2,s3,s4
	if(comp.is_maybe_equal(s2.f_lambda,s3.f_lambda) && comp.is_maybe_equal(s3.f_lambda,s4.f_lambda)){
		min_bounds.set_lower(s2.f_lambda);
		min_bounds.set_upper(s2.f_lambda);
		return true;
	}
	//check s1,s3,s4
	if(comp.is_maybe_equal(s1.f_lambda,s3.f_lambda) && comp.is_maybe_equal(s3.f_lambda,s4.f_lambda)){
		min_bounds.set_lower(s1.f_lambda);
		min_bounds.set_upper(s1.f_lambda);
		return true;
	}
	//check s1,s2,s4
	if(comp.is_maybe_equal(s1.f_lambda,s2.f_lambda) && comp.is_maybe_equal(s2.f_lambda,s4.f_lambda)){
		min_bounds.set_lower(s1.f_lambda);
		min_bounds.set_upper(s1.f_lambda);
		return true;
	}
	return false;
}
template<class scalar_type, template<typename > class functor>
void lb_search_opt<scalar_type, functor>::print_pivots(){
	std::cout << "prob_state:" << get_problem_state() <<std::endl ;
	std::cout << "sample points :" << std::endl;
	std::cout << std::setprecision(16);
	std::cout << "s1.x:" << s1.lambda << ", s1.y:" << s1.f_lambda <<std::endl;
	std::cout << "s2.x:" << s2.lambda << ", s2.y:" << s2.f_lambda <<std::endl;
	std::cout << "s3.x:" << s3.lambda << ", s3.y:" << s3.f_lambda <<std::endl;
	std::cout << "s4.x:" << s4.lambda << ", s4.y:" << s4.f_lambda <<std::endl;
	std::cout << "s5.x:" << s5.lambda << ", s5.y:" << s5.f_lambda <<std::endl;
}

}// end of namespace support_function
}// end of namespace continuous

#endif
