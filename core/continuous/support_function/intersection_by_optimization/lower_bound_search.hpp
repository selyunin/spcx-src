/*
 * lower_bound_search.cpp
 *
 *  Created on: Feb 21, 2010
 *      Author: ray
 */


namespace continuous {
	namespace support_function {

	template<class scalar_type, template<typename > class functor>
	void lower_bound_search<scalar_type, functor>::print_pivots(
			const sample_type& s1, const sample_type& s2,
			const sample_type& s3, const sample_type& s4,
			const sample_type& s5){
		std::cout << "sample points :" << std::endl;
		std::cout << std::setprecision(22);
		std::cout << "s1.x:" << s1.lambda << ", s1.y:" << s1.f_lambda <<std::endl;
		std::cout << "s2.x:" << s2.lambda << ", s2.y:" << s2.f_lambda <<std::endl;
		std::cout << "s3.x:" << s3.lambda << ", s3.y:" << s3.f_lambda <<std::endl;
		std::cout << "s4.x:" << s4.lambda << ", s4.y:" << s4.f_lambda <<std::endl;
		std::cout << "s5.x:" << s4.lambda << ", s5.y:" << s5.f_lambda <<std::endl;

	}

	template<class scalar_type, template<typename > class functor>
	bool lower_bound_search<scalar_type,functor>::get_pivots(
			const lower_bound_search<scalar_type, functor>::interval& my_interval,
			lower_bound_search<scalar_type, functor>::sample_type& t1,
			lower_bound_search<scalar_type, functor>::sample_type& t2,
			lower_bound_search<scalar_type, functor>::sample_type& t3,
			lower_bound_search<scalar_type, functor>::sample_type& t4,
			const unsigned int bound = 0) {

		sample_type s1,s2,s;
		math::numeric::approx_comparator<scalar_type> my_comparator;
		double GOLD = 1.618034;
		if(my_interval.lower().is_infinity() && my_interval.upper().is_infinity())
		{
			LOGGER(DEBUG7,"lower_bound_search:get_pivots[-infty, +infty]","Initiating pivots search");
			// sample at 2 arbitrary points, say 0 and 1.
	//		std::cout << "pivots: inf bound on both sides" << std::endl;
			s1 = sample(scalar_type(0));
			if(this->get_dynamics_map().is_empty())
				s2 = sample(scalar_type(1));
			else{
#ifdef OPT_INIT_SAMPLE__
				std::list<scalar_type> lambdas;
				lambdas = this->get_init_sample();
				s2 = sample(lambdas.front());
#else
				s2 = sample(scalar_type(1));
#endif
			}

			add_sample(s1.lambda, s1.f_lambda);
			add_sample(s2.lambda, s2.f_lambda);


			if(my_comparator.is_maybe_equal(s1.f_lambda,s2.f_lambda)){
				t1 = s1;
				t4 = s2;
				if(bound!=0 && this->get_size() >= bound)
					return false;
				t2 = sample(scalar_type(10/3));
				add_sample(t2.lambda, t2.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

				if(my_comparator.is_maybe_equal(t2.f_lambda,s2.f_lambda)){
					t1 = t2 = t3 = t4 = s2;
					LOGGER(DEBUG7,"lower_bound_search:get_pivots[-infty,+infty]","Pivots search terminated successfully with the minima found");

					return true;

				}
				t3 = sample(scalar_type(20/3));
				add_sample(t3.lambda, t3.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

				return false;
			}
			else if(my_comparator.is_definitely_strictly_smaller(s1.f_lambda, s2.f_lambda)){
	//			std::cout << "s1.f_lambda" << s1.f_lambda << std::endl;
	//			std::cout << "s2.f_lambda" << s2.f_lambda << std::endl;

				t4 = s2;
				t3 = s1;
				if(bound!=0 && this->get_size() >= bound)
					return false;

				// go left of s1.
				//t2 = sample(scalar_type(-10));
				// Downhill decend with an increase with a constant factor

				t2 = sample(t3.lambda - (t4.lambda - t3.lambda)*scalar_type(GOLD));
				add_sample(t2.lambda, t2.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;


	//			std::cout << "t2.lambda:" << t2.lambda << std::endl;
	//			std::cout << "t2.f_lambda:" << t2.f_lambda << std::endl;

				s = sample(t2.lambda - (t3.lambda - t2.lambda)*scalar_type(GOLD));
				add_sample(s.lambda, s.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;


				while(my_comparator.is_definitely_strictly_smaller(s.f_lambda,t2.f_lambda)){
					s1 = s;
					s = sample(s.lambda - (t2.lambda - s.lambda)*scalar_type(GOLD));
					add_sample(s.lambda, s.f_lambda);
					t2 = s1;
					if(bound!=0 && this->get_size() >= bound)
						return false;


					LOGGER(DEBUG7,"lower_bound_search:get_pivots[-infty,+infty]","Searching a pivot point towards -infty, currently at:"+to_string(s.lambda));

				};
				t1 = s;
			}
			else{ // s1.f_lambda > s2.f_lambda
				t1 = s1;
				t2 = s2;

				if(bound!=0 && this->get_size() >= bound)
					return false;

				// go right of s2
				//t3 = sample(scalar_type(20));
				t3 = sample(t2.lambda + (t2.lambda-t1.lambda)*scalar_type(GOLD)); // Downhill decend
				add_sample(t3.lambda, t3.f_lambda);

				if(bound!=0 && this->get_size() >= bound)
					return false;

				// now keep going right until sample is more that t3 sample.
				s = sample(t3.lambda + (t3.lambda - t2.lambda)*scalar_type(GOLD));
				add_sample(s.lambda, s.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;


				while(my_comparator.is_definitely_strictly_smaller(s.f_lambda,t3.f_lambda)){
					s1 = s;
					s = sample(s.lambda + (s.lambda - t3.lambda)*scalar_type(GOLD));
					add_sample(s.lambda, s.f_lambda);
					if(bound!=0 && this->get_size() >= bound)
						return false;

					t3 = s1;

					LOGGER(DEBUG7,"lower_bound_search:get_pivots[-infty,+infty]","Searching a pivot point towards +infty, currently at:"+to_string(s.lambda));

				};
				t4 = s;
			}
		}
		else if(my_interval.lower().is_infinity() && !my_interval.upper().is_infinity()){
			LOGGER(DEBUG7,"lower_bound_search:get_pivots[-infty,R]","Initiating pivots search");
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
				if(bound!=0 && this->get_size() >= bound)
					return false;

				// go to further left
				//t2 = sample(my_interval.upper().get_val()-scalar_type(2*MOVE_PARAM));
				t2 = sample(t3.lambda - (t4.lambda - t3.lambda)*scalar_type(GOLD));
				add_sample(t2.lambda, t2.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;


				s = sample(t2.lambda - (t3.lambda - t2.lambda)*scalar_type(GOLD));
				add_sample(s.lambda, s.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;


				while(my_comparator.is_definitely_strictly_smaller(s.f_lambda,t2.f_lambda)){
					s1 = s;
					s = sample(s.lambda - (t2.lambda - s.lambda)*scalar_type(GOLD));
					add_sample(s.lambda, s.f_lambda);
					if(bound!=0 && this->get_size() >= bound)
						return false;

					t2 = s1;


					LOGGER(DEBUG7,"lower_bound_search:get_pivots[-infty,R]","Searching a pivot point towards -infty, currently at:"+to_string(s.lambda));
				};
				t1 = s;
			}
			else if(my_comparator.is_definitely_strictly_smaller(s1.f_lambda, s2.f_lambda )){
				// Notice that we are sure of the min lying between t1 and t4 here.

	//			std::cout << "case 2" << std::endl;
				t1 = s2;
				t4 = s1;
				if(bound!=0 && this->get_size() >= bound)
					return false;

				t2 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(2));
				add_sample(t2.lambda, t2.f_lambda);

				if(bound!=0 && this->get_size() >= bound)
					return false;

				// Check for collinearity
				if(is_colinear(t1,t2,t4)){
					t1 = t2 = t3 = t4 = t4;
					LOGGER(DEBUG7,"lower_bound_search:get_pivots[-infty,R]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
					return true;
				}

				// move towards t4 until sample is less than t4 sample.

				s = sample(t2.lambda + (t4.lambda - t2.lambda)/scalar_type(2));
				add_sample(s.lambda, s.f_lambda);

				if(bound!=0 && this->get_size() >= bound)
					return false;

				sample_type s_prev = t2;
				while(math::maybe(my_comparator.is_LE(t4.f_lambda, s.f_lambda))){ // This while loop will get stuck if the min is close to t4.
					if(is_colinear(t4,s,s_prev)){
						t1 = t2 = t3 = t4 = t4;
						LOGGER(DEBUG7,"lower_bound_search:get_pivots[-infty,R]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
						return true;
					}
					s_prev = s;
					s = sample(s.lambda + (t4.lambda - s.lambda)/scalar_type(2));
					add_sample(s.lambda, s.f_lambda);

					if(bound!=0 && this->get_size() >= bound)
						return false;

	//				std::cout << "right bounded case: while loop" << std::endl;
	//				std::cout << "t4.x" << t1.lambda << "t4.f" << t1.f_lambda << std::endl;
	//				std::cout << "s.x:" << s.lambda << ", s.f:" << s.f_lambda << std::endl;
	//				std::cout << "get_pivots: while 2\n";
					LOGGER(DEBUG7,"lower_bound_search:get_pivots[-infty,R]","Searching a pivot point towards k, currently at:"+to_string(s.lambda));
				}
				t3 = s;

			}
			else{//s1.f_lambda == s2.f_lambda
				// Notice that we are sure of the min lying between t1 and t4 here.
	//			std::cout << "case 3" << std::endl;
				t1 = s2;
				t4 = s1;
				if(bound!=0 && this->get_size() >= bound)
					return false;

				t2 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(3));
				add_sample(t2.lambda, t2.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

				t3 = sample(t1.lambda + scalar_type(2)*(t4.lambda - t1.lambda)/scalar_type(3));
				add_sample(t3.lambda, t3.f_lambda);
					return false;


				if(my_comparator.is_maybe_equal(t2.f_lambda, t4.f_lambda) || my_comparator.is_maybe_equal(t3.f_lambda, t4.f_lambda)){ // meaning the function plot between t1 to t4 is a straight line
					t1 = t2 = t3 = t4;
					LOGGER(DEBUG7,"lower_bound_search:get_pivots[-infty,R]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
					return true;// true return means min is found already
				}
			}
		}
		else if(!my_interval.lower().is_infinity() && my_interval.upper().is_infinity()){
			// we have finite bound on the left side of the search domain

			LOGGER(DEBUG7,"lower_bound_search:get_pivots[L,+infty]","Initiating pivots search");
			s1 = sample(my_interval.lower().get_val());
			add_sample(s1.lambda, s1.f_lambda);
			s2 = sample(my_interval.lower().get_val()+scalar_type(1)); // lets call this value(1) as pivot movement parameter!
			add_sample(s2.lambda, s2.f_lambda);


			if(my_comparator.is_definitely_strictly_larger(s1.f_lambda, s2.f_lambda )){
				t1 = s1;
				t2 = s2;
				if(bound!=0 && this->get_size() >= bound)
					return false;

				// go to further right
				//t3 = sample(my_interval.upper().get_val()+scalar_type(2*10));
				t3 = sample(t2.lambda + (t2.lambda - t1.lambda)*scalar_type(GOLD));
				add_sample(t3.lambda, t3.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;


				s = sample(t3.lambda + (t3.lambda - t2.lambda)*scalar_type(GOLD));
				add_sample(s.lambda, s.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;


				while(my_comparator.is_definitely_strictly_smaller(s.f_lambda,t3.f_lambda)){
					s1 = s;
					s = sample(s.lambda + (s.lambda - t3.lambda)*scalar_type(GOLD));
					add_sample(s.lambda, s.f_lambda);
					if(bound!=0 && this->get_size() >= bound)
						return false;

					t3 = s1;


					LOGGER(DEBUG7,"lower_bound_search:get_pivots[L,+infty]","Searching a pivot point towards +infty, currently at:"+to_string(s.lambda));

				};
				t4 = s;
			}
			else if(my_comparator.is_definitely_strictly_smaller(s1.f_lambda, s2.f_lambda)){
				// Notice that we are sure of the min lying between t1 and t4 here.
				t1 = s1;
				t4 = s2;
				if(bound!=0 && this->get_size() >= bound)
					return false;

				t3 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(2));
				add_sample(t3.lambda, t3.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

	//			std::cout << "t3.x" << t3.lambda << "t3.f" << t3.f_lambda << std::endl;
				// move left towards t1 until sample is less than t1 sample.
				// check for colinearity
				if(is_colinear(t1,t3,t4)){
					LOGGER(DEBUG7,"lower_bound_search:get_pivots[L,+infty]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
					t1 = t2 = t3 = t4 = t1;
					return true;
				}
				else{
					s = sample(t1.lambda + (t3.lambda - t1.lambda)/scalar_type(2));
					add_sample(s.lambda, s.f_lambda);
					if(bound!=0 && this->get_size() >= bound)
						return false;

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
							LOGGER(DEBUG7,"lower_bound_search:get_pivots[k,+infty]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
							return true;
						}
						// Another check to see if the min is at t1
						if(my_comparator.is_maybe_equal(t1.f_lambda,s.f_lambda)){
							t1 = t2 = t3 = t4 = t1;
		//					std::cout << "min reached" << std::endl;
							LOGGER(DEBUG7,"lower_bound_search:get_pivots[k,+infty]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
							return true;
						}

						s_prev = s;
						s = sample(t1.lambda +(s.lambda - t1.lambda)/scalar_type(2));
						add_sample(s.lambda, s.f_lambda);

						if(bound!=0 && this->get_size() >= bound)
							return false;

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
				if(bound!=0 && this->get_size() >= bound)
					return false;

				t2 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(3));
				add_sample(t2.lambda, t2.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

				t3 = sample(t1.lambda + scalar_type(2)*(t4.lambda - t1.lambda)/scalar_type(3));
				add_sample(t3.lambda, t3.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

				if(my_comparator.is_maybe_equal(t2.f_lambda, t4.f_lambda) || my_comparator.is_maybe_equal(t3.f_lambda, t4.f_lambda)){
					// meaning the function plot between t1 to t4 is a straight line
					t1 = t2 = t3 = t4;
					LOGGER(DEBUG7,"lower_bound_search:get_pivots[k,+infty]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
					return true;// true return means min is found already
				}
			}
		}
		else{ //we have finite bounds on the search space for minimum
			LOGGER(DEBUG7,"lower_bound_search:get_pivots[L,R]","Initiating pivots search");
			sample_type s;
			t1 = sample(my_interval.lower().get_val());
			add_sample(t1.lambda, t1.f_lambda);
			t4 = sample(my_interval.upper().get_val());
			add_sample(t4.lambda, t4.f_lambda);
			if(bound!=0 && this->get_size() >= bound)
				return false;

			sample_type s_prev;
			if(my_comparator.is_definitely_strictly_larger(t1.f_lambda, t4.f_lambda)){
				t2 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(2));
				add_sample(t2.lambda, t2.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

				if(is_colinear(t1,t2,t4)){
					t1 = t2 = t3 = t4 = t4;
					LOGGER(DEBUG7,"lower_bound_search:get_pivots[k1,k2]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
					return true;
				}
				// move towards t4 until sample is less than t4 sample.

				s = sample(t2.lambda + (t4.lambda - t2.lambda)/scalar_type(2));
				add_sample(s.lambda, s.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

				s_prev = t2;
				while(my_comparator.is_LE(t4.f_lambda, s.f_lambda)){ // This while loop will get stuck if the min is at or close to t4.
					s_prev = s;
					s = sample(s.lambda + (t4.lambda - s.lambda)/scalar_type(2));
					add_sample(s.lambda, s.f_lambda);
					if(bound!=0 && this->get_size() >= bound)
						return false;

					if(is_colinear(t4,s,s_prev)){
						t1 = t2 = t3 = t4 = t4;
						LOGGER(DEBUG7,"lower_bound_search:get_pivots[k1,k2]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
						return true;
					}

				}
				t3 = s;
			}
			else if(my_comparator.is_definitely_strictly_smaller(t1.f_lambda, t4.f_lambda)){
				t3 = sample(t1.lambda + (t4.lambda - t1.lambda)/scalar_type(2));
				add_sample(t3.lambda, t3.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

				if(is_colinear(t1,t3,t4)){
					t1 = t2 = t3 = t4 = t1;
					LOGGER(DEBUG7,"lower_bound_search:get_pivots[k1,k2]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
					return true;
				}
				// move left towards t1 until sample is less than t1 sample.
				s = sample(t1.lambda + (t3.lambda - t1.lambda)/scalar_type(2));
				add_sample(s.lambda, s.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

				s_prev = t3;
				while(math::maybe(my_comparator.is_LE(t1.f_lambda, s.f_lambda))){// This while loop will get stuck if the min is at or close to t1
					s_prev = s;
					s = sample(t1.lambda +(s.lambda - t1.lambda)/scalar_type(2));
					add_sample(s.lambda, s.f_lambda);
					if(bound!=0 && this->get_size() >= bound)
						return false;

					if(is_colinear(t1,s,s_prev)){
						LOGGER(DEBUG7,"lower_bound_search:get_pivots[k1,k2]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
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
				if(bound!=0 && this->get_size() >= bound)
					return false;

				t3 = sample(t1.lambda + scalar_type(2)*(t4.lambda - t1.lambda)/scalar_type(3));
				add_sample(t3.lambda, t3.f_lambda);
				if(bound!=0 && this->get_size() >= bound)
					return false;

				if(my_comparator.is_maybe_equal(t2.f_lambda,t4.f_lambda) || my_comparator.is_maybe_equal(t3.f_lambda, t4.f_lambda)){ // meaning the function plot between t1 to t4 is a straight line
					t1 = t2 = t3 = t4;
					LOGGER(DEBUG7,"lower_bound_search:get_pivots[k1,k2]","Pivots search terminated successfully with the minima found:"+to_string(t1.f_lambda));
					return true;// true return means min is found already

				}
			}
		}
		return false;
}

template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type, functor>::min_interval lower_bound_search<scalar_type,functor>::section_search
		(lower_bound_search<scalar_type, functor>::sample_type& s1,
		 lower_bound_search<scalar_type, functor>::sample_type& s2,
		 lower_bound_search<scalar_type, functor>::sample_type& s3,
		 lower_bound_search<scalar_type, functor>::sample_type& s4){

	 // A 4 pivot search algorithm, specialized to merge to exact minimum
	 // for pwa functions, on finite number of iterations.
	 //

	//Lets have an automata implementation here:

	unsigned int STATE = 0;
	scalar_type minimum, min_x;

	const scalar_type epsilon = this->get_interval_tolerance();

	scalar_with_infinity<scalar_type> lower_bound(scalar_with_infinity<scalar_type>::neg_infty());
	scalar_with_infinity<scalar_type> upper_bound(scalar_with_infinity<scalar_type>::pos_infty());
	bool minimum_found = false,upper_bound_initialised = false, lower_bound_initialised = false;
	typename lower_bound_search<scalar_type, functor>::sample_type new_sample, s5, intersect; // sample for the 5 point case

	math::numeric::approx_comparator<scalar_type> my_comparator;
	unsigned int dbg_cnt = 0;

	while(upper_bound.is_infinity() || lower_bound.is_infinity() ||
			my_comparator.is_definitely_strictly_larger((upper_bound.get_val() - lower_bound.get_val()),epsilon))
	{
		//if(upper_bound.is_finite() && lower_bound.is_finite())
		IFLOGGER(DEBUG7){
			print_pivots(s1,s2,s3,s4,s5);
		}

/*
		if(upper_bound_initialised && lower_bound_initialised)
		{

			std::cout << "low:" << lower_bound.get_val() << std::endl;
			std::cout << "up:" << upper_bound.get_val() << std::endl;

			std::cout << "up - low = " << upper_bound.get_val() - lower_bound.get_val() << std::endl;
			std::cout << "Current State:" << STATE << std::endl;

		}
*/
		if(STATE == 0 && my_comparator.is_definitely_strictly_larger(s2.f_lambda,s3.f_lambda)){ // case 1
			//some renaming
			//std::cout << "STATE 0: 1\n";
			STATE = 1;
		}
		else if(STATE == 0 && my_comparator.is_definitely_strictly_smaller(s2.f_lambda,s3.f_lambda)){ // case 2
			//std::cout << "STATE 0: 2\n";
			STATE = 2;

		}
		else if( STATE == 0 && my_comparator.is_maybe_equal(s2.f_lambda,s3.f_lambda)){
			//std::cout << "STATE 0: 3\n";
			STATE = 3;
		}
		else if(STATE == 1){
			//std::cout << "STATE 1" << std::endl;
			/* code for maintaining the bounds on min */

			if(!upper_bound_initialised){
				//upper_bound = scalar_with_infinity<scalar_type>(
				//			lower_bound_search<scalar_type, functor>::get_minimum(s2.f_lambda,s3.f_lambda));
				//upper_bound = scalar_with_infinity<scalar_type>(get_minimum(s2.f_lambda,s3.f_lambda));
				upper_bound = scalar_with_infinity<scalar_type>(s3.f_lambda);
				upper_bound_initialised = true;
			}
			else{
				//upper_bound = scalar_with_infinity<scalar_type>(lower_bound_search<scalar_type, functor>::get_minimum(upper_bound.get_val(),
				//												lower_bound_search<scalar_type, functor>::get_minimum(s2.f_lambda,s3.f_lambda)));
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(), s3.f_lambda));
			}

			typename lower_bound_search<scalar_type, functor>::interval my_interval;
			typename lower_bound_search<scalar_type, functor>::sample_type min_1,min_2,s;
			bool conv_pts = false;
			s.lambda = s4.lambda;
			s.f_lambda = scalar_type(0);
			min_1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			min_2 = this->template line_intersection_stable<precise_float>(s2,s3,s,s4,conv_pts);

			if(conv_pts){// heuristic to speed up computation and to bypass computation error
				// minima reached heuristics
				minimum_found = true;
				min_x = s3.lambda;
				minimum = s3.f_lambda;
				STATE = 6;
				continue;
			}
/*
			std::cout << "STATE 1: " << "min_1.x=" << min_1.lambda << ", min_1.y=" << min_1.f_lambda << std::endl;
			std::cout << "STATE 1: " << "min_2.x=" << min_2.lambda << ", min_2.y=" << min_2.f_lambda << std::endl;
*/

			if(!lower_bound_initialised){
				lower_bound = scalar_with_infinity<scalar_type>(get_minimum(min_1.f_lambda,min_2.f_lambda));
				lower_bound_initialised = true;
			}
			else{
				lower_bound = scalar_with_infinity<scalar_type>(get_maximum(lower_bound.get_val(),get_minimum(min_1.f_lambda,min_2.f_lambda)));
			}
			/*
				Choose the subinterval which has a lower minima
			*/
			bool left=false, right=false;
			if(my_comparator.is_definitely_strictly_smaller(min_1.f_lambda, min_2.f_lambda)){
				my_interval = typename lower_bound_search<scalar_type, functor>::interval(
						scalar_with_infinity<scalar_type>(s2.lambda),
						scalar_with_infinity<scalar_type>(s3.lambda));
				left = true;
			}
			else
			{
				my_interval = typename lower_bound_search<scalar_type, functor>::interval(
										scalar_with_infinity<scalar_type>(s3.lambda),
										scalar_with_infinity<scalar_type>(s4.lambda));
				right = true;
			}
			LOGGER(DEBUG7,"lb search:section_search:"," STATE 1 requesting to sample");
			new_sample = sample_interval(my_interval,left);
			add_sample(new_sample.lambda, new_sample.f_lambda);


			/* Renaming Pivots */
			/* left flag is true if the left subinterval (s2,s3) is the sampling interval.
			   right flag is true if the right subinterval (s3,s3) is the sampling interval.
			*/
			if (left) {
				if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda,s3.f_lambda)){
					// move to 5 point state
					s5 = s4;
					s4 = s3;
					s3 = new_sample;
					STATE = 4;
				}
				else if(my_comparator.is_maybe_equal(new_sample.f_lambda,s3.f_lambda)){
					// move to state 3;
					s1 = s2;
					s2 = new_sample;
					STATE = 3;
				}
				else{// new_sample.f_lambda > s3.f_lambda
					s1 = s2;
					s2 = new_sample;
					STATE = 1;
				}
			}
			else { // right subinterval (s3,s4) is chosen.
				if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda,s3.f_lambda)){
					s1 = s2;
					s2 = s3;
					s3 = new_sample;
					STATE = 1;
				}
				else if(my_comparator.is_maybe_equal(new_sample.f_lambda,s3.f_lambda)){
					//renaming
					s1 = s2;
					s2 = s3;
					s3 = new_sample;
					STATE = 3;
				}
				else{ // new_sample.f_lambda > s3.f_lambda
					//renaming
					s5 = s4;
					s4 = new_sample;
					STATE = 4; // The 5 point state
				}
			}
//			std::cout << "STATE 1 EXIT\n";
		}
		else if(STATE == 2){
			//std::cout << "STATE 2" << std::endl;
			typename lower_bound_search<scalar_type, functor>::interval my_interval;

			/* code for maintaining the bounds on min*/

			if(!upper_bound_initialised){
				//upper_bound = scalar_with_infinity<scalar_type>(lower_bound_search<scalar_type, functor>::get_minimum(s2.f_lambda,s3.f_lambda));
				upper_bound = scalar_with_infinity<scalar_type>(s2.f_lambda);
				upper_bound_initialised = true;
			}
			else{
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(), s2.f_lambda));
			}

			typename lower_bound_search<scalar_type, functor>::sample_type min_1,min_2,s;
			s.lambda = s1.lambda;
			s.f_lambda = scalar_type(0);
			bool conv_pts = false;
			min_1 = this->template line_intersection_stable<precise_float>(s2,s3,s,s1,conv_pts);
			min_2 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			if(conv_pts){// heuristic to speed up computation and to bypass computation error
				// minima reached heuristics
				minimum_found = true;
				min_x = s2.lambda;
				minimum = s2.f_lambda;
				STATE = 6;
				continue;
			}

			if(!lower_bound_initialised){
				lower_bound = scalar_with_infinity<scalar_type>(get_minimum(min_1.f_lambda,min_2.f_lambda));
				lower_bound_initialised = true;
			}
			else{
				lower_bound = scalar_with_infinity<scalar_type>(get_maximum(lower_bound.get_val(),get_minimum(min_1.f_lambda,min_2.f_lambda)));
			}

			/*
				Choose the subinterval which has a lower minima
			*/
			bool left=false, right=false;

			if(my_comparator.is_definitely_strictly_smaller(min_1.f_lambda, min_2.f_lambda)){
				my_interval = typename lower_bound_search<scalar_type, functor>::interval(
						scalar_with_infinity<scalar_type>(s1.lambda),
						scalar_with_infinity<scalar_type>(s2.lambda));
				left = true;
			}
			else
			{
				my_interval = typename lower_bound_search<scalar_type, functor>::interval(
										scalar_with_infinity<scalar_type>(s2.lambda),
										scalar_with_infinity<scalar_type>(s3.lambda));
				right = true;
			}

			//LOGGER(DEBUG7,"lb search:section_search:"," STATE 2 requesting to sample");
			new_sample = sample_interval(my_interval,left);
			add_sample(new_sample.lambda, new_sample.f_lambda);

			/* Renaming Pivots */
			if(left){
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
					//print_pivots(s1,s2,s3,s4,s5);

					STATE = 3;
	//				std::cout << "STATE 2: 3\n";
				}
			}
			else{ //right subinterval is chosen.
				if(my_comparator.is_definitely_strictly_larger(new_sample.f_lambda, s2.f_lambda)){
					s4 = s3;
					s3 = new_sample;
					STATE = 2;
				}
				else if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda , s2.f_lambda)){
					s5 = s4;
					s4 = s3;
					s3 = new_sample;
					STATE = 4;
				}
				else{ // new_sample.f_lambda == s2.f_lambda
					s4 = s3;
					s3 = new_sample;
					//print_pivots(s1,s2,s3,s4,s5);
					STATE = 3;
				}
			}
//			std::cout << "STATE 2 EXIT\n";
		}
		else if(STATE == 3){
			//std::cout << "STATE 3" << std::endl;
			/* code for maintaining the bounds on min*/

			if(!upper_bound_initialised){
				upper_bound = scalar_with_infinity<scalar_type>(s2.f_lambda);
				upper_bound_initialised = true;
			}
			else{
//				upper_bound = scalar_with_infinity<scalar_type>(
//				lower_bound_search<scalar_type, functor>::get_minimum(upper_bound.get_val(), s2.f_lambda));
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(), s2.f_lambda));
			}

			bool conv_pts = false;
			intersect = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			if(conv_pts){
				minimum_found = true;
				min_x = s3.lambda;
				minimum = s3.f_lambda;
				STATE = 6;
				continue;
			}

			LOGGER(DEBUG7,"lb search:section_search:"," STATE 3 requesting to sample");
			new_sample = sample(intersect.lambda);
			add_sample(new_sample.lambda, new_sample.f_lambda);


			if(!lower_bound_initialised){
				lower_bound = scalar_with_infinity<scalar_type>(intersect.f_lambda);
				lower_bound_initialised = true;
			}
			else{
				lower_bound = scalar_with_infinity<scalar_type>(get_maximum(lower_bound.get_val(),intersect.f_lambda));
			}

			/*----*/

			if(my_comparator.is_maybe_equal(new_sample.f_lambda,s2.f_lambda)){
				STATE = 6; //stop state
				minimum_found = true;
				min_x = s2.lambda;
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
				std::cout << "lower_bound_search: State 3: Something is wrong. Should not reach here logically\n";
				print_pivots(s1,s2,s3,s4,s5);
				std::cout << std::setprecision(22) << "new_sample.x:" << new_sample.lambda <<", new_sample.f_lambda:" << new_sample.f_lambda << std::endl;
				this->print_samples();
				this->plot_graph();
				throw std::runtime_error("STOPPING LB SEARCH DUE TO UNEXPECTED EXCEPTION");
			}
//			std::cout << "STATE 3 EXIT\n";
		}
		else if(STATE == 4){
			//std::cout << "STATE 4" << std::endl;

			typename lower_bound_search<scalar_type, functor>::sample_type intersect1, intersect2;

			/* code for maintaining the bounds on min*/

			if(!upper_bound_initialised){
				upper_bound = scalar_with_infinity<scalar_type>(s3.f_lambda);
				upper_bound_initialised = true;
			}
			else{
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(),s3.f_lambda));
			}

			bool left = false, right = false, conv_pts = false;

			intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5,conv_pts);

			if(conv_pts){// heuristic to speed up computation and to bypass computation error
				// minima reached heuristics
				minimum_found = true;
				min_x = s3.lambda;
				minimum = s3.f_lambda;
				STATE = 6;
				continue;
			}

			if(!lower_bound_initialised){
				lower_bound = scalar_with_infinity<scalar_type>(get_minimum(intersect1.f_lambda,intersect2.f_lambda));
				lower_bound_initialised = true;
			}
			else{
				lower_bound = scalar_with_infinity<scalar_type>(get_maximum(lower_bound.get_val(),get_minimum(intersect1.f_lambda,intersect2.f_lambda)));
			}

			/*----*/

			/*DEBUG*/
//			std::cout << std::setprecision(16) << "intersect1.x:" << intersect1.lambda <<", intersect1.y:" << intersect1.f_lambda << std::endl;
//			std::cout << std::setprecision(16) << "intersect2.x:" << intersect2.lambda <<", intersect2.y:" << intersect2.f_lambda << std::endl;

			/* Choose the sampling interval */
			/** Condition of finding the minima */
//			if(my_comparator.is_maybe_equal(intersect1.f_lambda, intersect2.f_lambda) &&
//					my_comparator.is_maybe_equal(intersect1.lambda, intersect2.lambda)){
			if(math::numeric::is_MEQ(intersect1.f_lambda-intersect2.f_lambda,scalar_type(0)) &&
			  math::numeric::is_MEQ(intersect1.lambda-intersect2.lambda,scalar_type(0))){

				minimum_found = true;
				minimum = intersect1.f_lambda;
				min_x = intersect1.lambda;
//				std::cout << "Reached inside STATE 4 equality condition\n";
				STATE = 6;
			}
			else{
				if (my_comparator.is_definitely_strictly_smaller(intersect1.f_lambda , intersect2.f_lambda)) {
					LOGGER(DEBUG7,"lb search:section_search:"," STATE 4 requesting to sample");
					new_sample = sample(intersect1.lambda);
					add_sample(new_sample.lambda, new_sample.f_lambda);

					left = true;
				}
				else{
					LOGGER(DEBUG7,"lb search:section_search:"," STATE 4 requesting to sample");
					new_sample = sample(intersect2.lambda);
					add_sample(new_sample.lambda, new_sample.f_lambda);

					right = true;
				}
				/*Renaming Pivots */

				if(left){
					if(my_comparator.is_definitely_strictly_larger(new_sample.f_lambda , s3.f_lambda)){
						s1 = s2;
						s2 = new_sample;
						STATE = 4;
					}
					else if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda , s3.f_lambda)){
						s5 = s4;
						s4 = s3;
						s3 = new_sample;
						STATE = 4;
					}
					else{ // new_sample.f_lambda == s3.f_lambda
						s1 = s2;
						s2 = new_sample;
						STATE = 3;
					}
				}
				else{// right
					if(my_comparator.is_definitely_strictly_larger(new_sample.f_lambda , s3.f_lambda)){
						s5 = s4;
						s4 = new_sample;
						STATE = 4;
					}
					else if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda , s3.f_lambda)){
						s1 = s2;
						s2 = s3;
						s3 = new_sample;
						STATE = 4;
					}
					else{//equal
						s1 = s2;
						s2 = s3;
						s3 = new_sample;
						STATE = 3;
					}
				}
			}
//			std::cout << "STATE 4 Exit\n";
		}

		else if(STATE == 6){ //stop state.
			//std::cout << "STATE 6\n";
			lower_bound = scalar_with_infinity<scalar_type>(minimum);
			upper_bound = scalar_with_infinity<scalar_type>(minimum);
//			std::cout << "STATE 6 EXIT\n";
		}
		else{}
	}//end of while
//	std::cout << "MIN PRB SOLVED" << std::endl;
	LOGGER(DEBUG1,"Lower_bound_search:section_search:","lb_search:#Samples:"+to_string(this->get_size()));

	if(minimum_found){
		typename math::numeric::interval<scalar_type> my_y_interval = math::numeric::interval<scalar_type>(scalar_with_infinity<scalar_type>(minimum),scalar_with_infinity<scalar_type>(minimum));
		typename math::numeric::interval<scalar_type> my_x_interval = math::numeric::interval<scalar_type>(
				scalar_with_infinity<scalar_type>(min_x),scalar_with_infinity<scalar_type>(min_x));

		typename convex_opt<scalar_type, functor>::min_interval min_structure;
		min_structure.lambda_interval = my_x_interval;
		min_structure.f_lambda_interval = my_y_interval;
		return min_structure;
	}
	else{
		typename lower_bound_search<scalar_type, functor>::interval my_x_interval(
						scalar_with_infinity<scalar_type>(s1.lambda),scalar_with_infinity<scalar_type>(s4.lambda));

		typename lower_bound_search<scalar_type, functor>::interval my_y_interval(lower_bound,upper_bound);


		typename convex_opt<scalar_type, functor>::min_interval min_structure = {my_x_interval, my_y_interval};
		return min_structure;
	}

}

/**
 * Search bounded on number of function samples
 * @param s1
 * @param s2
 * @param s3
 * @param s4
 * @return
 */
template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type, functor>::min_interval lower_bound_search<scalar_type,functor>::section_search_bounded
		(lower_bound_search<scalar_type, functor>::sample_type& s1,
		 lower_bound_search<scalar_type, functor>::sample_type& s2,
		 lower_bound_search<scalar_type, functor>::sample_type& s3,
		 lower_bound_search<scalar_type, functor>::sample_type& s4,
		 const unsigned int bound){

	 // A 4 pivot search algorithm, specialized to merge to exact minimum
	 // for pwa functions, on finite number of iterations.
	 //

	//Lets have an automata implementation here:

	unsigned int STATE = 0;
	scalar_type minimum, min_x;

	const scalar_type epsilon = this->get_interval_tolerance();

	scalar_with_infinity<scalar_type> lower_bound(scalar_with_infinity<scalar_type>::neg_infty());
	scalar_type min = get_minimum(s1.f_lambda,s2.f_lambda);
	min = get_minimum(min,s3.f_lambda);
	min = get_minimum(min, s4.f_lambda);
	scalar_with_infinity<scalar_type> upper_bound(min);
	bool minimum_found = false,upper_bound_initialised = true, lower_bound_initialised = false;
	typename lower_bound_search<scalar_type, functor>::sample_type new_sample, s5, intersect; // sample for the 5 point case

	math::numeric::approx_comparator<scalar_type> my_comparator;
	unsigned int dbg_cnt = 0;

	while(upper_bound.is_infinity() || lower_bound.is_infinity() ||
			my_comparator.is_definitely_strictly_larger((upper_bound.get_val() - lower_bound.get_val()),epsilon))
	{
		//check the bound here
		if(this->get_size() >= bound)
			break;
		//if(upper_bound.is_finite() && lower_bound.is_finite())
		IFLOGGER(DEBUG7){
			print_pivots(s1,s2,s3,s4,s5);
		}

/*
		if(upper_bound_initialised && lower_bound_initialised)
		{

			std::cout << "low:" << lower_bound.get_val() << std::endl;
			std::cout << "up:" << upper_bound.get_val() << std::endl;

			std::cout << "up - low = " << upper_bound.get_val() - lower_bound.get_val() << std::endl;
			std::cout << "Current State:" << STATE << std::endl;

		}
*/
		if(STATE == 0 && my_comparator.is_definitely_strictly_larger(s2.f_lambda,s3.f_lambda)){ // case 1
			//some renaming
			//std::cout << "STATE 0: 1\n";
			STATE = 1;
		}
		else if(STATE == 0 && my_comparator.is_definitely_strictly_smaller(s2.f_lambda,s3.f_lambda)){ // case 2
			//std::cout << "STATE 0: 2\n";
			STATE = 2;

		}
		else if( STATE == 0 && my_comparator.is_maybe_equal(s2.f_lambda,s3.f_lambda)){
			//std::cout << "STATE 0: 3\n";
			STATE = 3;
		}
		else if(STATE == 1){
			//std::cout << "STATE 1" << std::endl;
			/* code for maintaining the bounds on min */

			if(!upper_bound_initialised){
				//upper_bound = scalar_with_infinity<scalar_type>(
				//			lower_bound_search<scalar_type, functor>::get_minimum(s2.f_lambda,s3.f_lambda));
				//upper_bound = scalar_with_infinity<scalar_type>(get_minimum(s2.f_lambda,s3.f_lambda));
				upper_bound = scalar_with_infinity<scalar_type>(s3.f_lambda);
				upper_bound_initialised = true;
			}
			else{
				//upper_bound = scalar_with_infinity<scalar_type>(lower_bound_search<scalar_type, functor>::get_minimum(upper_bound.get_val(),
				//												lower_bound_search<scalar_type, functor>::get_minimum(s2.f_lambda,s3.f_lambda)));
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(), s3.f_lambda));
			}

			typename lower_bound_search<scalar_type, functor>::interval my_interval;
			typename lower_bound_search<scalar_type, functor>::sample_type min_1,min_2,s;
			bool conv_pts = false;
			s.lambda = s4.lambda;
			s.f_lambda = scalar_type(0);
			min_1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			min_2 = this->template line_intersection_stable<precise_float>(s2,s3,s,s4,conv_pts);

			if(conv_pts){// heuristic to speed up computation and to bypass computation error
				// minima reached heuristics
				minimum_found = true;
				min_x = s3.lambda;
				minimum = s3.f_lambda;
				STATE = 6;
				continue;
			}
/*
			std::cout << "STATE 1: " << "min_1.x=" << min_1.lambda << ", min_1.y=" << min_1.f_lambda << std::endl;
			std::cout << "STATE 1: " << "min_2.x=" << min_2.lambda << ", min_2.y=" << min_2.f_lambda << std::endl;
*/

			if(!lower_bound_initialised){
				lower_bound = scalar_with_infinity<scalar_type>(get_minimum(min_1.f_lambda,min_2.f_lambda));
				lower_bound_initialised = true;
			}
			else{
				lower_bound = scalar_with_infinity<scalar_type>(get_maximum(lower_bound.get_val(),get_minimum(min_1.f_lambda,min_2.f_lambda)));
			}
			/*
				Choose the subinterval which has a lower minima
			*/
			bool left=false, right=false;
			if(my_comparator.is_definitely_strictly_smaller(min_1.f_lambda, min_2.f_lambda)){
				my_interval = typename lower_bound_search<scalar_type, functor>::interval(
						scalar_with_infinity<scalar_type>(s2.lambda),
						scalar_with_infinity<scalar_type>(s3.lambda));
				left = true;
			}
			else
			{
				my_interval = typename lower_bound_search<scalar_type, functor>::interval(
										scalar_with_infinity<scalar_type>(s3.lambda),
										scalar_with_infinity<scalar_type>(s4.lambda));
				right = true;
			}
			LOGGER(DEBUG7,"lb search:section_search:"," STATE 1 requesting to sample");
			new_sample = sample_interval(my_interval,left);
			add_sample(new_sample.lambda, new_sample.f_lambda);


			/* Renaming Pivots */
			/* left flag is true if the left subinterval (s2,s3) is the sampling interval.
			   right flag is true if the right subinterval (s3,s3) is the sampling interval.
			*/
			if (left) {
				if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda,s3.f_lambda)){
					// move to 5 point state
					s5 = s4;
					s4 = s3;
					s3 = new_sample;
					STATE = 4;
				}
				else if(my_comparator.is_maybe_equal(new_sample.f_lambda,s3.f_lambda)){
					// move to state 3;
					s1 = s2;
					s2 = new_sample;
					STATE = 3;
				}
				else{// new_sample.f_lambda > s3.f_lambda
					s1 = s2;
					s2 = new_sample;
					STATE = 1;
				}
			}
			else { // right subinterval (s3,s4) is chosen.
				if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda,s3.f_lambda)){
					s1 = s2;
					s2 = s3;
					s3 = new_sample;
					STATE = 1;
				}
				else if(my_comparator.is_maybe_equal(new_sample.f_lambda,s3.f_lambda)){
					//renaming
					s1 = s2;
					s2 = s3;
					s3 = new_sample;
					STATE = 3;
				}
				else{ // new_sample.f_lambda > s3.f_lambda
					//renaming
					s5 = s4;
					s4 = new_sample;
					STATE = 4; // The 5 point state
				}
			}
//			std::cout << "STATE 1 EXIT\n";
		}
		else if(STATE == 2){
			//std::cout << "STATE 2" << std::endl;
			typename lower_bound_search<scalar_type, functor>::interval my_interval;

			/* code for maintaining the bounds on min*/

			if(!upper_bound_initialised){
				//upper_bound = scalar_with_infinity<scalar_type>(lower_bound_search<scalar_type, functor>::get_minimum(s2.f_lambda,s3.f_lambda));
				upper_bound = scalar_with_infinity<scalar_type>(s2.f_lambda);
				upper_bound_initialised = true;
			}
			else{
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(), s2.f_lambda));
			}

			typename lower_bound_search<scalar_type, functor>::sample_type min_1,min_2,s;
			s.lambda = s1.lambda;
			s.f_lambda = scalar_type(0);
			bool conv_pts = false;
			min_1 = this->template line_intersection_stable<precise_float>(s2,s3,s,s1,conv_pts);
			min_2 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			if(conv_pts){// heuristic to speed up computation and to bypass computation error
				// minima reached heuristics
				minimum_found = true;
				min_x = s2.lambda;
				minimum = s2.f_lambda;
				STATE = 6;
				continue;
			}

			if(!lower_bound_initialised){
				lower_bound = scalar_with_infinity<scalar_type>(get_minimum(min_1.f_lambda,min_2.f_lambda));
				lower_bound_initialised = true;
			}
			else{
				lower_bound = scalar_with_infinity<scalar_type>(get_maximum(lower_bound.get_val(),get_minimum(min_1.f_lambda,min_2.f_lambda)));
			}

			/*
				Choose the subinterval which has a lower minima
			*/
			bool left=false, right=false;

			if(my_comparator.is_definitely_strictly_smaller(min_1.f_lambda, min_2.f_lambda)){
				my_interval = typename lower_bound_search<scalar_type, functor>::interval(
						scalar_with_infinity<scalar_type>(s1.lambda),
						scalar_with_infinity<scalar_type>(s2.lambda));
				left = true;
			}
			else
			{
				my_interval = typename lower_bound_search<scalar_type, functor>::interval(
										scalar_with_infinity<scalar_type>(s2.lambda),
										scalar_with_infinity<scalar_type>(s3.lambda));
				right = true;
			}

			//LOGGER(DEBUG7,"lb search:section_search:"," STATE 2 requesting to sample");
			new_sample = sample_interval(my_interval,left);
			add_sample(new_sample.lambda, new_sample.f_lambda);

			/* Renaming Pivots */
			if(left){
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
					//print_pivots(s1,s2,s3,s4,s5);

					STATE = 3;
	//				std::cout << "STATE 2: 3\n";
				}
			}
			else{ //right subinterval is chosen.
				if(my_comparator.is_definitely_strictly_larger(new_sample.f_lambda, s2.f_lambda)){
					s4 = s3;
					s3 = new_sample;
					STATE = 2;
				}
				else if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda , s2.f_lambda)){
					s5 = s4;
					s4 = s3;
					s3 = new_sample;
					STATE = 4;
				}
				else{ // new_sample.f_lambda == s2.f_lambda
					s4 = s3;
					s3 = new_sample;
					//print_pivots(s1,s2,s3,s4,s5);
					STATE = 3;
				}
			}
//			std::cout << "STATE 2 EXIT\n";
		}
		else if(STATE == 3){
			//std::cout << "STATE 3" << std::endl;
			/* code for maintaining the bounds on min*/

			if(!upper_bound_initialised){
				upper_bound = scalar_with_infinity<scalar_type>(s2.f_lambda);
				upper_bound_initialised = true;
			}
			else{
//				upper_bound = scalar_with_infinity<scalar_type>(
//				lower_bound_search<scalar_type, functor>::get_minimum(upper_bound.get_val(), s2.f_lambda));
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(), s2.f_lambda));
			}

			bool conv_pts = false;
			intersect = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			if(conv_pts){
				minimum_found = true;
				min_x = s3.lambda;
				minimum = s3.f_lambda;
				STATE = 6;
				continue;
			}

			LOGGER(DEBUG7,"lb search:section_search:"," STATE 3 requesting to sample");
			new_sample = sample(intersect.lambda);
			add_sample(new_sample.lambda, new_sample.f_lambda);


			if(!lower_bound_initialised){
				lower_bound = scalar_with_infinity<scalar_type>(intersect.f_lambda);
				lower_bound_initialised = true;
			}
			else{
				lower_bound = scalar_with_infinity<scalar_type>(get_maximum(lower_bound.get_val(),intersect.f_lambda));
			}

			/*----*/

			if(my_comparator.is_maybe_equal(new_sample.f_lambda,s2.f_lambda)){
				STATE = 6; //stop state
				minimum_found = true;
				min_x = s2.lambda;
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
				std::cout << "lower_bound_search: State 3: Something is wrong. Should not reach here logically\n";
				print_pivots(s1,s2,s3,s4,s5);
				std::cout << std::setprecision(22) << "new_sample.x:" << new_sample.lambda <<", new_sample.f_lambda:" << new_sample.f_lambda << std::endl;
				this->print_samples();
				this->plot_graph();
				throw std::runtime_error("STOPPING LB SEARCH DUE TO UNEXPECTED EXCEPTION");
			}
//			std::cout << "STATE 3 EXIT\n";
		}
		else if(STATE == 4){
			//std::cout << "STATE 4" << std::endl;

			typename lower_bound_search<scalar_type, functor>::sample_type intersect1, intersect2;

			/* code for maintaining the bounds on min*/

			if(!upper_bound_initialised){
				upper_bound = scalar_with_infinity<scalar_type>(s3.f_lambda);
				upper_bound_initialised = true;
			}
			else{
				upper_bound = scalar_with_infinity<scalar_type>(get_minimum(upper_bound.get_val(),s3.f_lambda));
			}

			bool left = false, right = false, conv_pts = false;

			intersect1 = this->template line_intersection_stable<precise_float>(s1,s2,s3,s4,conv_pts);
			intersect2 = this->template line_intersection_stable<precise_float>(s2,s3,s4,s5,conv_pts);

			if(conv_pts){// heuristic to speed up computation and to bypass computation error
				// minima reached heuristics
				minimum_found = true;
				min_x = s3.lambda;
				minimum = s3.f_lambda;
				STATE = 6;
				continue;
			}

			if(!lower_bound_initialised){
				lower_bound = scalar_with_infinity<scalar_type>(get_minimum(intersect1.f_lambda,intersect2.f_lambda));
				lower_bound_initialised = true;
			}
			else{
				lower_bound = scalar_with_infinity<scalar_type>(get_maximum(lower_bound.get_val(),get_minimum(intersect1.f_lambda,intersect2.f_lambda)));
			}

			/*----*/

			/*DEBUG*/
//			std::cout << std::setprecision(16) << "intersect1.x:" << intersect1.lambda <<", intersect1.y:" << intersect1.f_lambda << std::endl;
//			std::cout << std::setprecision(16) << "intersect2.x:" << intersect2.lambda <<", intersect2.y:" << intersect2.f_lambda << std::endl;

			/* Choose the sampling interval */
			/** Condition of finding the minima */
//			if(my_comparator.is_maybe_equal(intersect1.f_lambda, intersect2.f_lambda) &&
//					my_comparator.is_maybe_equal(intersect1.lambda, intersect2.lambda)){
			if(math::numeric::is_MEQ(intersect1.f_lambda-intersect2.f_lambda,scalar_type(0)) &&
			  math::numeric::is_MEQ(intersect1.lambda-intersect2.lambda,scalar_type(0))){

				minimum_found = true;
				minimum = intersect1.f_lambda;
				min_x = intersect1.lambda;
//				std::cout << "Reached inside STATE 4 equality condition\n";
				STATE = 6;
			}
			else{
				if (my_comparator.is_definitely_strictly_smaller(intersect1.f_lambda , intersect2.f_lambda)) {
					LOGGER(DEBUG7,"lb search:section_search:"," STATE 4 requesting to sample");
					new_sample = sample(intersect1.lambda);
					add_sample(new_sample.lambda, new_sample.f_lambda);

					left = true;
				}
				else{
					LOGGER(DEBUG7,"lb search:section_search:"," STATE 4 requesting to sample");
					new_sample = sample(intersect2.lambda);
					add_sample(new_sample.lambda, new_sample.f_lambda);

					right = true;
				}
				/*Renaming Pivots */

				if(left){
					if(my_comparator.is_definitely_strictly_larger(new_sample.f_lambda , s3.f_lambda)){
						s1 = s2;
						s2 = new_sample;
						STATE = 4;
					}
					else if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda , s3.f_lambda)){
						s5 = s4;
						s4 = s3;
						s3 = new_sample;
						STATE = 4;
					}
					else{ // new_sample.f_lambda == s3.f_lambda
						s1 = s2;
						s2 = new_sample;
						STATE = 3;
					}
				}
				else{// right
					if(my_comparator.is_definitely_strictly_larger(new_sample.f_lambda , s3.f_lambda)){
						s5 = s4;
						s4 = new_sample;
						STATE = 4;
					}
					else if(my_comparator.is_definitely_strictly_smaller(new_sample.f_lambda , s3.f_lambda)){
						s1 = s2;
						s2 = s3;
						s3 = new_sample;
						STATE = 4;
					}
					else{//equal
						s1 = s2;
						s2 = s3;
						s3 = new_sample;
						STATE = 3;
					}
				}
			}
//			std::cout << "STATE 4 Exit\n";
		}

		else if(STATE == 6){ //stop state.
			//std::cout << "STATE 6\n";
			lower_bound = scalar_with_infinity<scalar_type>(minimum);
			upper_bound = scalar_with_infinity<scalar_type>(minimum);
//			std::cout << "STATE 6 EXIT\n";
		}
		else{}

	}//end of while
//	std::cout << "MIN PRB SOLVED" << std::endl;
	LOGGER(DEBUG1,"Lower_bound_search:section_search:","lb_search:#Samples:"+to_string(this->get_size()));

	if(minimum_found){
		typename math::numeric::interval<scalar_type> my_y_interval = math::numeric::interval<scalar_type>(scalar_with_infinity<scalar_type>(minimum),scalar_with_infinity<scalar_type>(minimum));
		typename math::numeric::interval<scalar_type> my_x_interval = math::numeric::interval<scalar_type>(
				scalar_with_infinity<scalar_type>(min_x),scalar_with_infinity<scalar_type>(min_x));

		typename convex_opt<scalar_type, functor>::min_interval min_structure;
		min_structure.lambda_interval = my_x_interval;
		min_structure.f_lambda_interval = my_y_interval;
		return min_structure;
	}
	else{
		typename lower_bound_search<scalar_type, functor>::interval my_x_interval(
						scalar_with_infinity<scalar_type>(s1.lambda),scalar_with_infinity<scalar_type>(s4.lambda));

		typename lower_bound_search<scalar_type, functor>::interval my_y_interval(lower_bound,upper_bound);


		typename convex_opt<scalar_type, functor>::min_interval min_structure = {my_x_interval, my_y_interval};
		return min_structure;
	}

}
template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type, functor>::min_interval lower_bound_search<scalar_type,functor>::section_search(
		const lower_bound_search<scalar_type,functor>::interval& my_interval) {

	LOGGERSW(DEBUG1,"lower_bound_search: section_search:","lower_bound_search executing with error bound:"+to_string(this->get_interval_tolerance()));

	sample_type t1,t2,t3,t4;
	bool min_found_flag = false;
	typename lower_bound_search<scalar_type,functor>::interval min_bounds;

	if(this->get_minbrak_type() == "gold_desc"){
		if(!my_interval.upper().is_finite() && !my_interval.lower().is_finite()){
			minbrak_simple(t1,t2,t3,t4);
			min_found_flag = is_min_found(t1,t2,t3,t4);
		}
		else
			min_found_flag = get_pivots(my_interval, t1,t2,t3,t4);
	}
	else if(this->get_minbrak_type() == "parab_desc"){
		minbrak_para_ext(t1,t2,t3,t4);
	}
	else{
		throw std::runtime_error("lower_bound_search: minima bracketing method not known");
	}

	typename convex_opt<scalar_type, functor>::min_interval min_structure;

	LOGGER(DEBUG7,"lower_bound_search:section_search:","End of Pivots selection");

	if(min_found_flag){ // min is already found;
//		std::cout << "min found after pivot selection\n";
		typename lower_bound_search<scalar_type,functor>::interval my_x_interval(scalar_with_infinity<scalar_type>(t1.lambda),
																 	            scalar_with_infinity<scalar_type>(t2.lambda));
		typename lower_bound_search<scalar_type,functor>::interval my_y_interval(scalar_with_infinity<scalar_type>(t1.f_lambda),
																		 	            scalar_with_infinity<scalar_type>(t2.f_lambda));

		min_structure.lambda_interval = my_x_interval;
		min_structure.f_lambda_interval = my_y_interval;

		LOGGER(DEBUG1,"Lower_bound_search:section_search:","min found at pivot search with #Samples:"+to_string(this->get_size()));

#ifdef PLOT_INTERSECTION_SAMPLES
			this->plot_graph();
#endif

		return  min_structure;
	}

//	print_pivots(t1,t2,t3,t4,t4);

	min_structure = section_search(t1,t2,t3,t4);
#ifdef PLOT_INTERSECTION_SAMPLES
		this->plot_graph();
#endif

	return min_structure;
}
/**
 * Search bounded on the number of function samples.
 *
 * @param my_interval
 * @param bound
 * @return
 */
template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type, functor>::min_interval lower_bound_search<scalar_type,functor>::section_search_bounded(
		const lower_bound_search<scalar_type,functor>::interval& my_interval, const unsigned int bound) {

	LOGGERSW(DEBUG1,"lower_bound_search: section_search:","lower_bound_search executing with error bound:"+to_string(this->get_interval_tolerance()));

	sample_type t1,t2,t3,t4;
	t1.f_lambda = t2.f_lambda = t3.f_lambda = t4.f_lambda = INT_MAX;
	bool min_found_flag = false;

	if(this->get_minbrak_type() == "gold_desc"){
		if(!my_interval.upper().is_finite() && !my_interval.lower().is_finite()){
			minbrak_simple(t1,t2,t3,t4, bound);
			min_found_flag = is_min_found(t1,t2,t3,t4);
		}
		else
			min_found_flag = get_pivots(my_interval, t1,t2,t3,t4, bound);
	}
	else if(this->get_minbrak_type() == "parab_desc"){
		minbrak_para_ext(t1,t2,t3,t4, bound);
	}
	else{
		throw std::runtime_error("lower_bound_search: minima bracketing method not known");
	}

	typename convex_opt<scalar_type, functor>::min_interval min_structure;

	LOGGER(DEBUG7,"lower_bound_search:section_search:","End of Pivots selection");

	if(min_found_flag){ // min is already found;
//		std::cout << "min found after pivot selection\n";
		typename lower_bound_search<scalar_type,functor>::interval my_x_interval(scalar_with_infinity<scalar_type>(t1.lambda),
																 	            scalar_with_infinity<scalar_type>(t2.lambda));
		typename lower_bound_search<scalar_type,functor>::interval my_y_interval(scalar_with_infinity<scalar_type>(t1.f_lambda),
																		 	            scalar_with_infinity<scalar_type>(t2.f_lambda));

		min_structure.lambda_interval = my_x_interval;
		min_structure.f_lambda_interval = my_y_interval;

		LOGGER(DEBUG1,"Lower_bound_search:section_search:","min found at pivot search with #Samples:"+to_string(this->get_size()));

#ifdef PLOT_INTERSECTION_SAMPLES
			this->plot_graph();
#endif

		return  min_structure;
	}

//	print_pivots(t1,t2,t3,t4,t4);

	min_structure = section_search_bounded(t1,t2,t3,t4,bound);
#ifdef PLOT_INTERSECTION_SAMPLES
		this->plot_graph();
#endif

	return min_structure;
}

template<class scalar_type, template<typename > class functor>
typename convex_opt<scalar_type, functor>::min_interval lower_bound_search<scalar_type,functor>::section_search(){

	typename lower_bound_search<scalar_type,functor>::interval my_interval(scalar_with_infinity<scalar_type>::neg_infty(),
																  scalar_with_infinity<scalar_type>::pos_infty());
	return section_search(my_interval);
}
template<class scalar_type, template<typename > class functor>
bool lower_bound_search<scalar_type,functor>::is_min_found(lower_bound_search<scalar_type, functor>::sample_type& s1,
		lower_bound_search<scalar_type, functor>::sample_type& s2,
		lower_bound_search<scalar_type, functor>::sample_type& s3,
		lower_bound_search<scalar_type, functor>::sample_type& s4){
	math::numeric::approx_comparator<scalar_type> comp;
	//check s1,s2,s3
	if(comp.is_equal(s1.f_lambda,s2.f_lambda) && comp.is_equal(s2.f_lambda,s3.f_lambda)){
		s1=s2=s3=s4=s1;
		return true;
	}
	//check s2,s3,s4
	if(comp.is_equal(s2.f_lambda,s3.f_lambda) && comp.is_equal(s3.f_lambda,s4.f_lambda)){
		s1=s2=s3=s4=s2;
		return true;
	}
	//check s1,s3,s4
	if(comp.is_equal(s1.f_lambda,s3.f_lambda) && comp.is_equal(s3.f_lambda,s4.f_lambda)){
		s1=s2=s3=s4=s1;
		return true;
	}
	//check s1,s2,s4
	if(comp.is_equal(s1.f_lambda,s2.f_lambda) && comp.is_equal(s2.f_lambda,s4.f_lambda)){
		s1=s2=s3=s4=s1;
		return true;
	}
	return false;
}

} // end of namespace support_function
} // end of namespace continuous
