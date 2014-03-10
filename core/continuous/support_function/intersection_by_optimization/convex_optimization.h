#ifndef _CONVEX_OPT
#define _CONVEX_OPT


#include "core/continuous/support_function_provider.h"
#include "core/continuous/support_function/sf_derived/sf_sum.h"
#include "core/continuous/support_function/sf_base/sf_unary.h"
#include "math/numeric/approx_comparator.h"
#include "math/type_conversion.h"
#include "math/scalar_types/scalar_with_infinity.h"
#include "math/matrix.h"
#include "math/tribool.h"
#include "sample_plotter.h"
#include "math/numeric/interval.h"
#include "utility/logger_stopwatch.h"

namespace continuous {
namespace support_function {

/**
 * Different functions
 *
 */
template<typename scalar_type>
class sup_functor{
public:
	typedef math::vdom_vector<scalar_type> direction;

	/**
	 * Semantics of the members in the context of minimizing the
	 * intersection of a support function providing set with an hyperplane
	 * are as follows:
	 *
	 * @param l the direction in which to compute the support function
	 * @param c the normal to the hyperplane
	 * @param d the inh_term of the hyperplane
	 * @param S the Support function provider set
	 *
	 */
	sup_functor(const direction l,const  direction c, const scalar_type d,
			const continuous::support_function_provider::const_ptr S) : my_v(l), my_c(c), my_d(d), my_S(S){};

	scalar_type compute(const scalar_type& arg) const {
		bool is_empty, is_bounded;
//		std::cout << "my_v" << my_v << ", my_c:" << my_c <<std::endl;
		//math::vdom_vector<scalar_type> le = (my_v - arg * my_c);
		math::vdom_vector<scalar_type> le(my_c);
		le*=(-arg);
		le+=my_v;

		LOGGER_OS(DEBUG7,"sup_functor:compute") << "direction:"  << my_v << " with guard normal:" << my_c;

		math::vdom_vector<scalar_type> sv;
		scalar_type sopt;
		my_S->compute_support(le, sopt, sv, is_empty, is_bounded);
		if (!is_bounded)
			throw std::runtime_error("convex_opt_sample: unbounded set\n");
		if (is_empty)
			throw std::runtime_error("convex_opt_sample: empty set\n");

		LOGGER_OS(DEBUG7,"sup_functor:compute") << "l - lambda. n:" << le << " with sup = " << sopt;

		return sopt;
	};
	direction get_direction() const {
		return my_v;
	}
	direction get_cons_normal() const {
		return my_c;
	}
	continuous::support_function_provider::const_ptr get_S() const {
		return my_S;
	}

private:
	direction my_v; // The direction in which to compute the intersection
	direction my_c; // Normal to the hyperplane with which to compute the intersection
	scalar_type my_d;
					  /* Semantically, it is the inhomogenous term of the equation of a hyperplane. The Hyperplane
					   * under consideration is the one with which we are interested in computing the intersection
					   * precisely.
					   */
	continuous::support_function_provider::const_ptr my_S;
};

template<typename scalar_type>
class theta_functor{
public:
	typedef math::vdom_vector<scalar_type> direction;

	theta_functor(const direction l,const  direction c, const scalar_type d,
			const continuous::support_function_provider::const_ptr S) : my_v(l), my_c(c), my_d(d), my_S(S){};

	scalar_type compute(const scalar_type& theta) const {
		bool is_empty, is_bounded;
		math::vdom_vector<scalar_type> sv;

		unsigned int n = my_v.size();
		math::matrix<scalar_type> pi_transpose(n,2); // pi is the n-dim space to 2-dim space linear transformation matrix.

		// Initialize pi_transpose
		for(unsigned int i=0;i<n;i++){
			pi_transpose(i,0) = my_c[i];
			pi_transpose(i,1) = my_v[i];
		}

		math::vector<scalar_type> direction(2);
		direction[0] = std::cos(theta);
		direction[1] = std::sin(theta);

		// to be continued from here
		math::vector<scalar_type> trans_direction = pi_transpose.multiply_vector(direction);
		math::vdom_vector<scalar_type> my_vdom_v(trans_direction, my_v.get_index_to_variable_id_map());

		scalar_type sopt;
		my_S->compute_support(my_vdom_v, sopt, sv, is_empty, is_bounded);

		return (sopt - my_d * std::cos(theta) )/std::sin(theta);
	};


private:
	direction my_v; // The direction in which to compute the intersection
	direction my_c; // Normal to the hyperplane with which to compute the intersection
	scalar_type my_d;
					  /* Semantically, it is the inhomogenous term of the equation of a hyperplane. The Hyperplane
					   * under consideration is the one with which we are interested in computing the intersection
					   * precisely.
					   */
	continuous::support_function_provider::const_ptr my_S;

};

/**
 * Another sup functor with arguments as theta.
 */
template<class scalar_type>
class support_theta{
public:
	typedef math::vdom_vector<scalar_type> direction;

	support_theta(const continuous::support_function_provider::const_ptr S) : my_S(S){};

	scalar_type compute(const scalar_type& theta) const {
		bool is_empty, is_bounded;

		math::vector<scalar_type> direction(2);
		direction[0] = std::cos(theta);
		direction[1] = std::sin(theta);

		variable_id_set vis;
		index_to_variable_id_map_ptr iimap =
				index_to_variable_id_map::get_map_with_ids(my_S->get_variable_ids());


		math::vdom_vector<scalar_type> my_vdom_v(direction,iimap);
		math::vdom_vector<scalar_type> sv;

		scalar_type sopt;
		my_S->compute_support(my_vdom_v, sopt, sv, is_empty, is_bounded);
		if (!is_bounded)
			throw std::runtime_error("convex_opt_sample: unbounded set\n");
		if (is_empty)
			throw std::runtime_error("convex_opt_sample: empty set\n");

		return sopt;
	};
private:
	continuous::support_function_provider::const_ptr my_S;
};

/**
 * This class to provide interfaces required to compute precise intersection operation between
 * a convex set and a hyperplane c.x=d.
 *
 * @author Rajarshi
 */


template<class scalar_type, template<typename > class functor>
class convex_opt : public sample_plotter<scalar_type>{
public:
	/**
	 * This structure is a function and argument pair, where the function is from Reals to Reals.
	 */
	struct sample_type {
			scalar_type lambda;
			scalar_type f_lambda;
	};

	struct min_interval{
		math::numeric::interval<scalar_type> lambda_interval;
		math::numeric::interval<scalar_type> f_lambda_interval;
	};

	typedef boost::shared_ptr<convex_opt<scalar_type, functor> > ptr;
	typedef boost::shared_ptr<const convex_opt<scalar_type, functor> > const_ptr;
//	typedef typename convex_opt<scalar_type, functor>::sample sample_type;
	typedef typename convex_opt<scalar_type, functor>::min_interval min_interval;
//	typedef typename boost::numeric::interval<scalar_with_infinity<scalar_type> > interval;
	typedef typename math::numeric::interval<scalar_type> interval;

	/* Constructors */

	convex_opt(functor<scalar_type>& f, const std::string minbrak_type,
			const math::affine_map<scalar_type> dynamics_map,
			const support_function_provider::const_ptr U_ptr,
			const scalar_type eps) :
		my_functor(f), my_minbrak_type(minbrak_type), my_dynamics_map(dynamics_map), my_U(U_ptr), my_eps(eps){};

	convex_opt(functor<scalar_type>& f, const std::string minbrak_type,
			const scalar_type eps) : my_functor(f), my_minbrak_type(minbrak_type), my_eps(eps){
		my_dynamics_map = math::affine_map<scalar_type>(); //empty map created.

	};

	convex_opt(functor<scalar_type>& f, const scalar_type eps) : my_functor(f), my_eps(eps){
		my_dynamics_map = math::affine_map<scalar_type>(); // empty map created.
		my_minbrak_type = "gold_desc";
	};

	virtual ~convex_opt(){};
	/*
	 * Samples the objective function at lopt
	 * @param lopt
	 */
	 sample_type sample(const scalar_type& lopt) const;
	/**
	 * Given a tolerance epsilon, it returns the interval containing
	 * the min of a unimodal function such that the
	 * size of the interval <= epsilon.
	 *
	 * @return The structure with the minimum containing and the minimum(exact or approx)
	 * of a convex function
	 */
	virtual min_interval section_search() = 0;
	/**
	 * Given an initial search interval and a tolerance epsilon, it returns the
	 * interval containing the min of a unimodal function such that the
	 * size of the interval <= epsilon.
	 *
	 * @return The structure with the minimum containing and the minimum(exact or approx)
	 * of a convex function	 */
	virtual min_interval section_search(const interval& search_interval)  = 0;

	/**
	 * Returns the epsilon tolerance of the interval.
	 */
	scalar_type get_interval_tolerance() const
	{return my_eps;	};

	/**
	 * Return the minima bracketing method identifier string
	 * @return
	 */
	std::string get_minbrak_type() const {
		return my_minbrak_type;
	}

	/**
	 * Return the dynamics map associated. If not initialized, an empty
	 * map is returned.
	 *
	 * @return
	 */
	math::affine_map<scalar_type> get_dynamics_map() const {
		return my_dynamics_map;
	}

	/*
	 * Computes the intersection of 2 straight lines.
	 * @return The point of intersection. Throws runtime exception if the lines are parallel
	 */
	static sample_type line_intersection(sample_type s1, sample_type s2, sample_type p1, sample_type p2);
	/*
	 * More numerically stable computation of the intersection of 2 straight lines.
	 * @return The point of intersection. Throws runtime exception if the lines are parallel
	 */
	template<typename float_type>
	static sample_type line_intersection_stable(
			sample_type s1, sample_type s2,
			sample_type p1, sample_type p2,
			bool& convergent = false);

	/**
	 * Computes the equation of the line passing through s1 and s2 and samples
	 * the line at lambda.
	 * @param s1 a 2 dimensional point
	 * @param s2 a 2 dimensional point
	 * @param lambda
	 * @return
	 */
	static sample_type sample_line(sample_type s1, sample_type s2, scalar_type lambda);
	/**
	 * Samples the objective function inside the interval given by interval my_interval
	 *
	 * @param my_interval The interval within which to sample.
	 * @param side True meaning the interval is on the left of the middle sample
	 * of the 3 min-bracketing points. False, meaning the interval is on the right
	 * of the middle sample of the 3 min-bracketing points.
	 * @return sample point inside the interval
	 */
	sample_type sample_interval(const interval& my_interval, bool side) const;
	/**
	 * given an interval, it returns a sampling request inside the interval. Current implementation
	 * chooses the mid point of the passed interval as the sample request.
	 *
	 * @param my_interval The interval within which to request the sampling point
	 * @param side True meaning the interval is on the left of the middle sample
	 * of the 3 min-bracketing points. False, meaning the interval is on the right
	 * of the middle sample of the 3 min-bracketing points.
	 *
	 * @return The sampling point requested
	 */
	scalar_type sample_request(const interval& my_interval, bool side) const;

	static bool is_colinear(sample_type s1, sample_type s2, sample_type s3);

	/**
	 * Right shifts the samples a->b, b->c and c->d.
	 *
	 * @param a
	 * @param b
	 * @param c
	 * @param d
	 */
	static inline void shift3( sample_type& a, sample_type& b,
			sample_type& c,
			sample_type& d ){

		a.lambda = b.lambda;
		b.lambda = c.lambda;
		c.lambda = d.lambda;

		a.f_lambda = b.f_lambda;
		b.f_lambda = c.f_lambda;
		c.f_lambda = d.f_lambda;

	};
	/**
	 * Right shifts the samples a->b and b->c.
	 *
	 * @param a
	 * @param b
	 * @param c
	 */
	static inline void shift2(sample_type& a, sample_type& b,
			sample_type& c){
		a.lambda = b.lambda;
		b.lambda = c.lambda;

		a.f_lambda = b.f_lambda;
		b.f_lambda = c.f_lambda;
	}
	/**
	 * Swaps 2 samples.
	 *
	 * @param a
	 * @param b
	 */
	static inline void swap(sample_type& a, sample_type& b){
		scalar_type temp;
		temp = a.lambda;
		a.lambda = b.lambda;
		b.lambda = temp;

		temp = a.f_lambda;
		a.f_lambda = b.f_lambda;
		b.f_lambda = temp;
	};


	void minbrak_simple(sample_type& t1,sample_type& t2,sample_type& t3,sample_type& t4, const unsigned int bound);

	void minbrak_para_ext(sample_type& t1,sample_type& t2,sample_type& t3,sample_type& t4, const unsigned int bound);
	/**
	 * Guesses an initial sample of the function to start with the search algorithm.
	 *
	 * @return Initial pivot lambda
	 */
	std::list<scalar_type> get_init_sample() const;

	/**
	 * Plots the function under minimization. This is to get an idea visually as
	 * where is the minimum reached.
	 * @param a Initial sampling point
	 * @param b End sampling point
	 * @param delta Movement parameter from a to b.
	 */
	void plot_function(double a, double b, double delta) const {
		std::ofstream temp_file;
		temp_file.open("/tmp/func_samples.txt");
		for(double i=a; i<=b;i+=delta){
			temp_file << i << " " << my_functor.compute(scalar_type(i)) << std::endl;
		}
		temp_file.close();
		int res=std::system("graph -TX -BC /tmp/func_samples.txt");
	}

private:

	functor<scalar_type> my_functor;
	std::string my_minbrak_type;
	math::affine_map<scalar_type> my_dynamics_map; // The dynamics maps of a location.
	support_function_provider::const_ptr my_U; // input set
	scalar_type my_eps;
	// We add the dynamics map and the input set to convex opt class
	// becuase the min bracketing of the class depends on them.
};

} // end of namespace support_function
} // end of namespace continuous

/**
 * Finds the minimum of the two passed arguments.
 *
 * @param argument 1
 * @param argument 2
 * @return The minimum of the two arguments.
 */
template<typename scalar_type>
scalar_type get_minimum(scalar_type arg1,scalar_type arg2){
	math::numeric::approx_comparator<scalar_type> my_comparator;
	if(my_comparator.is_definitely_strictly_larger(arg1,arg2))
		return arg2;
	else
		return arg1;
}
/**
 * Finds the maximum of the two passed arguments.
 *
 * @param argument 1
 * @param argument 2
 * @return The max of the two arguments.
 */
template<typename scalar_type>
scalar_type get_maximum(scalar_type arg1,scalar_type arg2){
	math::numeric::approx_comparator<scalar_type> my_comparator;
	if(my_comparator.is_definitely_strictly_larger(arg1,arg2))
		return arg1;
	else
		return arg2;
}

template<typename scalar_type>
scalar_type SIGN(scalar_type a, scalar_type b){
	return b >= scalar_type(0) ? (a >= scalar_type(0) ? a : -a):(a >= scalar_type(0) ? -a : a);
}

template<typename scalar_type>
scalar_type MAX(const scalar_type& a, const scalar_type& b){
	return a > b? a : b;
}
#include "convex_optimization.hpp"


#endif /* _CONVEX_OPT */
