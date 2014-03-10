/**
 * Functions to compute the support function of a set, on which support function is defined intersected, with a
 * hyperplane
 *
 * @author Rajarshi
 */
#ifndef SF_PLANE_INTERSECTION_H_
#define SF_PLANE_INTERSECTION_H_

#include "core/continuous/support_function_provider.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "math/vdom/affine_map.h"
#include "math/vdom/vdom_vector.h"
#include "math/vdom/lin_constraint.h"
#include "math/numeric/interval.h"
#include "core/continuous/support_function/intersection_by_optimization/lower_bound_search.h"
#include "core/continuous/support_function/intersection_by_optimization/lb_search_opt.h"
#include <boost/pointer_cast.hpp>
#include "core/continuous/support_function/sf_base/sf_set.h"
#include "io/GEN_format/GEN_formatter.h"
#include "core/continuous/support_function/sfm/sfm_guard_intersection.h"

//#include "core/continuous/support_function/sf_base/sf_set.h"

/**
 * Computes the support function of support_function_provider set S intersection hyperplane G.
 * This function assumes that the intersection is non-empty.
 *
 * 1. Create a convex_opt object with support_functor
 * 2. Find out the transformation to shift the given hyperplane to origin. Apply the same transformation
 *    to the support function provider set.
 * 3. Apply the minimization algorithm to the shifted polytope and hyperplane. Get the minimum.
 * 4. Get the correct min for un shifted set and hyperplane
 *
 *
 * @param S continuous set on which support function operation is defined
 * @param l the direction in which to compute the support function
 * @param max_value the support function value at direction l
 * @param support_vec the support vector at direction l
 * @param is_empty flag to denote if S intersect G in empty
 * @param is_bounded flag to denote if S intersect G is bounded
 * @param G Hyperplane given as a linear constraint
 */
namespace continuous {
namespace support_function {


//------------------------------------
inline void compute_support_intersection(
		const support_function_provider::const_ptr& S,
		const math::vdom_vector<double>& l_dom, double& max_value,
		bool& is_empty, bool& is_bounded, const math::lin_constraint<double>& G,
		const std::string& minbrak_type, const double& intersection_error) {
	// This method assumes that the intersection is non-empty and bounded
	is_empty = false;
	is_bounded = true;

	positional_vdomain set_dom =  positional_vdomain(S->get_variable_ids());

	double G_inh_term = G.get_canonic_inh_coeff();
	comparison_operator my_sign = G.get_sign();

	math::vdom_vector<double> G_n = G.get_normal();
	G_n.reorder(set_dom);

	//translate S to S'
	math::vdom_vector<double> shift_vector = get_shift_vector<double>(G,set_dom);
	continuous::support_function_provider::ptr shifted_set_ptr = shift_set<double>(S,shift_vector);

	// check on the optimizing conditions
	bool to_activate = false;
	check_opt_conditions<double>(S,G,l_dom,is_empty,is_bounded,to_activate, max_value);

	if(is_empty || !is_bounded)
		return;
	if(to_activate){
		// Create a convex_opt object and pass the shifted set and guard as arguments
		continuous::support_function::sup_functor<double> my_sup_function(l_dom,
				G_n, G_inh_term, shifted_set_ptr);
		continuous::support_function::convex_opt<double,
				support_function::sup_functor>* my_convex_prb =
				new continuous::support_function::lower_bound_search<double,
						support_function::sup_functor>(my_sup_function, minbrak_type, intersection_error);

		math::numeric::interval<double> my_interval;

		if(my_sign == EQ){
			my_interval = math::numeric::interval<double>(
					scalar_with_infinity<double>::neg_infty(),
					scalar_with_infinity<double>::pos_infty());
		}
		else if(my_sign == LE){
			my_interval = math::numeric::interval<double>(
							scalar_with_infinity<double>(0),
							scalar_with_infinity<double>::pos_infty());
		}
		else{
			std::runtime_error("compute_support_intersection: Expects guard in canonical form\n");
		}

		support_function::convex_opt<double, support_function::sup_functor>::min_interval
				min_structure;
		min_structure = my_convex_prb->section_search(my_interval);
//		my_convex_prb->plot_function(0.2,0.8,0.001);
//		my_convex_prb->plot_graph();
		adjust_shifting<double>(min_structure.f_lambda_interval,l_dom,shift_vector);
		max_value = min_structure.f_lambda_interval.upper().get_val();
	}

	LOGGER(DEBUG7,"compute_support_intersection" ,"sup for direction:"+to_string(l_dom)+" is = "+to_string(max_value));
//		std::cout << "lambda interval with min: " << min_structure.lambda_interval.lower().get_val() << ", " << min_structure.lambda_interval.upper().get_val() << std::endl;
//		std::cout << "f_lambda interval with min: " << min_structure.f_lambda_interval.lower().get_val() << ", " << min_structure.f_lambda_interval.upper().get_val() << std::endl;
//		std::cout << "trans vector:" << trans_vector << ", " << "l_dom:" << l_dom << ", scal_prod:" << scalar_product(trans_vector, l_dom) << std::endl;

}

inline void compute_support_intersection(
		const support_function_provider::const_ptr& S,
		const math::vdom_vector<double>& l_dom, double& max_value,
		bool& is_empty, bool& is_bounded, const math::lin_constraint<double>& G,
		const unsigned int sample_bound,
		const std::string& minbrak_type, const double& intersection_error) {
	// This method assumes that the intersection is non-empty and bounded
	is_empty = false;
	is_bounded = true;

	positional_vdomain set_dom =  positional_vdomain(S->get_variable_ids());

	double G_inh_term = G.get_canonic_inh_coeff();
	comparison_operator my_sign = G.get_sign();

	math::vdom_vector<double> G_n = G.get_normal();
	G_n.reorder(set_dom);

	//translate S to S'
	math::vdom_vector<double> shift_vector = get_shift_vector<double>(G,set_dom);
	continuous::support_function_provider::ptr shifted_set_ptr = shift_set<double>(S,shift_vector);

	// check on the optimizing conditions
	bool to_activate = false;
	double max_val;
	check_opt_conditions(S,G,l_dom,is_empty,is_bounded,to_activate, max_value);

	if(is_empty || !is_bounded)
		return;
	if(to_activate){
		// Create a convex_opt object and pass the shifted set and guard as arguments
		continuous::support_function::sup_functor<double> my_sup_function(l_dom,
				G_n, G_inh_term, shifted_set_ptr);
		continuous::support_function::lower_bound_search<double,
				support_function::sup_functor>* my_convex_prb =
				new continuous::support_function::lower_bound_search<double,
						support_function::sup_functor>(my_sup_function, minbrak_type, intersection_error);

		math::numeric::interval<double> my_interval;

		if(my_sign == EQ){
			my_interval = math::numeric::interval<double>(
					scalar_with_infinity<double>::neg_infty(),
					scalar_with_infinity<double>::pos_infty());
		}
		else if(my_sign == LE){
			my_interval = math::numeric::interval<double>(
							scalar_with_infinity<double>(0),
							scalar_with_infinity<double>::pos_infty());
		}
		else{
			std::runtime_error("compute_support_intersection: Expects guard in canonical form\n");
		}
		//std::cout << "lbs tester: before call\n";
		support_function::convex_opt<double, support_function::sup_functor>::min_interval
				min_structure;
		min_structure = my_convex_prb->section_search_bounded(my_interval,sample_bound);

		adjust_shifting<double>(min_structure.f_lambda_interval,l_dom,shift_vector);
		max_value = min_structure.f_lambda_interval.upper().get_val();
	}
}

inline void compute_support_intersection(
		const support_function_provider::const_ptr& S,
		const math::vdom_vector<Rational>& l_dom, Rational& max_value,
		bool& is_empty, bool& is_bounded, const math::lin_constraint<Rational>& G,
		const std::string minbrak_type, const Rational& intersection_error) {
	// This method assumes that the intersection is non-empty and bounded

	Rational G_inh_term = G.get_canonic_inh_coeff();
	Rational trans_vec_norm;
	Rational G_n_norm(0);
	comparison_operator my_sign = G.get_sign();

	positional_vdomain set_dom = positional_vdomain(S->get_variable_ids());

	math::lin_expression<Rational> le = G.get_canonic_l();

	math::vdom_vector<Rational> G_n = G.get_normal();
	G_n.reorder(set_dom);


	/* Optimization:
	 *
	 */
	math::vdom_vector<Rational> sup_vec;

	/* Check the condition:
	 * 1: G is a hyperplane
	 * 2: l is a multiple o G_n
	 * then return -G_inh_term as the support function value
	 */
	Rational n_norm = G_n.infinity_norm();
	Rational l_norm = l_dom.infinity_norm();
	Rational f = n_norm/l_norm;

	math::vdom_vector<Rational> d = f * l_dom;
	math::vdom_vector<Rational> d_neg = -d;


	if (!math::numeric::is_MEQ(f,Rational(0))){
		if (math::numeric::is_MEQ(G_n.begin(), G_n.end(), d.begin(),d.end())) {// l and n are in the same direction
			max_value = -G_inh_term / f;
			LOGGER(DEBUG7,"compute_support_intersection:" ,"sup for direction:"+to_string(l_dom)+" is = "+to_string(max_value));
			return;
		}
		else if (math::numeric::is_MEQ(G_n.begin(), G_n.end(), d_neg.begin(),d_neg.end())) { // l and n are in the opposite direction
			if(my_sign == EQ){
				max_value =  G_inh_term / f ; // return - inh_term of the guard
				LOGGER(DEBUG7,"compute_support_intersection:" ,"sup for direction:"+to_string(l_dom)+" is = "+to_string(max_value));
				return;
			}
			else{
				S->compute_support(l_dom,max_value,sup_vec,is_empty,is_bounded);
				LOGGER(DEBUG7,"compute_support_intersection" ,"sup for direction:"+to_string(l_dom)+" is = "+to_string(max_value));
				return;
			}
		}
		else{} // Continue with the lower bound search below
	}

	for (math::vector<Rational>::const_iterator it = G_n.begin(); it != G_n.end(); it++) {
		G_n_norm = G_n_norm + (*it) * (*it);
	}


	trans_vec_norm = G_inh_term / G_n_norm;
	math::vdom_vector<Rational> trans_vector(G_n);
	trans_vector = trans_vec_norm * trans_vector;

	// Create affine map for the translation

	math::affine_map<Rational> my_map(set_dom, trans_vector);
	//	std::cout << "my map:" << my_map << std::endl;

	//debug: Lets print the set before applying shift
//	std::ofstream myfile;
//	myfile.open("/tmp/set_orig.txt");
//	io::GEN_formatter formatter(myfile,0.0);
//	formatter.output(*S);
//	myfile.close();

	continuous::support_function_provider::const_ptr shifted_set_ptr(
			new continuous::support_function::sf_unary<Rational>(S,
					my_map));

	// Create a convex_opt object and pass the shifted set and guard as arguments
	continuous::support_function::sup_functor<Rational> my_sup_function(l_dom,
			G_n, G_inh_term, shifted_set_ptr);
	continuous::support_function::convex_opt<Rational,
			support_function::sup_functor>* my_convex_prb =
			new continuous::support_function::lower_bound_search<Rational,
					support_function::sup_functor>(my_sup_function, minbrak_type, intersection_error);


	math::numeric::interval<Rational> my_interval;
	if(my_sign == EQ){
		my_interval = math::numeric::interval<Rational>(
				scalar_with_infinity<Rational>::neg_infty(),
				scalar_with_infinity<Rational>::pos_infty());
	}
	else if(my_sign == LE){
/*		my_interval = math::numeric::interval<double>(
				scalar_with_infinity<double>::neg_infty(),
				scalar_with_infinity<double>(0)); */
		my_interval = math::numeric::interval<Rational>(
						scalar_with_infinity<Rational>(0),
						scalar_with_infinity<Rational>::pos_infty());
	}
	else{
		std::runtime_error("compute_support_intersection: Expects guard in canonical form\n");
	}
	//std::cout << "lbs tester: before call\n";
	support_function::convex_opt<Rational, support_function::sup_functor>::min_interval
			min_structure;
	min_structure = my_convex_prb->section_search(my_interval);

	max_value = min_structure.f_lambda_interval.upper().get_val()
			- scalar_product(trans_vector,l_dom);
	LOGGER(DEBUG7,"compute_support_intersection" ,"sup for direction:"+to_string(l_dom)+" is = "+to_string(max_value));

}

/**
 * Checks if the intersection between the support_function_provider set S and Hyperplane
 * given as linear contraint G is empty or not
 *
 * @param S continuous set on which support function operation is defined.
 * @param G A hyperplane given as a linear constraint.
 * @return true if S \cap G is empty. False otherwise.
 */
template<class T>
bool is_intersection_empty(const support_function_provider::const_ptr& S,
		const math::lin_constraint<T>& G) {
	using namespace math;
	using namespace math::numeric;

	bool is_empty, is_max_bounded, is_min_bounded;
	T max_val_pos, max_val_neg, G_inh_term;

	/* Check the constraint sign for equality */
	comparison_operator my_sign = G.get_sign();
	// Lets check the trivialities first
	if (S->is_empty() || !math::maybe(G.is_satisfiable())) {
		//debug
		//std::cout << "S is empty and G is unsat\n";
		return true; // intersection with empty set is empty always!
		//		std::cout << "intersection with empty set" << std::endl;
	}
	if (math::definitely(G.is_always_satisfied())){ // If one or both is universe and the other is not empty, the intersection is
		// not empty.
		return false;
	}
	if(my_sign != EQ)
		return inequality_empty_check(S,G);

	math::vdom_vector<T> sv; // support vector (with iimap of S)
	math::vdom_vector<T> le = G.get_normal();

	S->compute_support(le, max_val_pos, sv, is_empty, is_max_bounded);
	S->compute_support(-le, max_val_neg, sv, is_empty, is_min_bounded);
//		std::cout << "max_val_pos=" << max_val_pos << std::endl;
//		std::cout << "max_val_neg=" << max_val_neg << std::endl;
//		std::cout << "G inh term=" << -G.get_canonic_inh_coeff() << std::endl;

	T min_val = -max_val_neg;
	G_inh_term = -G.get_canonic_inh_coeff();

	// canonic form for LEQ side : le^T x LEQ G_inh_term

	// here equality: it's empty if min > G_inh_term or max < G_inh_term
	if (is_min_bounded && definitely(is_GT(min_val, G_inh_term))) {
		return true;
	}
	if (is_max_bounded && definitely(is_LT(max_val_pos, G_inh_term))) {
		return true;
	}
	return false;
}
/**
 * Checks if the intersection between the support_function_provider set S and Guard constraint
 * G in canonical form as well as not equality is empty or not.
 *
 * @param S continuous set on which support function operation is defined.
 * @param G A linear constraint not equality.
 * @return true if S \cap G is empty. False otherwise.
 */
template<class T>
bool inequality_empty_check(const support_function_provider::const_ptr& S,
			const math::lin_constraint<T>& G){
	using namespace math;

	bool is_empty, is_max_bounded, is_min_bounded;
	T max_val_pos, max_val_neg, G_inh_term;

	math::vdom_vector<T> sv; // support vector (with iimap of S)
	math::vdom_vector<T> le = G.get_normal();
	comparison_operator my_sign = G.get_canonic_sign();
	G_inh_term = -G.get_canonic_inh_coeff();

	S->compute_support(-le, max_val_neg, sv, is_empty, is_min_bounded);
	T min_val = -max_val_neg;

	// canonic form for LEQ side : le^T x LEQ G_inh_term
	// here equality: it's empty if min > G_inh_term

	// note: we don't care about strict inequality, since compute support returns inf anyway
	if (definitely(numeric::is_GT(min_val, G_inh_term))) {
		return true;
	} else {
		return false;
	}
}

} // end of namespace support_function
} // end of namespace continuous

#endif /* SF_PLANE_INTERSECTION_H_ */
