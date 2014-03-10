/*
 * affine_interpolation.h
 *
 *  Created on: Aug 12, 2011
 *      Author: frehse
 */

#ifndef AFFINE_INTERPOLATION_H_
#define AFFINE_INTERPOLATION_H_

#include "core/continuous/continuous_dynamics/typed_dynamics.h"
#include "core/continuous/polyhedra/simplex.h"
#include "core/continuous/polyhedra/hyperbox/bounding_box.h"
#include "math/vdom/affine_map.h"
#include "math/vdom/vdom_vector.h"
#include "math/vdom/lin_constraint_system.h"
#include "core/continuous/support_function/template_directions/choose_directions.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/predicates/valuation_function_tree_ops.h"
#include "math/numeric/boost_interval_utility.h"

#ifdef USE_NLTOOLBOX
	#include <nltoolbox/HybridizationInterface.h>
#endif

#include "extern/aaflib_loader.h"

namespace continuous {

/** Evaluates the dynamics of index i at point pt according to dynamics
 * in data.
 *
 * data must be a pointer to typed_dynamics<double>.
 */
inline
double dynamicsAdapter(const int i, const std::vector<double> &pt, void* data) {
	typed_dynamics<double>* dyn = static_cast<typed_dynamics<double>*>(data);
	const positional_vdomain& dom = dyn->codom();

	typed_dynamics<double>::state x(dom, math::vector<double>(pt));
	return dyn->compute_deriv(i, x);
}

/** Evaluates the Hessian of index i wrt to x_j and x_k at point pt according to dynamics
 * in data.
 *
 * data must be a pointer to typed_dynamics<double>.
 */
inline
double hessianAdapter(const int i, const std::vector<double> &pt, const int j,
		const int k, void* data) {
	typed_dynamics<double>* dyn = static_cast<typed_dynamics<double>*>(data);
	const positional_vdomain& dom = dyn->codom();

	typed_dynamics<double>::state x(dom, math::vector<double>(pt));
	return dyn->compute_hessian(i, j, k, x);
}

/** Transforms a vector to a std::vector */
inline std::vector<double> convert_to_std_vector(
		const math::vector<double>& v) {
	std::vector<double> b(v.size());
	for (math::vector<double>::size_type i = 0; i < v.size(); ++i) {
		b[i] = v[i];
	}
	return b;
}

/** Transforms a std::vector to a vector */
inline math::vector<double> convert_from_std_vector(
		const std::vector<double>& v) {
	math::vector<double> b(v.size());
	for (std::vector<double>::size_type i = 0; i < v.size(); ++i) {
		b[i] = v[i];
	}
	return b;
}

/** Transforms a matrix M to a vector of vectors */
inline std::vector<std::vector<double> > convert_to_vector_vector(
		const math::matrix<double>& M) {
	math::matrix<double>::size_type N1, N2;
	N1 = M.size1();
	N2 = M.size2();
	std::vector<std::vector<double> > A = std::vector<std::vector<double> >(N1,
			std::vector<double>(N2));
	for (math::matrix<double>::size_type i = 0; i < N1; ++i) {
		for (math::matrix<double>::size_type j = 0; j < N2; ++j) {
			A[i][j] = M(i, j);
		}
	}
	return A;
}

/** Transforms a vector of vectors to a matrix */
inline math::matrix<double> convert_to_matrix(
		const std::vector<std::vector<double> > &A) {
	math::matrix<double>::size_type N1, N2;
	N1 = A.size();
	math::matrix<double> M;
	if (N1 > 0) {
		N2 = A[0].size();
		M = math::matrix<double>(N1, N2);
		for (math::matrix<double>::size_type i = 0; i < N1; ++i) {
			for (math::matrix<double>::size_type j = 0; j < N2; ++j) {
				M(i, j) = A[i][j];
			}
		}
	} else {
		M = math::matrix<double>();
	}
	return M;
}

/** Transform a set of constraints to NLToolbox vectors */
inline
void convert_to_vectors(const math::lin_constraint_system<double>& cons,
		unsigned& dim, unsigned& nbConstraint,
		std::vector<std::vector<double> > &A, std::vector<double> &b) {
	// get matrices from constraints
	math::matrix<double> MA;
	math::vector<double> Mb;
	positional_vdomain dom;
	canonic_matrix_form(MA, Mb, dom, cons);

	// get vectors from matrices
	A = convert_to_vector_vector(MA);
	b = convert_to_std_vector(Mb);
	dim = dom.size();
	nbConstraint = MA.size1();
	if (dim != MA.size2()) {
		std::cerr << "Domain: " << dom << ", dimensions: " << dim;
		throw std::runtime_error(
				"non-matching dimensions in " + std::string(__FUNCTION__));
	}
}

/** Transform from NLToolbox vectors to a set of constraints */
inline math::lin_constraint_system<double> constraints_from_vectors(
		const positional_vdomain& dom,
		const std::vector<std::vector<double> > &A,
		const std::vector<double> &b) {
	math::lin_constraint_system<double> cons;

	math::matrix<double>::size_type N1, N2;
	N1 = A.size();
	if (N1 > 0) {
		N2 = A[0].size();
		if (dom.size() != N2) {
			throw std::runtime_error(
					"non-matching dimensions in " + std::string(__FUNCTION__));
		}
		// for each constraint
		for (math::matrix<double>::size_type i = 0; i < N1; ++i) {
			math::vector<double> v(A[i]);
			math::vdom_vector<double> v_dom(dom, v);
			double c = b[i];
			math::lin_constraint<double> con(v_dom, -c, LE);
			cons.push_back(con);
		}
	}
	return cons;
}

/** Transform from NLToolbox dynamics to dynamics */
inline ode_affine_dynamics<double> dynamics_from_vectors(const positional_vdomain& dom,const std::vector<std::vector<double> > &A,
		const std::vector<double> &b) {
	math::affine_map<double> M(dom,convert_to_matrix(A),math::vector<double>(b));
	return ode_affine_dynamics<double>(M);
}

/** Try to retrieve constraints from a support_function_provider */
template<class scalar_type>
math::lin_constraint_system<scalar_type> compute_constraints(
		const support_function_provider& X0, const positional_vdomain& dom) {
	using namespace support_function;
	bool remove_redundant = true;
	if (const polyhedron<scalar_type>* p = dynamic_cast<const polyhedron<
			scalar_type>*>(&X0)) {
		if (remove_redundant) {
			constr_polyhedron<scalar_type> poly;
			poly.add_constraints(*p->get_constraints());
			poly.remove_redundant_constraints();
			return *poly.get_constraints();
		}
		return *p->get_constraints();
	} else {
		typedef std::list<math::vector<double> > direction_list_type;
		direction_list_type dir_list =
				direction_chooser::get_directions<double>(dom);
		math::vdom_vector<scalar_type> d(dom);
		typedef std::set<math::vdom_vector<scalar_type>,
				math::numeric::lex_comp_less<scalar_type, math::vdom_vector> > direction_set;
		direction_set dirs;
		for (direction_list_type::const_iterator it = dir_list.begin();
				it != dir_list.end(); ++it) {
			d.set_vector(*it);
			dirs.insert(d);
		}
		constr_polyhedron<scalar_type> poly = compute_outer_poly(X0, dirs);
		if (remove_redundant)
			poly.remove_redundant_constraints();
		return *poly.get_constraints();
	}
}

/** Evaluate a bounding box of the given ode system */
template<class scalar_type>
continuous::finite_hyperbox<scalar_type> derivative_bounds(const continuous::typed_dynamics<scalar_type>& dyn, const support_function_provider& domain) {
	using namespace continuous;
	finite_hyperbox<scalar_type> domain_box=finite_bounding_box<scalar_type>(domain);

	if (domain_box.is_empty()) {
		return finite_hyperbox<scalar_type>::empty_box(domain_box.domain());
	}

	if (const ode_dynamics<scalar_type>* ode_dyn = dynamic_cast<const ode_dynamics<scalar_type>*>(&dyn)){
		// use interval arithmetic

		//typedef boost::numeric::interval<scalar_type> bound_interval;
		//typedef boost::numeric::interval<scalar_type> interval_type;
		//typedef typename boost::numeric::interval_lib::unprotect<interval_type >::type bound_interval;
//		using namespace boost::numeric::interval_lib;
//		typedef boost::numeric::interval<double,
//		                   policies<save_state<rounded_transc_opp<double> >,
//		                            checking_base<double> > > bound_interval;
		typedef AAF bound_interval;

		// convert the domain_box to an interval vector
		const positional_vdomain& dom = domain_box.get_g_dom().domain();
		math::vdom_vector<bound_interval> arg_vec(dom),deriv_vec(dom);
		{
			// create domain intervals
			const math::vector<scalar_type> c = domain_box.get_c();
			const math::vector<scalar_type> g = domain_box.get_g();
			for (unsigned int i = 0; i < dom.size(); ++i) {
//				arg_vec[i] = bound_interval(c[i] - g[i], c[i] + g[i]);
				arg_vec[i] = bound_interval(AAInterval(c[i] - g[i], c[i] + g[i]));
			}
		}

		typename math::vdom_vector<bound_interval>::iterator dit = deriv_vec.begin();
		for (typename ode_dynamics<scalar_type>::const_iterator it = ode_dyn->begin(); it!= ode_dyn->end(); ++it,++dit) {
			bound_interval intv = valuation_functions::arithmetic_eval_node<bound_interval>(*it,arg_vec);
			*dit = intv; // store in deriv_vec
		}
		// convert deriv_vec to hyperbox
		math::vector<scalar_type> c(dom.size());
		math::vector<scalar_type> g(dom.size());
		for (unsigned int i = 0; i < dom.size(); ++i) {
			const bound_interval& d = deriv_vec[i];
//			const scalar_type& l = d.lower();
//			const scalar_type& u = d.upper();
			const scalar_type& l = d.getMin();
			const scalar_type& u = d.getMax();

			c[i] = (l + u) / scalar_type(2);
			g[i] = u - c[i];
		}
		return finite_hyperbox<scalar_type>(c, g, dom);
	} else if (const ode_affine_dynamics<scalar_type>* ode_dyn = dynamic_cast<const ode_affine_dynamics<scalar_type>*>(&dyn)) {
		return finite_bounding_box(domain,*ode_dyn);
	} else {
		throw std::runtime_error("derivative_bounds: missing implementation");
	}
	return finite_hyperbox<scalar_type>();
}

/** Evaluate a bounding box of the given equation system */
template<class scalar_type>
continuous::finite_hyperbox<scalar_type> derivative_err_bounds(const continuous::typed_dynamics<scalar_type>& dyn, const math::affine_map<scalar_type>& M, const support_function_provider& domain) {
	using namespace continuous;
	using namespace valuation_functions;

	finite_hyperbox<scalar_type> domain_box=finite_bounding_box<scalar_type>(domain);

	if (domain_box.is_empty()) {
		return finite_hyperbox<scalar_type>::empty_box(domain_box.domain());
	}

	if (const ode_dynamics<scalar_type>* ode_dyn = dynamic_cast<const ode_dynamics<scalar_type>*>(&dyn)){
		// use interval arithmetic

		//typedef boost::numeric::interval<scalar_type> bound_interval;
		//typedef boost::numeric::interval<scalar_type> interval_type;
		//typedef typename boost::numeric::interval_lib::unprotect<interval_type >::type bound_interval;
//		using namespace boost::numeric::interval_lib;
//		typedef boost::numeric::interval<double,
//		                   policies<save_state<rounded_transc_opp<double> >,
//		                            checking_base<double> > > bound_interval;
		typedef AAF bound_interval;

		// convert the domain_box to an interval vector
		const positional_vdomain& dom = domain_box.get_g_dom().domain();
		math::vdom_vector<bound_interval> arg_vec(dom),deriv_vec(dom);
		{
			// create domain intervals
			const math::vector<scalar_type> c = domain_box.get_c();
			const math::vector<scalar_type> g = domain_box.get_g();
			for (unsigned int i = 0; i < dom.size(); ++i) {
				scalar_type low = c[i] - g[i];
				scalar_type high = c[i] + g[i];
//				arg_vec[i] = bound_interval(low,high);
				arg_vec[i] = bound_interval(AAInterval(low,high));
			}
		}

		typename math::vdom_vector<bound_interval>::iterator dit = deriv_vec.begin();
		size_t index = 0;
		for (typename ode_dynamics<scalar_type>::const_iterator it = ode_dyn->begin(); it!= ode_dyn->end(); ++it,++dit,++index) {
			// obtain the function M[i] - *it
			tree::node::ptr aff_node = scalar_product_node<scalar_type>(M.get_A().vector_from_row(index),M.get_b()[index],dom);
			tree::node::ptr err_node = node_subtraction(aff_node,*it);
			bound_interval intv = arithmetic_eval_node<bound_interval>(err_node,arg_vec);
			*dit = intv; // store in deriv_vec
		}
		// convert deriv_vec to hyperbox
		math::vector<scalar_type> c(dom.size());
		math::vector<scalar_type> g(dom.size());
		for (unsigned int i = 0; i < dom.size(); ++i) {
			const bound_interval& d = deriv_vec[i];
//			std::cout << "derivative error for " << dom.get_variable(i) <<": "<< d << std::endl;
//			const scalar_type& l = d.lower();
//			const scalar_type& u = d.upper();
			const scalar_type& l = d.getMin();
			const scalar_type& u = d.getMax();
			c[i] = (l + u) / scalar_type(2);
			g[i] = u - c[i];
		}
		return finite_hyperbox<scalar_type>(c, g, dom);
	} else {
		throw std::runtime_error("derivative_bounds: missing implementation");
	}
	return finite_hyperbox<scalar_type>();
}

struct affine_interpolation_result {
	typedef enum {
		SUCCESS, REQUIRESPLITTING
	} type;
};

/** Return the min of the two, not counting negative values */
template<typename scalar_type>
const scalar_type& pos_min(const scalar_type& x, const scalar_type& y) {
	if (x<scalar_type(0) || (y>=scalar_type(0) && y<x)) {
		return y;
	}
	return x;
}

/** Return the min of the two, not counting negative values */
template<typename scalar_type>
math::vector<scalar_type> pos_min(const math::vector<scalar_type>& x, const math::vector<scalar_type>& y) {
	if (x.size()<=0 || x.size()!=y.size())
		return math::vector<scalar_type>();

	math::vector<scalar_type> z(x);
	for (typename math::vector<scalar_type>::size_type i = 0; i<x.size(); ++i) {
		if (x[i]<scalar_type(0) || (y[i]>=scalar_type(0) && y[i]<x[i])) {
			z[i] = y[i];
		}
	}
	return z;
}

/** An estimate of max dwell time of states X0 inside bounds based on bounding boxes
 *
 * The dynamics are given by dyn + U. If no bounds can be established,
 * the result is negative. */
template<typename scalar_type>
scalar_type estimate_max_dwell_time(const typed_dynamics<scalar_type>& dyn,
		const support_function_provider::const_ptr& U,
		const support_function_provider& X0,
		const support_function_provider& bounds) {
	using namespace math;
	using namespace math::numeric;
	using namespace continuous;

	typedef typename hyperbox<scalar_type>::point_type point;
	typedef typename hyperbox<scalar_type>::value_type scalar_inf;
	typedef math::vdom_vector<scalar_inf> vdom_point;

	scalar_type T(-1); // bound on dwell time

	// bounds on the derivatives
	hyperbox<scalar_type> deriv_bounds;
	if (const ode_affine_dynamics<scalar_type>* ode_dyn =
			dynamic_cast<const ode_affine_dynamics<scalar_type>*>(&dyn)) {
		deriv_bounds = compute_bounding_box<scalar_type>(bounds, *ode_dyn);
//		std::cout << "derivative bounds: " << deriv_bounds << " from " << bounds << std::endl;
	} else {
		deriv_bounds = construct_hyperbox<scalar_type>(derivative_bounds(dyn, bounds));
	}

	// add offset
	if (U) {
//		std::cout << "computing U bounds: " << std::endl;
		deriv_bounds += compute_bounding_box<scalar_type>(*U);
//		std::cout << "deriv + U bounds: " << deriv_bounds << std::endl;
	}

	// bounds on the domain
	hyperbox<scalar_type> X0_bounds = compute_bounding_box<scalar_type>(X0);
//	std::cout << "X0 bounds: " << X0_bounds << std::endl;
	hyperbox<scalar_type> inv_bounds = compute_bounding_box<scalar_type>(bounds);
//	std::cout << "inv bounds: " << inv_bounds << std::endl;

	// we need to bring deriv_bounds, X0_bounds and inv_bounds to the same domain
	vdom_point lD(deriv_bounds.domain(), deriv_bounds.get_l());
	vdom_point uD(deriv_bounds.domain(), deriv_bounds.get_u());
	vdom_point lX(X0_bounds.domain(), X0_bounds.get_l());
	vdom_point uX(X0_bounds.domain(), X0_bounds.get_u());
	vdom_point lI(inv_bounds.domain(), inv_bounds.get_l());
	vdom_point uI(inv_bounds.domain(), inv_bounds.get_u());
	// get the common domain
	positional_vdomain all_dom;
	all_dom = compose(deriv_bounds.domain(),compose(X0_bounds.domain(),inv_bounds.domain()));
	// remap everything, using the appropriate default values
	// for newly embedded variables
	lD.reorder(all_dom,scalar_inf::neg_infty());
	uD.reorder(all_dom,scalar_inf::pos_infty());
	lX.reorder(all_dom,scalar_inf::neg_infty());
	uX.reorder(all_dom,scalar_inf::pos_infty());
	lI.reorder(all_dom,scalar_inf::neg_infty());
	uI.reorder(all_dom,scalar_inf::pos_infty());
//	std::cout << "lX:" << lX << "uI:" << uI << std::endl;

	// size of all vectors
	typename math::vector<scalar_type>::size_type N = lX.size();

	scalar_inf zero(0);
	scalar_inf Tinf = scalar_inf::pos_infty();

	// in each coordinate, find if there's a bound
	// a) for positive deriv, it's uI-lX divided by lower absolute speed
	vdom_point pos_width=uI-lX;
	for (typename math::vector<scalar_type>::size_type i = 0; i<N; ++i) {
//		std::cout << "for " << all_dom.get_variable(i) << pos_width[i] << "," << lD[i] << "->" <<  pos_width[i]/lD[i] << std::endl;
		if (pos_width[i].is_finite() && pos_width[i].is_finite() && definitely(is_GT(lD[i],zero))) {
			Tinf = std::min(Tinf,pos_width[i]/lD[i]);
		}
	}
	// b) for negative deriv, it's lI-uX divided by lower absolute speed
	vdom_point neg_width=lI-uX;
	for (typename math::vector<scalar_type>::size_type i = 0; i<N; ++i) {
		if (definitely(is_LT(uD[i],zero))) {
			Tinf = std::min(Tinf,neg_width[i]/uD[i]);
		}
	}
	if (Tinf.is_finite()) {
		T = Tinf.get_val();
	}

	LOGGER_OS(DEBUG5,__FUNCTION__) << "estimated max dwell time: " << T;
	return T;
}

/** Compute an affine interpolation of a nonlinear function
 *
 *	Uses the NLtoolbox by Romain Testylier
 *
 *	Takes nonlinear dynamics dyn and an initial set of states X0.
 *	Returns the approximating affine dynamics M plus a
 *	bound on the approximation error U, and the
 *	validity domain as domain_cons.
 */
template<class scalar_type>
affine_interpolation_result::type affine_interpolation(
		const typed_dynamics<scalar_type>& dyn,
		const support_function_provider& X0,
		math::affine_map<scalar_type>& M,
		support_function_provider::const_ptr& U,
		math::lin_constraint_system<scalar_type>& domain_cons,
		double precision
		) {
	affine_interpolation_result::type result;

#ifdef USE_NLTOOLBOX
	const positional_vdomain& dom = dyn.codom();
	size_t dim = dom.size();
	nltool::HybridizationInterface hybInterface(dim);
	void* dyn_pointer = (void*)(&dyn);
	hybInterface.setSystemDynamic(dim, dynamicsAdapter, hessianAdapter, dyn_pointer);

	// initial set
	size_t nbConstraint;
	math::lin_constraint_system<double> cons = compute_constraints<scalar_type>(X0,dom);
	std::vector<std::vector<double> > initialA;
	std::vector<double> initialB;
	convert_to_vectors(cons,dim,nbConstraint,initialA,initialB);
	// set
	hybInterface.setInitialSet(dim, nbConstraint, initialA, initialB);

	// compute approximation
	std::vector<std::vector<double> > approxA;
	std::vector<double> approxB;
	std::vector<std::vector<double> > domainA;
	std::vector<double> domainB;
	// bounds on the approximation error (hyperbox vertices)
	std::vector<double> err_low;
	std::vector<double> err_up;

    // redirect  NLToolbox output
    std::stringstream ob;
    nltool::HybridizationResult nltoolres;
	{
		stream_redirector redirect(std::cout, ob);
//	nltool::HybridizationResult nltoolres = hybInterface.computeApproximation(precision, approxA, approxB,
//			domainA, domainB);
		nltoolres = hybInterface.computeApproximation(precision, approxA,
				approxB, domainA, domainB, err_low, err_up);
		LOGGER_OS(DEBUG6,__FUNCTION__) << ob.str();
	}

	if (nltoolres==nltool::SUCCESS) {
		result = affine_interpolation_result::SUCCESS;

		// obtain the dynamics and domain
		continuous::ode_affine_dynamics<double> approx_dyn = dynamics_from_vectors(dom,approxA,approxB);
//		std::cout << "approximation dynamics: " << approx_dyn << std::endl;

		M = approx_dyn;

		// the following version is obsolete with NLtoolbox update from 2013-05-29
//		typename finite_hyperbox<scalar_type>::point_type center(dim,scalar_type(0));
//		typename finite_hyperbox<scalar_type>::point_type generator(dim,scalar_type(precision));
//		U = support_function_provider::const_ptr(new finite_hyperbox<scalar_type>(center,generator,dom));
//		std::cout << "input set: " << finite_hyperbox<scalar_type>(center,generator,dom) << std::endl;

		math::vector<double> e_low =  convert_from_std_vector(err_low);
		math::vector<double> e_up =  convert_from_std_vector(err_up);
		typename finite_hyperbox<scalar_type>::point_type center = (e_low+e_up)/scalar_type(2);
		typename finite_hyperbox<scalar_type>::point_type generator = e_up - center;;
	    LOGGER_OS(DEBUG6,__FUNCTION__) << "requested precision "+to_string(precision)+ ", got "+to_string(generator);

		U = support_function_provider::const_ptr(new finite_hyperbox<scalar_type>(center,generator,dom));

		// obtain the approximation domain
		domain_cons = constraints_from_vectors(dom,domainA,domainB);
//		std::cout << "approximation domain: " << domain_cons << std::endl;

	} else if (nltoolres==nltool::REQUIRESPLITTING) {
		result = affine_interpolation_result::REQUIRESPLITTING;
	} else {
		throw std::runtime_error("unkown result from NLToolbox");
	}
#else
	throw std::runtime_error("missing NLToolbox (recompile required)");
#endif

	return result;
}
/** Compute an affine interpolation of a nonlinear function on the vertices of a simplex
 *
 */
/*
template<class scalar_type>
ode_affine_dynamics<scalar_type> affine_interpolation(
		const simplex<scalar_type>& s,
		const typed_dynamics<scalar_type>& nonLinearSystem) {
	// asserts:
	// The number of equations needs to be the same as the dimension of the simplex.

	int dim = s.get_dim();
	bool singular;
	math::matrix<scalar_type> A((dim + 1) * dim, (dim + 1) * dim);
	const typename simplex<scalar_type>::point_set_type & vertices =
			s.get_vertices();
	math::vector<scalar_type> b((dim + 1) * dim);
	int i = 0;
	//point_set_type it=vertices.begin();
	typename simplex<scalar_type>::point_set_type::const_iterator it;
	// for each vertex (treat vertex number i)
	for (it = vertices.begin(); it != vertices.end(); it++) {
		math::vdom_vector<scalar_type> state_i(s.domain(), it->first);
		math::vdom_vector<scalar_type> deriv_i = nonLinearSystem.compute_deriv(
				state_i);
		// fill dim rows of the matrix A and dim elements of b
		for (int j = 0; j < dim; j++) {
			b[i * dim + j] = deriv_i[j];
			for (int k = 0; k < dim; k++) {
				A(i * dim + j, k + j * (dim + 1)) = state_i[k];
			}
			A(i * dim + j, dim + j * (dim + 1)) = 1;
		}
		++i;
	}

	// Compute A^-1*b
	math::vector<scalar_type> result(A.inverse(singular) * b);
	if (singular) {
		// throw exception
		std::stringstream ss;
		ss << A;
		throw std::runtime_error(
				"Singular matrix cannot be inversed:\n" + ss.str());
	}

	// parse the result
	math::matrix<scalar_type> B(dim, dim);
	math::vector<scalar_type> c(dim);
	for (int j = 0; j < dim; j++) {
		for (int k = 0; k < dim; k++) {
			B(j, k) = result[j * (dim + 1) + k];
		}
		c[j] = result[j * (dim + 1) + dim];
	}

	// create the affine interpolation

	math::affine_map<scalar_type> affineMap(s.domain(), B, c);
	ode_affine_dynamics<scalar_type> affineSystem(affineMap);
	return affineSystem;
}
;
*/

}

#endif /* AFFINE_INTERPOLATION_H_ */
