/*
 * vertice_approx_2D.h
 *
 *  Created on: Jan 19, 2010
 *      Author: frehse
 */

#ifndef VERTICE_APPROX_2D_H_
#define VERTICE_APPROX_2D_H_

#include "core/continuous/support_function_provider.h"
#include "math/numeric/basic_geometry.h"
#include "math/type_conversion.h"

namespace continuous {
namespace support_function {

/** This class computes 2D vertex approximations of a convex set S;
 * in the 2D subspace spanned by given orthonormal vectors v and w, and with
 * a maximal error bound eps.
 *
 * This version does not require S to provide support vectors.
 * They are inferred by interpolation.
 *
 * A vertex approximation of a convex set S consists of two sets of points,
 * X and X_{over}, such that chull(X) \subseteq S \subseteq chull(X \cup X_{over}).
 * The error bound signifies that for any x in chull(X) there exists an x' in
 * chull(X \cup X_{over}) such that for all z that are convex combinations of v
 * and w, |z^T (x'-x)| < eps.
 * @todo Check error bound.
 */

template<typename scalar_type>
class vertice_approx_2D {
public:
	struct sample {
		scalar_type l;
		scalar_type a;
		scalar_type b;
		bool is_over;
	};
	typedef std::list<sample> sample_list;
	typedef math::vdom_vector<scalar_type> direction;

	/** Can approximate the support function of S in directions
	 * that are linear combinations of v and w up to accuracy eps>=0.
	 * Terminates if eps>0 or S is a bounded polyhedron.
	 */
	vertice_approx_2D(direction v, direction w,
			const support_function_provider& S, scalar_type eps) :
		my_v(v), my_w(w), my_S(S), my_eps(eps), my_approx_eps(scalar_type(1)
				/ scalar_type(100000)) {
		assert(math::numeric::approx_comparator<scalar_type>::is_maybe_zero(math::scalar_product(v,w)));
	}
	;

	/** Get the sf and infer the sv by interpolating around lopt. */
	sample sample(const scalar_type& x1) {
		bool is_empty, is_bounded;

		math::vdom_vector<scalar_type> l1 = my_v + x1 * my_w;
		//std::cout << "opt " << my_v << "+"<< lopt << "*"<< my_w << "=" << le;
		math::vdom_vector<scalar_type> sv;
		scalar_type y1;
		my_S.compute_support(l1, y1, sv, is_empty, is_bounded);
		if (!is_bounded)
			throw std::runtime_error("sv_2D_approximator: unbounded set\n");
		if (is_empty)
			throw std::runtime_error("sv_2D_approximator: empty set\n");

		bool found = false;
		scalar_type x2, x3, y2, y3;
		do {
			x2 = x1 - my_approx_eps;
			math::vdom_vector<scalar_type> l2 = my_v + x2 * my_w;
			my_S.compute_support(l2, y2, sv, is_empty, is_bounded);
			x3 = x1 + my_approx_eps;
			math::vdom_vector<scalar_type> l3 = my_v + x3 * my_w;
			my_S.compute_support(l3, y3, sv, is_empty, is_bounded);
			found = math::numeric::is_colinear(x2, y2, x1, y1, x3, y3);
			if (!found) {
				my_approx_eps = my_approx_eps / scalar_type(2);
			}
		} while (!found);

		scalar_type b = (y3 - y2) / (x3 - x2);
		scalar_type a = y3 - (b * x2);

		std::cout << " at " << sv << "(" << a << "," << b << ")" << " eps " << my_approx_eps << std::endl;
		// insert sample after i == before j
		typename vertice_approx_2D<scalar_type>::sample new_s = { x1, a, b,
				false };
		return new_s;
	}
	;

	void approx(sample_list& pts, typename sample_list::iterator& i,
			const typename sample_list::const_iterator& iend) {
		typename sample_list::iterator j = i;
		++j;

		scalar_type lopt, sn, sp, ap, bp, errfact(1);
		scalar_type a1, a2, b1, b2, l1, l2;
		while (j != iend) {
			a1 = i->a;
			b1 = i->b, l1 = i->l;
			a2 = j->a;
			b2 = j->b, l2 = j->l;
			//std::cout << "checking l1=" << l1 << ", b1=" << b1  << ", l2=" << l2  << ", b2=" << b2 << "; " << std::endl;
			bool refine = true;
			bool within_err = false;
			refine
					= math::numeric::approx_comparator<scalar_type>::is_definitely_strictly_larger(
							b2, b1);
			if (refine) { // no. 1
				// optimal lambda
				lopt = (a1 - a2) / (b2 - b1);
				// lopt needs to be between l1 and l2
				refine
						= math::numeric::approx_comparator<scalar_type>::is_definitely_strictly_larger(
								lopt, l1) && math::numeric::approx_comparator<
								scalar_type>::is_definitely_strictly_larger(l2,
								lopt);
							//std::cout << "b2>b1, l1<" << lopt << "<l2:" << refine << "; " << std::endl;
			}
			if (refine) { // no. 2
				// lower bound
				sn = a1 + b1 * lopt;
				// upper bound
				bp = (a2 - a1 + b2 * l2 - b1 * l1) / (l2 - l1);
				ap = a1 + (b1 - bp) * l1;
				sp = ap + bp * lopt;
				// check whether in bounds
				// correct so that angles are evenly distributed
				double lopt_double = convert_element<double> (lopt);
				// The factor 2 below comes from the coordinate transform
				errfact = scalar_type(2.0 * std::sqrt((1.0 + lopt_double
						* lopt_double)));
				refine
						= math::numeric::approx_comparator<scalar_type>::is_definitely_strictly_larger(
								sp - sn, my_eps * errfact);
				// if 0 < err < eps we don't refine, but we must insert another vertice
				// (otherwise we get an underapproximation)
				within_err = !refine && math::numeric::approx_comparator<
						scalar_type>::is_definitely_strictly_larger(sp, sn);
				// std::cout << "err:" << sp - sn << "; ";
			}
			// either obtain a new sample,
			// insert an overapproximation vertice,
			// or skip to the next point
			if (refine) {
				// obtain sample at lopt
				typename vertice_approx_2D<scalar_type>::sample new_s =
						sample(lopt);
				typename sample_list::iterator tmp_j = j;
				//				j = pts.insert(j, new_s);
				// if new_s.b == b2 then l2 is obsolete; however, for now don't delete
				// since otherwise we have to convert the lambdas below properly

				// has a bug
				if (math::numeric::approx_comparator<scalar_type>::is_maybe_equal(
						b2, new_s.b) || math::numeric::approx_comparator<
						scalar_type>::is_maybe_equal(b1, new_s.b)) {
					i = j;
					++j;
				} else {
					j = pts.insert(j, new_s);
					//std::cout << "sample:" << new_s.a << "," << new_s.b << std::endl;
				}

				//				std::cout << "inserting:" << pts.size() << " pts, err " << sp-sn << std::endl;
			} else if (within_err) {
				//std::cout << "over:" << ap << "," << bp << std::endl;
				typename vertice_approx_2D<scalar_type>::sample new_s = {
						lopt, ap, bp, true }; // it's an overapproximation point, so "true"
				i = pts.insert(j, new_s);
			} else {
				// within error bounds, move on to the next one
				i = j;
				++j;
				//				std::cout << "advancing." << std::endl;
			}
		}
	}
	;

	/** Sample in variables x_id and y_id up to accuracy eps.
	 * The returned points are sorted ccwise.
	 *
	 * @Note: There are redundant points if the support vector in
	 * two axis directions is the same. This could be filtered out,
	 * but then again it's only max. 4 points. */
	static sample_list approx(const support_function_provider& S,
			variable_id x_id, variable_id y_id, scalar_type eps) {
		sample_list pts;

		// Prepare variable stuff
		variable_id_set vis;
		vis.insert(x_id);
		vis.insert(y_id);
		index_to_variable_id_map_ptr iimap =
				index_to_variable_id_map::get_map_with_ids(vis);

		// Initialize in the axis directions
		direction x(iimap);
		direction y(iimap);
		//		// 45 degree directions
		//		x.set_existing_coeff_with_id(x_id, scalar_type(1));
		//		y.set_existing_coeff_with_id(y_id, scalar_type(1));
		x.set_existing_coeff_with_id(x_id, scalar_type(1));
		x.set_existing_coeff_with_id(y_id, scalar_type(1));
		y.set_existing_coeff_with_id(x_id, scalar_type(-1));
		y.set_existing_coeff_with_id(y_id, scalar_type(1));

		typename sample_list::iterator i;

		// Sample 4 sectors in ccw order with lambda in [-1,1]
		typename vertice_approx_2D<scalar_type>::sample spl;
		// 1. v=x,w=y
		{
			vertice_approx_2D<scalar_type> A(x, y, S, eps);
			spl = A.sample(scalar_type(-1));
			//std::cout << "seed1:" << spl.a << "," << spl.b << std::endl;
			pts.push_back(spl);
			spl = A.sample(scalar_type(1));
			//std::cout << "seed2:" << spl.a << "," << spl.b << std::endl;
			pts.push_back(spl);
			i = pts.begin();
			A.approx(pts, i, pts.end());
		}
		{
			//			std::cout << "2." << std::endl;
			// 2. v=y,w=-x
			vertice_approx_2D<scalar_type> A(y, -x, S, eps);
			// reuse the last sample -> set lambda to -1
			sample_list pts2;
			spl = *pts.rbegin();
			typename vertice_approx_2D<scalar_type>::sample spl2 = {
					scalar_type(-1) / spl.l, spl.b, -spl.a };
			pts2.push_back(spl2);
			// get the second point
			spl = A.sample(scalar_type(1));
			pts2.push_back(spl);
			i = pts2.begin();
			A.approx(pts2, i, pts2.end());

			// copy the corrected values back to the original list
			// start with the second entry since the first was copied from the
			// previous round
			for (typename sample_list::const_iterator it = ++pts2.begin(); it
					!= pts2.end(); ++it) {
				spl = *it;
				typename vertice_approx_2D<scalar_type>::sample spl2 = {
						spl.l, -spl.b, spl.a };
				pts.push_back(spl2);
			}
		}
		{
			//			std::cout << "3." << std::endl;
			// 3. v=-x,w=-y
			vertice_approx_2D<scalar_type> A(-x, -y, S, eps);
			// reuse the last sample -> set lambda to -1
			sample_list pts2;
			spl = *pts.rbegin();
			typename vertice_approx_2D<scalar_type>::sample spl2 = {
					scalar_type(-1) / spl.l, -spl.a, -spl.b };
			pts2.push_back(spl2);
			// get the second point
			spl = A.sample(scalar_type(1));
			pts2.push_back(spl);
			i = pts2.begin();
			A.approx(pts2, i, pts2.end());

			// copy the corrected values back to the original list
			for (typename sample_list::const_iterator it = ++pts2.begin(); it
					!= pts2.end(); ++it) {
				spl = *it;
				typename vertice_approx_2D<scalar_type>::sample spl2 = {
						spl.l, -spl.a, -spl.b };
				pts.push_back(spl2);
			}
		}
		{
			//			std::cout << "4." << std::endl;
			// 4. v=-y,w=x
			vertice_approx_2D<scalar_type> A(-y, x, S, eps);
			// reuse the last sample -> set lambda to -1
			sample_list pts2;
			spl = *pts.rbegin();
			typename vertice_approx_2D<scalar_type>::sample spl2 = {
					scalar_type(-1) / spl.l, -spl.b, spl.a };
			pts2.push_back(spl2);
			// get the second point
			// we're full way around the circle, so
			// instead of spl = A.sample(scalar_type(1));
			// we can just use the first one
			spl = *pts.begin();
			typename vertice_approx_2D<scalar_type>::sample spl3 = {
					scalar_type(-1) / spl.l, -spl.b, spl.a };
			pts2.push_back(spl3);
			i = pts2.begin();
			A.approx(pts2, i, pts2.end());

			// copy the corrected values back to the original list
			// don't copy the last one, since it was taken from the first round
			for (typename sample_list::const_iterator it = ++pts2.begin(); &(*it)
					!= &(*pts2.rbegin()); ++it) {
				spl = *it;
				typename vertice_approx_2D<scalar_type>::sample spl2 = {
						spl.l, spl.b, -spl.a };
				pts.push_back(spl2);
			}
		}

		// correct for the initial coordinate transform
		// copy the corrected values back to the original list
		sample_list pts2;
		for (typename sample_list::const_iterator it = pts.begin(); it
				!= pts.end(); ++it) {
			spl = *it;
			typename vertice_approx_2D<scalar_type>::sample spl2 = { spl.l,
					(spl.a - spl.b) / scalar_type(2), (spl.a + spl.b)
							/ scalar_type(2) };

			pts2.push_back(spl2);
		}
		return pts2;
	}
	;

private:
	direction my_v;
	direction my_w;
	const support_function_provider& my_S;
	scalar_type my_eps;
	scalar_type my_approx_eps;
};

}
}

#endif /* VERTICE_APPROX_2D_H_ */
