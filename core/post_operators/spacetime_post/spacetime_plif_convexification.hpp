/*
 * spacetime_plif_convexification.hpp
 *
 *  Created on: Oct 22, 2012
 *      Author: notroot
 */

#ifndef SPACETIME_PLIF_CONVEXIFICATION_HPP_
#define SPACETIME_PLIF_CONVEXIFICATION_HPP_

#include "spacetime_plif.h"

namespace spacetime {

/** Measure the approximation error over a given set of evolutions and scopes */
class measure_concave_approx_error {
public:
	typedef std::pair<spacetime_plif::evolution_cache::const_iterator, plif::interval>
			evolution_with_scope;
	typedef std::list<evolution_with_scope> evo_scopes;
	typedef spacetime_plif::scalar_type scalar_type;
	typedef spacetime_plif::plf_type plf_type;
	typedef std::list<plf_type> plf_list;
	typedef spacetime_plif::time_interval time_interval;

	measure_concave_approx_error(time_interval domain=time_interval::whole()) {
		my_domain = domain;
	};

	measure_concave_approx_error(const evo_scopes& es,time_interval domain=time_interval::whole()) {
		my_domain = domain;
		// construct the list of plfs to measure by restricting them to their scopes
		for (evo_scopes::const_iterator it=es.begin();it!=es.end();++it) {
			const plf_type& f = it->first->second.first.get_upper();
//std::cout << "restricting measure to [" << it->second.lower() << "," << it->second.upper() << "]" << std::endl;
			time_interval dom = plif::intersect(my_domain,it->second);
			plf_type scoped_f = plif::restrict(f,dom);
//scoped_f.display(std::cout);
			// we need at least 3 points for the convex hull to make a difference
			// @todo test also convexity here
			if (scoped_f.get_list().size()>2) {
				my_plfs.push_back(scoped_f);
			}
		}
	};

	/** Returns the approximation error for the given inflection point position */
	scalar_type operator()(scalar_type x) {
//std::cout << "chull_err(" << x << "):";
		scalar_type max_err(0);
		// For each of the evolutions
		std::vector<spacetime_plif::duration> cut_points;
		cut_points.push_back(x);
		for (plf_list::const_iterator it = my_plfs.begin(); it != my_plfs.end(); ++it) {
			// test if the inflection point is inside the interval;
			// we only need to measure if it's strictly inside
			const spacetime_plif::time_interval& itv = it->get_domain();
			// restrict x to the domain
			scalar_type x_here =
					std::max(itv.lower(), std::min(x, itv.upper()));
			// dissect at the inflection point
			std::vector<plf_type> split_f;
			if (x>itv.lower() && x<itv.upper()) {
				split_f = dissect(*it, cut_points);
			} else {
				split_f.push_back(*it);
			}
			for (std::vector<plf_type>::const_iterator jt = split_f.begin(); jt
					!= split_f.end(); ++jt) {
//				jt->display(std::cout);
				if (jt->get_list().size()>1) {
					plf_type ch_of_jt = plif::concave_hull(*jt);
//					plif::piecewise_linear_interval_function margin(*jt, ch_of_jt);
//					plif::plif_graph(
//							margin,
//							"X",
//							"/tmp/test_spacetime_ch");

					scalar_type err = -plif::infimum(*jt - ch_of_jt);
//					std::cout << err << ",";
					max_err = std::max(max_err, err);
				}
			}

		}
		return max_err;
	};

private:
	plf_list my_plfs;
	time_interval my_domain;
};

inline
spacetime_plif::time_point_set spacetime_plif::compute_cut_points(
		const cut_point_method& m) const {
	LOGGERSWOC(DEBUG,"compute_cut_points","computing cut points");

	time_point_set cut_points;

	// choose cutpoints according to given method:
	// FIXED_PIECE_COUNT, TIME_STEP, MIN_CONCAVE_PIECES

	double t_first = get_time_domain().lower();
	double t_last = get_time_domain().upper();
	double T = t_last - t_first;
	// Cut into a single piece: no cutpoints
	if (m.type == cut_point_method::FIXED_COUNT_SAME_SIZE_PIECES) {
		size_t N = m.piece_count;
		if (N < 1) {
			throw basic_exception("can't produce less than 1 piece");
		}
		// cut points at regular intervals
		for (size_t i = 1; i < N; ++i) {
			cut_points.insert(t_first + T / N * i);
		}
	} else if (m.type == cut_point_method::TIME_STEP) {
		size_t N = ceil(T / m.time_step);
		for (size_t i = 1; i < N; ++i) {
			cut_points.insert(t_first + m.time_step * i);
		}
	} else if (m.type == cut_point_method::MIN_CONCAVE_PIECES) {
		// do sophisitacted stuff
		bool optim_cut_points = false;

		// for each direction, get the cut point intervals
		std::vector<plif::interval> all_intervals;
		typedef std::pair<evolution_cache::const_iterator,plif::interval> evolution_with_scope;
		std::vector<evolution_with_scope> evo_with_neighb_intv;

		for (evolution_cache::const_iterator it = my_evolutions.begin(); it
				!= my_evolutions.end(); ++it) {
			plf_type f = it->second.first.get_lower();
			const plf_type& g = it->second.first.get_upper();

//			std::cout << std::endl << "lower:" << f;
//			std::cout << std::endl << "upper:" << g;

			// can't approximate with less than the actual error
			// @todo do a proper measurement here
			// error_type conv_err = std::max(m.approx_error,it->second.second);
			error_type conv_err = m.approx_error;

			if (!m.lower_is_error_reference) {
				// instead of the lower bound compute the error from the upper bound
				//f = g;

				// we want to add to the approx error
				conv_err = it->second.second + m.approx_error;
			}
			plf_type f_shifted;
			f_shifted = augment_by_error(f,conv_err);

			// to make sure things are not infeasible, take the max with the upper bound
			f_shifted = pointwise_maximum(f_shifted,g);

			std::pair<std::vector<plif::interval>,std::vector<plif::interval> > infl_info = get_inflection_intervals(g,f_shifted,my_numeric_error.rel(),my_numeric_error.abs());
			std::vector<plif::interval>& ivec = infl_info.first;
			// the inflection intervals of the current evolution to the rest
			all_intervals.reserve(all_intervals.size()+ivec.size());
//			all_intervals.insert(all_intervals.end(),ivec.begin(),ivec.end());
			// trying to make sure things stay in the same order
			for (std::vector<plif::interval>::const_iterator nit =
					infl_info.first.begin(); nit != infl_info.first.end(); ++nit) {
				all_intervals.push_back(*nit);
			}
			// remember for each interval the evolution and the scope interval
			if (optim_cut_points) {
				evo_with_neighb_intv.reserve(
						evo_with_neighb_intv.size() + ivec.size());
				for (std::vector<plif::interval>::const_iterator nit =
						infl_info.second.begin(); nit != infl_info.second.end();
						++nit) {
					evo_with_neighb_intv.push_back(
							evolution_with_scope(it, *nit));
				}
			}
		}
		using plif::operator<<;
//		std::cout << "all intervals: " << all_intervals << std::endl;
		std::pair<std::vector<plif::interval>,std::vector<std::vector<size_t> > > overlap_info =
		 plif::partition(all_intervals);
		std::vector<plif::interval>& overlap_intervals = overlap_info.first;
		std::vector<std::vector<size_t> >& group_indices = overlap_info.second;

//		std::cout << "overlapping intervals: " << overlap_intervals << std::endl;
		//std::cout << "group indices: " << group_indices << std::endl;

		// put all the

		std::vector<plif::interval> scope_intervals(overlap_intervals.size(),plif::interval::empty());
		// choose cut_points in overlap intervals
		size_t i = 0;
		for (std::vector<plif::interval>::const_iterator it =
				overlap_intervals.begin(); it != overlap_intervals.end(); ++it) {
//std::cout << "overlapping interval: [" << it->lower() << "," << it->upper() << "]" << std::endl;
			// The involved evolutions and scope intervals
			const std::vector<size_t>& group = group_indices[i];

			scalar_type x_low = it->lower();
			scalar_type x_upper = it->upper();
			scalar_type x_middle = (x_low + x_upper) / scalar_type(2);
			scalar_type x_min;

			if (optim_cut_points) {
				time_interval measure_intv(get_time_domain());
				/** @attention This is assuming cut_points go from left to right,
				 * so they always restrict the lower bound.
				 */
//			if (!cut_points.empty())
//				measure_intv = time_interval(*cut_points.rbegin(),measure_intv.upper());
				//std::cout << "measure interval: [" << measure_intv.lower() << "," << measure_intv.upper() << "]" << std::endl;
				measure_concave_approx_error measure(measure_intv);
				if (true) {
					// Get the list of evolutions and scopes involved in this group
					//std::cout << "Group:";
					std::list<evolution_with_scope> evo_scopes;
					for (size_t j = 0; j < group.size(); ++j) {
						size_t all_intervals_index = group[j];
						//std::cout << "group interval: [" << all_intervals[all_intervals_index].lower() << "," << all_intervals[all_intervals_index].upper() << "] from ";
						evo_scopes.push_back(
								evo_with_neighb_intv[all_intervals_index]);
						//std::cout << evo_with_neighb_intv[all_intervals_index].first->first << std::endl;
						// compute scope of group
						scope_intervals[i] =
								plif::hull(scope_intervals[i],
										evo_with_neighb_intv[all_intervals_index].second);
					}
					//std::cout << "scope interval: [" << scope_intervals[i].lower() << "," << scope_intervals[i].upper() << "]" << std::endl;
					// Optimize the cut_point in the interval *it by measuring the error in evo_scopes
					measure = measure_concave_approx_error(evo_scopes,
							measure_intv);
				}

				// measure the error at lower bound
				scalar_type err_low = measure(x_low);
				scalar_type err_upper = measure(x_upper);
				scalar_type err_mid = measure(x_middle);

				std::cout << "at lower : " << err_low;
				std::cout << "at upper : " << err_upper;
				std::cout << "at middle: " << err_mid << std::endl;

				/*			size_t K=10;
				 for (size_t k=1;k<K;++k) {
				 scalar_type x = x_low + (x_upper-x_low) / scalar_type(K) * k;
				 scalar_type y=measure(x);
				 if (y<err_mid) {
				 err_mid = y;
				 x_middle = x;
				 }
				 }
				 */

				if (err_low < err_upper)
					if (err_low < err_mid)
						x_min = x_low;
					else
						x_min = x_middle;
				else if (err_upper < err_mid)
					x_min = x_upper;
				else
					x_min = x_middle;

				std::cout << "at search : " << err_mid << std::endl;
			} else {
				x_min = x_middle;
			}

			cut_points.insert(x_min);
			++i;
		}
	} else if (m.type == cut_point_method::ALL_PIECES) {
		// for each direction, get the cut points, and take their union
		for (evolution_cache::const_iterator it = my_evolutions.begin(); it
				!= my_evolutions.end(); ++it) {
			const plf_type& g = it->second.first.get_upper();

			// get the inflection points
			std::pair<std::vector<scalar_type>, std::vector<scalar_type> > a_infl_info =
					plif::inflection_points(g);

			// add all inflection points to the cut points
			LOGGER(DEBUG7, __FUNCTION__,
									"adding "+to_string(a_infl_info.first.size())+" cutpoints");
			cut_points.insert(a_infl_info.first.begin(),a_infl_info.first.end());
		}
	} else {
		throw basic_exception("unknown cut_point_method type");
	}
	LOGGER(DEBUG5, __FUNCTION__,
							"total "+to_string(cut_points.size())+" cutpoints");
	return cut_points;
}

inline std::vector<spacetime_plif> spacetime_plif::convexify(
		const cut_point_method& m) const {
	LOGGERSWOC(DEBUG, __FUNCTION__, "convexifying flowpipe");

	std::vector<spacetime_plif> new_splifs;

	// Find a inflection points for cutting into convex parts
	time_point_set cut_points = compute_cut_points(m);

	size_t M = cut_points.size() + 1; // number of pieces
	size_t N = size(); // number of evolutions

	// no evolutions, so nothing to do
	if (N == 0) {
		// return a universe set
		new_splifs.push_back(*this);
		return new_splifs;
	}

	// convert the cutpoints to intervals
	std::vector<duration> cut_point_vec(cut_points.begin(),cut_points.end());
	std::vector<time_interval> subdomain_vector = plif::cut_points_to_intervals(cut_point_vec,get_time_domain());
	M = subdomain_vector.size();

	// make virgin splifs, each for the corresponding time interval
	new_splifs = std::vector<spacetime_plif>(M,*this);

	// for each piece, restrict to the corresponding subdomain, then take convex hull
	for (size_t i = 0; i < M; ++i) {
		new_splifs[i].restrict_to_subdomain(subdomain_vector[i]);
		new_splifs[i].assign_convex_hull();
	}

	return new_splifs;
}

inline polyhedron_collection<spacetime_plif::scalar_type> spacetime_plif::compute_outer_polyhedra(
		const cut_point_method& m) const {
	LOGGERSWOC(DEBUG, "compute_outer_polyhedra", "computing outer polyhedra");

	// Find a inflection points for cutting into convex parts
	time_point_set cut_points = compute_cut_points(m);

	size_t M = cut_points.size() + 1; // number of pieces
	size_t N = size(); // number of evolutions

	polyhedron_collection<scalar_type> coll;
	{ // scope for the stopwatch
//		LOGGERSW(DEBUG3, "compute_outer_polyhedra", "constructing polyhedra");
		logger::logger_id logid=LOGGER(DEBUG6, "compute_outer_polyhedra",
				"computing "+to_string(M)+" polytopes for "+to_string(N)+ " directions");

		// Prepare a vector of M universe polyhedra
		// We must by convention return at least one polyhedron
		if (N == 0) {
			// return a universe set
			polyhedron<scalar_type>::ptr poly = polyhedron<scalar_type>::ptr(
					new constr_polyhedron<scalar_type>());
			return polyhedron_collection<scalar_type>(poly);
		}

		std::vector<constr_polyhedron<scalar_type>::ptr> polys(M);
		// assign to each pointer a universe polyhedron
		// (can't do it with vector constructor since it's got to be a different pointer each time)
		for (size_t i = 0; i < M; ++i) {
			polys[i] = constr_polyhedron<scalar_type>::ptr(
					new constr_polyhedron<scalar_type>());
		}

		// For each of the N evolutions, create a concave plf for each of the M pieces,
		// then add it to the corresponding polyhedron
		spacetime_plif* nonconst_this = const_cast<spacetime_plif*>(this);
		for (evolution_cache::iterator it =
				nonconst_this->my_evolutions.begin(); it != nonconst_this->my_evolutions.end(); ++it) {
			const direction& dir = it->first;

			// dissect the evolution according to the cut points
			std::vector<duration> cut_point_vector(cut_points.begin(),
					cut_points.end());
			std::vector<plif_type> plif_pieces = dissect(it->second.first,
					cut_point_vector);

			if (plif_pieces.size() != M)
				throw basic_exception(
						"number of dissected pieces doesn't correspond to cut points");

			// for each piece, compute it's concave hull, then obtain its constraints and add them
			// to the corresponding polyhedron
//std::cout << "cut points: " << cut_point_vector << std::endl;
//std::cout << "direction: " << dir << std::endl;
			for (size_t i = 0; i < M; ++i) {

				plf_type cchull = simplify(concave_hull(plif_pieces[i].get_upper()));
				error_type conv_err = m.approx_error;
				/*
				if  (conv_err.abs()<it->second.second.abs()) {
					conv_err = error_type(conv_err.rel(),it->second.second.abs());
				}
				if  (conv_err.rel()<it->second.second.rel()) {
					conv_err = error_type(it->second.second.rel(),conv_err.abs());
				}
				*/

				if (my_simplify_convex && plif_pieces[i].get_upper().size()>4 && m.approx_error.abs()>1e-12) {
//				if (m.type==cut_point_method::FIXED_COUNT_SAME_SIZE_PIECES && M==1 && my_simplify_convex && plif_pieces[i].get_upper().size()>4 && m.approx_error.abs()>1e-12) {
//				if (my_simplify_concave && cchull.size()>4) {					// simplify and reassign
					plf_type f = plif_pieces[i].get_lower();
					const plf_type& g = cchull;
					size_t before_size = cchull.size();

					/** Relative error <= r <=> g <= max(f/(1+r),f*(1+r)) */
					plf_type f_shifted;
					if (!m.lower_is_error_reference) {
						// instead of the lower bound compute the error from the upper bound
						f = cchull;
					}
					if (conv_err.rel()>1e-12) {
						scalar_type a = conv_err.rel()+scalar_type(1);
						plf_type f_plus_rel_err = pointwise_maximum(a*f,(scalar_type(1)/a)*f);
						f_shifted = pointwise_maximum(f_plus_rel_err,f+conv_err.abs())+conv_err.abs();
					} else {
						f_shifted = f+conv_err.abs();
					}
					// to make sure things are not infeasible, take the max with the upper bound
					f_shifted = simplify(pointwise_maximum(f_shifted,g));

					//plf_type f_shifted = g+conv_err.abs();
					//plf_type g_approx = plif::min_link_greedy_left(g,f_shifted);
					//plf_type g_approx = shortest_path(g,f_shifted,g.get_list().front(),g.get_list().back(),true);
					cchull = simplify(min_link_concave_relax(cchull,f_shifted));
//					plf_type g_approx = plif::pw_convex_approx_greedy_left(g,f_shifted);

//					plif::plf_graph(g_approx,"gif","/tmp/test_greedy");
//					plif::piecewise_linear_interval_function margin(g,f_shifted);
//					plif::plif_graph(margin,"X","/tmp/test_spacetime_margin_greedy","-m3 -W 0.006 -S 4 /tmp/test_greedy.txt -s -W 0 ");

//					cchull = concave_hull(g_approx);

					// put it in the evolution so future containment can be made
					plif_pieces[i].set_upper(cchull);
					IFLOGGER(DEBUG6) {
						LOGGER(DEBUG6, "compute_outer_polyhedra",
								"with err "+to_string(conv_err) + " reduced "+to_string(before_size)+" breakpoints to "+to_string(cchull.size()));
					}
				} else {

//					cchull = plif_pieces[i].get_upper();

//					plif::piecewise_linear_interval_function margin(plif_pieces[i].get_lower(),plif_pieces[i].get_upper());
//					plif::plif_graph(margin,"X","/tmp/test_spacetime_margin_greedy","-S 4 -W 0 ");
				}

//cchull.display(std::cout);
				lin_constraint_system<scalar_type> lin_cons;
				{
//					LOGGERSW(DEBUG5, "compute_outer_polyhedra","computing convex hull");
					lin_cons = plf_to_constraints(it->first, cchull);
				}
//				if (lin_cons.size()==0) {
//				std::cout << "upper: "; plif_pieces[i].get_upper().display(std::cout);
//				std::cout << "chull: "; cchull.display(std::cout);
//				}

				polys[i]->add_constraints(lin_cons);
			}

		}

		// insert resulting polys in collection
		// get the new domain
		positional_vdomain dom_t = my_aff.dyn.domain();
		dom_t.add_variable(my_time_variable);
		positional_vdomain::size_type t_pos = dom_t.pos(my_time_variable);

		// add constraints to limit time
		direction time_dir(dom_t);
		duration t_beg = get_time_domain().lower();
		duration t_end = get_time_domain().lower();

		time_point_set::const_iterator cut_it = cut_points.begin();
		unsigned int nb_cons_max = 0;
		for (size_t i = 0; i < M; ++i) {
			// add constraints restricting time
			lin_constraint_system<scalar_type> lin_cons;
			t_beg = t_end;
			if (cut_it != cut_points.end()) {
				t_end = *cut_it;
				++cut_it;
			} else
				t_end = get_time_domain().upper();

			// first: t <= t_end <=> t - t_end <= 0
			time_dir[t_pos] = scalar_type(1);
			lin_constraint<scalar_type> con(time_dir, -t_end, LE);
			lin_cons.insert(con);
			// second: t >= t_beg <=> -t + t_beg <= 0
			time_dir[t_pos] = scalar_type(-1);
			con = lin_constraint<scalar_type>(time_dir, t_beg, LE);
			lin_cons.insert(con);
//				std::cout << std::endl << "time constraints: " << lin_cons << std::endl;
			polys[i]->add_constraints(lin_cons);

			// for stats, compute the nb of constraints
			nb_cons_max = std::max(nb_cons_max,polys[i]->get_constraints()->size());

			coll.insert(polys[i]);
			IFLOGGER(DEBUG7) {
				LOGGER(DEBUG7, "compute_outer_polyhedra",
						"got outer polyhedron with "+to_string(polys[i]->get_constraints()->size())+" constraints");
			}
		}
		LOGGER_ATTACH(DEBUG5, "compute_outer_polyhedra",
					" got up to "+to_string(nb_cons_max)+" constraints", logid);
	}
	return coll;
}

}

#endif /* SPACETIME_PLIF_CONVEXIFICATION_HPP_ */
