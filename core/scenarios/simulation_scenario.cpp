#include "core/scenarios/simulation_scenario.h"
#include "math/ode_solving/traj_simu/postc_simulation.h"
#include "math/ode_solving/traj_simu/postd_simulation.h"
#include "core/scenarios/support_fun_scenario.h"
#include "core/hybrid_automata/adaptors/support_function_adaptors.h"

#include "utility/stl_helper_functions.h"

namespace hybrid_automata{


static const global_types::coefficient_type simulation_number_type =
		global_types::STD_DOUBLE;

namespace {
template<typename T>
struct rel_error_setter {
	static bool implement(double c) {
		math::numeric::approx_comparator<T>::set_rel_error_bound(
				convert_element<T> (c));
		return true;
	}
	;
};
template<typename T>
struct abs_error_setter {
	static bool implement(double c) {
		math::numeric::approx_comparator<T>::set_abs_error_bound(
				convert_element<T> (c));
		return true;
	}
	;
};
}

void simulation_apply_options(reachability_scenario* sp,
		options::options_processor::variables_map& vmap) {
	/** @todo This should somehow not be necessary at this point. Global policy? */

	std::string s;
	double gth = -1;
	unsigned int algo = 0;
	unsigned int unisampl = 0;

	if (options::options_processor::get_string_option(vmap, "time-horizon", s)) {
		sp->set_time_horizon(from_string<double> (s));
	}
	if (options::options_processor::get_string_option(vmap, "global-time-horizon", s)) {
		gth=from_string<double> (s);
	}
	if (options::options_processor::get_string_option(vmap, "sampling-time", s)) {
		sp->set_sampling_time(from_string<double> (s));
	}
	if (options::options_processor::get_string_option(vmap, "simu-intern-algo", s)) {
		algo = from_string<double> (s);
	}
	if (options::options_processor::get_string_option(vmap, "simu-init-sampling-points", s)) {
			unisampl = from_string<double> (s);
		}

	double lth = sp->get_time_horizon();
	double sampling = sp->get_sampling_time();

	postc_simulation::ode_parameters par = postc_simulation::ode_parameters();
	par.max_timestep = sampling;
	//TODO set tolerance ??
	postc_simulation::ptr cpost = postc_simulation::ptr(new postc_simulation() );
	sp->set_continuous_post_operator(cpost);

	// choose the discrete_post accordingly
	postd_simulation::ptr  dpost = postd_simulation::ptr(new postd_simulation(cpost));

	cpost->set_time_horizon(lth);
	cpost->set_global_time_horizon(gth);
	cpost->set_ode_parameters(par);
	cpost->set_init_uniform_sampling(unisampl);
	//TODO cpost->set_bounding_box();
	dpost->reset_options();
	dpost->set_simu_internal_algo(algo);

	// set the global tolerance values so they match the ode solver accuracy
	double global_tol_rel = 10*math::ode::ode_defaults::get_rel_tol();
	double global_tol_abs = 10*math::ode::ode_defaults::get_abs_tol();

	global_types::coefficient_type_caller<bool, double,
			rel_error_setter>::call(global_tol_rel, sp->get_number_type());
	global_types::coefficient_type_caller<bool, double,
			abs_error_setter>::call(global_tol_abs, sp->get_number_type());

	sp->set_discrete_post_operator(dpost);

}

reachability_scenario get_simulation_scenario() {

	reachability_scenario s = get_support_fun_scenario();
	s.set_name("simu");

	//	s.set_continuous_post_operator(continuous_post::ptr(
	//			new continuous_post_sfm<global_types::type_selector<
	//					supp_f_number_type>::type> ));
	//	s.set_discrete_post_operator(discrete_post::ptr(
	//			new discrete_post_sfm_inters<global_types::type_selector<
	//					supp_f_number_type>::type> ));

	// let continuous and discrete post be chosen according to default options

	//TODO write new class????
	automaton_to_supp_f_adaptor_ptr aut_ad = automaton_to_supp_f_adaptor_ptr(
			new automaton_to_supp_f_adaptor(global_types::STD_BOOL,
					simulation_number_type,true));
	aut_ad->init();
	s.set_adapt_automaton_visitor(aut_ad);
	s.set_option_handler(&simulation_apply_options);
	return s;
}

}
