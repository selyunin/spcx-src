/*
 * INTV_formatter.cpp
 *
 */

#include "INTV_formatter.h"

#include "core/discrete/discrete_set.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/continuous/support_function/sf_derived/sf_chull.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"
#include "core/hybrid_automata/automaton_name_formatter.h"
#include "io/common_input/input_options.h"
#include "math/vdom/context_variable_formatter.h"
#include "core/continuous/polyhedra/hyperbox/finite_hyperbox_utility.h"
#include "core/continuous/support_function/spacetime_flowpipe.h"
#include "math/ode_solving/traj_simu/continuous_set_simulation.h"
// for computing bounding box
#include "core/continuous/support_function/spacetime_flowpipe_utility.h"

namespace io {

void INTV_formatter::output(const hybrid_automata::symbolic_state& sstate) {
	if (sstate.get_discrete_set()->is_empty()) {
		// generate a warning that the set is empty
		basic_warning("INTV format output",
				"cannot represent empty set in INTV format",
				basic_warning::AMBIGUOUS_OUTPUT);
	} else {
		base_class::output(sstate);
	}
}

continuous::finite_hyperbox<double> get_bounding_box(
		const continuous::spacetime_flowpipe<global_types::float_type>& flowp) {
	using namespace continuous;

	typedef global_types::float_type scalar_type;
	typedef spacetime_flowpipe<scalar_type> flowpipe;

	const positional_vdomain& dom = flowp.domain();
	flowpipe::direction d(dom);
	size_t N = dom.size();

	if (flowp.get_time_domain().is_empty()) {
		return finite_hyperbox<double>::empty_box(dom);
	}

	finite_hyperbox<scalar_type>::vdom_vector_type l(dom),u(dom),c,g;
	if (true) {
		// compute bounding box directions
		// This is just to make sure; if already present no computational overhead.

		// fix const cast
		flowpipe& nonconst_flowp = const_cast<flowpipe&>(flowp);

		for (size_t i = 0; i < N; ++i) {
			// don't compute flowpipes, just obtain what's already computed
			flowpipe::error_type err(1e10, 1e10);
			LOGGER(DEBUG5, "get_bounding_box",
					"computing evolution of variable "+dom.get_variable(i).get_name()+" for bounding box");
			d = flowpipe::direction(dom);
			d[i] = 1.0;
			const flowpipe::spacetime_plif::annotated_plif& evo_ann_pos =
					nonconst_flowp.get_or_compute_evolution(d, err);
			u[i] = plif::supremum(evo_ann_pos.first.get_upper());

//				plif_graph(evo_ann_pos.first, "X", "/tmp/test_gen");
			d[i] = -1.0;
			const flowpipe::spacetime_plif::annotated_plif& evo_ann_neg =
					nonconst_flowp.get_or_compute_evolution(d, err);
//				plif_graph(-(evo_ann_neg.first), "X", "/tmp/test_gen");
			l[i] = -plif::supremum(evo_ann_neg.first.get_upper());
		}
	}

	c = (l+u)*scalar_type(0.5);
	g = u-c;
	finite_hyperbox<scalar_type> bbox(c,g);

	return bbox;
}

/** Obtain a support function provider, doing some conversions if necessary.
 *
 * For sfm_cont_set, the outer poly collection is constructed.
 */
continuous::support_function_provider::ptr convert_to_support_function_provider(
		continuous::continuous_set_ptr c) {

	using namespace continuous;

	support_function_provider::ptr s=support_function_provider::ptr()	;

	support_function::sfm_cont_set<double>::ptr sfm = boost::dynamic_pointer_cast<support_function::sfm_cont_set<double> >(c);
	if (sfm) {
		// convert sfm_cont_set to outer poly collection
		if (!math::definitely(sfm->is_empty())) {
			continuous::polyhedron_collection<double> polys = sfm->get_outer_polytope_collection();
			continuous::support_function::sf_unary_ref<double,continuous::polyhedron_collection<double> > chull(polys);
			s = support_function_provider::ptr(
					new finite_hyperbox<double>(finite_bounding_box<double>(chull)));
//			std::cout << "creating outer polys" << std::endl;
//			std::cout << " with bbox " << finite_bounding_box<double>(*s) << " empty: " << std::boolalpha << s->is_empty() << std::endl;
		}
	} else if (const continuous::spacetime_flowpipe<global_types::float_type> * p =
			dynamic_cast<const continuous::spacetime_flowpipe<
					global_types::float_type> *>(c.get())) {
		s = support_function_provider::ptr(
				new finite_hyperbox<double>(get_bounding_box(*p)));
	} else {
		// otherwise try to obtain a sf provider
		s = boost::dynamic_pointer_cast<support_function_provider>(c);
		if (!s) {
			// otherwise, try if it's a polyhedron_collection
			continuous::polyhedron_collection<double>::ptr polys_ptr =
					boost::dynamic_pointer_cast<
							continuous::polyhedron_collection<double> >(c);
			if (polys_ptr) {
			continuous::support_function::sf_unary_ref<double,
						continuous::polyhedron_collection<double> > chull(
						*polys_ptr);
				s = support_function_provider::ptr(
						new finite_hyperbox<double> (
								finite_bounding_box<double> (chull)));
				//			std::cout << "creating outer polys" << std::endl;
				//			std::cout << " with bbox " << finite_bounding_box<double>(*s) << " empty: " << std::boolalpha << s->is_empty() << std::endl;
			}
		}
	}
	if (!s) {
		throw std::runtime_error("INTV_FORMAT:output: Cannot print cont set\n");
	}
	return s;
}
;

void INTV_formatter::output(const hybrid_automata::symbolic_state_collection& sstates){
	using namespace continuous;

	// set variable output to be context dependent
	// note: deactivated this because streaming the formatter below voids any previous formatter.
	//       when output is written to this stream outside of this function, the formatter no longer
	//       exists and results in a call to an obsolete pointer.
	//context_variable_formatter form(get_context());
	//get_os() << form;

	support_function_provider::ptr global_set_ptr;
	support_function_provider::ptr sf_provider_set_ptr, sf_provider_set_ptr_1;

	/* Create a map of location id to interval */
	typedef std::map<variable_id, bound_interval> var_bounds_map;
	typedef std::map<hybrid_automata::location_constraint_set, var_bounds_map> loc_var_bounds_map;

	typedef std::map<hybrid_automata::location_constraint_set, support_function_provider::ptr> loc_cont_set_ptr_map;
	loc_cont_set_ptr_map loc_set_map = loc_cont_set_ptr_map();

	loc_var_bounds_map loc_bounds_map = loc_var_bounds_map();

	/* Initialize variable-bounds map*/

	bool flag = false;
	discrete::discrete_set_ptr d_set_ptr;
	unsigned int sym_states_cnt = 0;
	for (hybrid_automata::symbolic_state_collection::const_iterator it =
			sstates.begin(); it != sstates.end(); ++it) {
		sym_states_cnt++;
		if (!flag) {
			global_set_ptr = convert_to_support_function_provider((*it)->get_continuous_set());
			if (global_set_ptr && !math::definitely(global_set_ptr->is_empty())) {
				flag = true;
				sf_provider_set_ptr = global_set_ptr;
			}
		} else {
			if (sf_provider_set_ptr) {
				sf_provider_set_ptr_1 = convert_to_support_function_provider((*it)->get_continuous_set());
				if (sf_provider_set_ptr_1 && !math::definitely(sf_provider_set_ptr_1->is_empty())) {
					sf_provider_set_ptr
							= support_function_provider::ptr(
									new support_function::sf_chull<double>(
											sf_provider_set_ptr,
											sf_provider_set_ptr_1));
					global_set_ptr = sf_provider_set_ptr;
				}
			}
		}
		d_set_ptr = (*it)->get_discrete_set();
		std::pair<hybrid_automata::location_constraint_set, var_bounds_map>
				loc_bounds_pair;

		support_function_provider::ptr cont_set_ptr =
				convert_to_support_function_provider(
						(*it)->get_continuous_set());
		support_function_provider::ptr sf_prov_set;

		if (d_set_ptr && cont_set_ptr && !math::definitely(cont_set_ptr->is_empty()))
			for (discrete::discrete_set::const_iterator iter =
					d_set_ptr->begin(); iter != d_set_ptr->end(); ++iter) {
				loc_cont_set_ptr_map::iterator loc_iter = loc_set_map.find(
						*iter); // check *iter

				if (loc_iter != loc_set_map.end() && cont_set_ptr) {
					sf_prov_set = support_function_provider::ptr(
							new support_function::sf_chull<double>(
									loc_iter->second, cont_set_ptr));
					loc_set_map.erase(loc_iter);
					loc_set_map.insert(std::pair<
							hybrid_automata::location_constraint_set,
							support_function_provider::ptr>(*iter, sf_prov_set));
				} else {
					loc_set_map.insert(
							std::pair<hybrid_automata::location_constraint_set,
									support_function_provider::ptr>(*iter,
									cont_set_ptr));
				}
			}
	}
	/* Construct the variable bounds on locations using the sf_chull cont_set created for each location_constraint */

	if (!loc_set_map.empty() && global_set_ptr) {
		typedef double scalar_type;
		scalar_type t_s, t_e;
		bool is_empty_l, is_bounded_l, is_empty_u, is_bounded_u;
		variable_id v_id;
		math::vdom_vector<scalar_type> l_s, l_e, support_vec;
		positional_vdomain set_dom(global_set_ptr->get_variable_ids());
		support_function_provider::ptr cont_set_ptr;
		bound_interval v_bound;

		variable_id_list vars = get_output_variables();
		if (vars.empty()) {
			variable_id_set vars_set = global_set_ptr->get_variable_ids();
			/* Make a variable id list from variable id set */
			for (variable_id_set::const_iterator it = vars_set.begin(); it
					!= vars_set.end(); ++it) {
				vars.push_back(*it);
			}
		}

		for (loc_cont_set_ptr_map::iterator it = loc_set_map.begin(); it
				!= loc_set_map.end(); ++it) {
			cont_set_ptr = it->second;

			var_bounds_map loc_var_bounds = var_bounds_map();

			for (variable_id_list::const_iterator v_it = vars.begin(); v_it
					!= vars.end(); ++v_it) {
				v_id = *v_it;
				l_e = math::vdom_vector<scalar_type>(set_dom);
				l_s = math::vdom_vector<scalar_type>(set_dom);

				l_e.set_coeff_with_id(v_id, scalar_type(1));
				cont_set_ptr->compute_support(l_e, t_e, support_vec,
						is_empty_u, is_bounded_u);
				l_s.set_coeff_with_id(v_id, scalar_type(-1));
				cont_set_ptr->compute_support(l_s, t_s, support_vec,
						is_empty_l, is_bounded_l);

				v_bound = bound_interval();
				if (is_bounded_l)
					v_bound.set_lower(-t_s);
				if (is_bounded_u)
					v_bound.set_upper(t_e);
				if (is_empty_l || is_empty_u)
					v_bound.set_empty();

				loc_var_bounds.insert(std::pair<variable_id, bound_interval>(
						v_id, v_bound));
			}
			loc_bounds_map.insert(std::pair<
					hybrid_automata::location_constraint_set, var_bounds_map>(
					it->first, loc_var_bounds));
		}
		/* Print the Global bounds on the system variables */
		get_os() << "\nBounds on the variables over the entire set:" << std::endl;
		// base_class::output(*global_set_ptr);
		// @todo don't know why the call to base class doesn't work any more
		INTV_continuous_set_formatter<continuous::support_function_provider>::output(*this,*global_set_ptr);

		/* Print the local bounds on the system variables */

		get_os() << "Location-wise bounds on the variables:" << std::endl
				<< std::endl;
		for (loc_var_bounds_map::const_iterator it = loc_bounds_map.begin(); it
				!= loc_bounds_map.end(); ++it) {
			hybrid_automata::context_automaton_name_formatter lform("",
					get_context());
			get_os() << "Location: " << lform << it->first << std::endl; // prints the location constraint

			var_bounds_map v_bounds = it->second;
			for (variable_id_list::const_iterator v_it = vars.begin(); v_it
					!= vars.end(); ++v_it) {
				get_os() << variable(*v_it) << ": " << v_bounds[*v_it]
						<< std::endl;
			}
			get_os() << std::endl;
		}
	} else {
		get_os() << "empty set" << std::endl;
	}
}

void INTV_continuous_set_formatter<continuous::support_function_provider>::output(
		output_formatter& of, const continuous::support_function_provider& c) {
	typedef double scalar_type;
	typedef math::numeric::interval<scalar_type> bound_interval;

	bool found = true;
	variable_id_list vars;

	if (of.get_output_variables().empty()) {
		variable_id_set vars_set = c.get_variable_ids();
		/* Make a variable id list from variable id set */
		for(variable_id_set::const_iterator it = vars_set.begin(); it!=vars_set.end(); ++it){
			vars.push_back(*it);
		}
	} else {
		vars = of.get_output_variables();
		found = true;
	}
	if (!found) {
		throw std::runtime_error(
				"can't handle these dimensions for support_function_provider");
	}
	else{
		scalar_type t_s, t_e;
		bool is_empty_l, is_bounded_l, is_empty_u, is_bounded_u;
		math::vdom_vector<scalar_type> l_s,l_e,support_vec;

		positional_vdomain set_dom(c.get_variable_ids());


		bound_interval v_bound;
		variable_id v_id;

		for(variable_id_list::const_iterator it = vars.begin();it!=vars.end(); ++it){
			v_id = *it;
			l_e = math::vdom_vector<scalar_type>(set_dom);
			l_s = math::vdom_vector<scalar_type>(set_dom);

			l_e.set_coeff_with_id(v_id,scalar_type(1));
			c.compute_support(l_e,t_e,support_vec,is_empty_u,is_bounded_u);
			l_s.set_coeff_with_id(v_id,scalar_type(-1));
			c.compute_support(l_s,t_s,support_vec,is_empty_l,is_bounded_l);

			v_bound = bound_interval();
			if(is_bounded_l)
				v_bound.set_lower(-t_s);
			if(is_bounded_u)
				v_bound.set_upper(t_e);
			if (is_empty_l || is_empty_u)
				v_bound.set_empty();

			of.get_os() << variable(v_id) << ": " << v_bound << std::endl;
		}
	}

	of.get_os() << std::endl;
}

void INTV_formatter::output(const continuous::continuous_set& c) {
//	context_variable_formatter form(get_context());
//	get_os() << form;

	if (const continuous::spacetime_flowpipe<global_types::float_type> * p =
			dynamic_cast<const continuous::spacetime_flowpipe<
					global_types::float_type> *>(&c)) {
		output(*p);
	} else {
		// try converting to a support_function_provider
		continuous::continuous_set_ptr c_ptr(c.clone());
		continuous::support_function_provider::ptr sup =
				convert_to_support_function_provider(c_ptr);

		base_class::output(*sup);
	}
}

void INTV_formatter::output(
		const continuous::spacetime_flowpipe<global_types::float_type>& flowp) {
	using namespace continuous;

	typedef global_types::float_type scalar_type;

	finite_hyperbox<scalar_type> bbox = compute_finite_bounding_box<scalar_type>(flowp);

	output(bbox);
}


}
