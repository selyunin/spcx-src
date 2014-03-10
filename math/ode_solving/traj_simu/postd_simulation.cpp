#include "math/ode_solving/traj_simu/postd_simulation.h"
#include "core/discrete/discrete_set.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/transition.h"
#include "core/hybrid_automata/location.h"
#include "core/discrete/singleton_set.h"
#include "core/hybrid_automata/hybrid_automaton_utility.h"
#include "core/symbolic_states/symbolic_state_collection.h"

#include "core/continuous/continuous_dynamics/typed_dynamics.h"
#include "core/continuous/continuous_dynamics/ode_dynamics.h"

#include "core/continuous/polyhedra/polyhedron.h"
#include "math/vdom/state_functor_utility.h"

#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"

#include "math/ode_solving/traj_simu/continuous_set_simulation.h"

//debug
#include "utility/stl_helper_functions.h"
#include "utility/logger.h"
#include <iostream>
#include <sstream>


/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
}


namespace hybrid_automata {


discrete_post* postd_simulation::clone(){
	return  new postd_simulation(*this);
}

continuous::continuous_set_collection postd_simulation::post(
		const jump_constraints& trans,
		continuous::continuous_set_const_ptr source_inv,
		continuous::continuous_set_const_ptr target_inv,
		continuous::continuous_set_const_ptr cset) const {

		throw std::runtime_error("feature not implemented");
}


void postd_simulation::add_post_states_trans(const hybrid_automaton_const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set,
		const transition_id& trans, const symbolic_state_ptr& sstate) const{
	throw std::runtime_error("feature not implemented");
}


void postd_simulation::add_post_states(const hybrid_automaton_const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set,
		const label_id& lab,
		const symbolic_state_ptr& sstate) const{
	throw std::runtime_error("feature not implemented");
}


void postd_simulation::add_post_states(const hybrid_automaton_const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set,
		const label_id_set& lab_set, const symbolic_state_ptr& sstate) const{

	typedef simple_iterators::collection_const_base<transition_id>::const_iterator
			transition_const_iterator;

	typedef simple_iterators::collection_const_base<location_id>::const_iterator
			location_const_iterator;

	if(sstate){
		discrete::discrete_set::const_ptr d=sstate->get_discrete_set();

		std::pair<transition_const_iterator, transition_const_iterator> pr;
		for (discrete::discrete_set::const_iterator d_it = d->begin(); d_it != d->end(); ++d_it) {
			std::list<transition_id> out_list;
			location_id_set locs = aut->get_locations(*d_it);
			for (location_id_set::const_iterator l_it = locs.begin(); l_it != locs.end(); ++l_it) {
				for (label_id_set::const_iterator lab_it = lab_set.begin(); lab_it
								!= lab_set.end(); ++lab_it) {

					pr = aut->get_outgoing_transitions(*l_it, *lab_it);
					copy(pr.first, pr.second, std::back_inserter(out_list));
				}
				add_post_states(aut,passed_result_set,waiting_result_set,out_list,sstate);

			}
		}
	}
}

void postd_simulation::add_post_states(const hybrid_automaton_const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set,
		const std::list<transition_id>& trans_set, const symbolic_state_ptr& sstate) const{
	using namespace continuous;

	if (sstate) {

		std::vector <jump_constraints> trans=  std::vector<jump_constraints >();
		trans.reserve( trans_set.size() );
		std::vector<continuous_set_const_ptr> targets_inv = std::vector<continuous_set_const_ptr> ();
		targets_inv.reserve( trans_set.size() );

		continuous_set_const_ptr cset;


		// extraction of source invariant + dynamic
		//TODO check if it is ok in general
		//std::list<transition_id>::const_iterator it = trans_set.begin(); OLD VERSION, PB WITH NO TRANSISTION
		//transition::ptr trans_ptr = aut->get_transition(*it);
		//const location_id& s_id = trans_ptr->get_source();
		const location_id_set locs = aut->get_locations(*(sstate->get_discrete_set()->begin()));
		if (locs.empty()) {
			throw std::runtime_error("postd_simulation::add_post_states(l.121) : no location to apply postd!");
		}
		// All locs in dset are supposed to have the same dynamics,
		// so let's just choose the dynamics of the first location
		// in the set
		const location_id& s_id = *locs.begin();

		LOGGER_OS(HIGH,__FUNCTION__) << "Computing trajectories in location " << canonicalize_location_constraint(aut, s_id);

		continuous_set::const_ptr s_inv = aut->get_location(s_id)->get_time_constraints().get_invariant();
		continuous_set::const_ptr source_cset = sstate->get_continuous_set();
		continuous_dynamics::const_ptr s_dyn = aut->get_location(s_id)->get_time_constraints().get_dynamics();

		std::vector<discrete::discrete_set::ptr> dset;
		dset.reserve(trans_set.size());



		for (std::list<transition_id>::const_iterator it=trans_set.begin(); it != trans_set.end(); ++it) {
			// Get the target location of the transition.

			transition::ptr trans_ptr = aut->get_transition(*it);

			const location_id& t_id = trans_ptr->get_target();

			location_constraint_set lcs = location_constraint_set(aut->get_id(),
							t_id);
			discrete::discrete_set::ptr ds = discrete::discrete_set::ptr(
							new discrete::singleton_set(lcs));
			dset.push_back(ds);

			try {

				// Get the invariant of the target location
				continuous_set::const_ptr
						t_inv =	aut->get_location(t_id)->get_time_constraints().get_invariant();
				const jump_constraints& jump =
						trans_ptr->get_jump_constraints();

				trans.push_back(jump);
				targets_inv.push_back(t_inv);

			} catch (std::exception& e) {
				std::stringstream s;
				logger::copyfmt_to(s);
				s << canonicalize_location_constraint(aut, s_id);
				s << " to location ";
				s << canonicalize_location_constraint(aut, t_id);
				throw basic_exception(
						"could not apply discrete post for transition with id "
								+ to_string((*it)) + ", label "
								+ named_label::get_name(trans_ptr->get_label())
								+ ", from location " + s.str() + ".", e);
			}

		}

		// a structure for storing the trajectories computed for the roots in cset_roots
		continuous_set_simulation<scalar_type>::ptr source_trajs;
		std::vector< continuous_set::ptr > csets = post(trans, s_inv, s_dyn, targets_inv,
				source_cset, source_trajs);

		// add the trajectories in the source location to the passed result list
		passed_result_set->add(
				symbolic_state::ptr(
						new symbolic_state(sstate->get_discrete_set(),
								source_trajs)));

		int i = 0;
		for (std::vector< continuous_set::ptr  >::iterator it = csets.begin(); it
				!= csets.end(); ++it) {
			continuous_set::ptr cset = *it;
			if (cset && !cset->is_empty()){
				waiting_result_set->add(symbolic_state::ptr(new symbolic_state(dset[i],
										cset)));
			//	printf("DEBUG DEBUG : adding sstate for trans %i \n",i);
			}
			//printf("DEBUG DEBUG : boucle : i= %i \n",i);
			i++;
		}

	}
}

/**
 *
 * Post-operator for simulation, does both continuous and discrete post (because they
 * are highly tied together in the simulation algorithm).
 * Takes as input states, and handles them as root of the trajectory to be computed.
 * Trajectory is obtained by integration of the ODE of the dynamics. While the trajectory is computed,
 * root in guards and invariant are searched and possible jump intervals are sampled to obtain new roots for computation.
 * For now, handles only deterministic dynamics and guard resets
 */
std::vector<continuous::continuous_set::ptr >  postd_simulation::post(
				const std::vector<jump_constraints> & trans,
				const continuous::continuous_set_const_ptr & source_inv,
				const continuous::continuous_dynamics::const_ptr & source_dyn,
				const std::vector<continuous::continuous_set_const_ptr> & targets_inv,
				const continuous::continuous_set_const_ptr & cset,
				continuous::continuous_set_simulation<scalar_type>::ptr& cset_trajs) const{
	using namespace continuous;
	typedef continuous_set_simulation<scalar_type>::trajectory trajectory;
	typedef continuous_set_simulation<scalar_type>::state state;
	typedef math::lin_constraint_system<scalar_type>::ptr lin_cons_sys_ptr;

	LOGGERSW(DEBUG1,"postd_simulation::post","Simulating trajectories for state set");

	assert(trans.size() == targets_inv.size() );

	std::vector<simu::lin_constraint_system::const_ptr> guards = std::vector<simu::lin_constraint_system::const_ptr>();
	guards.reserve(trans.size());

	std::vector<math::state_functor<scalar_type>::const_ptr > g_assign = std::vector<math::state_functor<scalar_type>::const_ptr>();
	g_assign.reserve(trans.size());

	simu::state_jump_results res;
	simu::roots_and_trajectory randt;

	// Casting all necessary informations from generic classes to classes needed by simulation algorithm

	for (int i = 0; i < trans.size(); i++) {

		math::affine_state_functor<scalar_type>::const_ptr state_fun;
		//functors for reset functions
		const reset_affine_transform<scalar_type>* resetg = NULL;
		if (trans[i].get_transform()) {

			resetg = dynamic_cast<const reset_affine_transform<	scalar_type>*> (trans[i].get_transform().get());
			if (trans[i].get_transform().get() && !resetg) {
				std::stringstream ss;
				logger::copyfmt_to(ss);
				ss << "offending transform:" << trans[i].get_transform();
				throw std::runtime_error(
						"postd_simulation::post : transform type not supported\n"
								+ ss.str());
			}

			state_fun = math::affine_state_functor<scalar_type>::const_ptr(
					new math::affine_state_functor<scalar_type>(*resetg));
			g_assign.push_back(state_fun);
		}
		// if inst.reset_trans_ptr null -> identity map
		else {
			state_fun = math::affine_state_functor<scalar_type>::const_ptr(
					new math::affine_state_functor<scalar_type>(
							math::affine_map<scalar_type>()));
			g_assign.push_back(state_fun);
		}

		// Get guards constraints
		const polyhedron<scalar_type>::const_ptr g = boost::dynamic_pointer_cast<
				const polyhedron<scalar_type> >(trans[i].get_guard() );

		if (trans[i].get_guard() && !g) {
			throw std::runtime_error(
					"postd_simulation::post : continuous_set type not supported as invariant");
		}

		math::lin_constraint_system<scalar_type>::const_ptr lcsg;

		if(g){
			lcsg = g->get_constraints();
		}
		else{
			printf("WARNING : null constraints");
			lcsg = math::lin_constraint_system<scalar_type>::const_ptr(
					new math::lin_constraint_system<scalar_type>()  );
		}

		lin_cons_sys_ptr lcsg_i=lin_cons_sys_ptr(new  math::lin_constraint_system<scalar_type>(*lcsg) );


		// Note: the backtransformed invariant is already added during preprocessing

		//		constr_polyhedron<scalar_type>::const_ptr inv_poly =
		//				boost::dynamic_pointer_cast<const constr_polyhedron<scalar_type> >(targets_inv[i]);
//		// Get Invariant constraints
//		if(inv_poly){
//
//			if(resetg){
////			std::cout << "DEBUG INV REVERSE MAP\n BEFORE \n resetg : \n" << *resetg << "\n and constraints : \n"
////							<< *(inv_poly->get_constraints()) << std::endl;
//				constr_polyhedron<scalar_type> inv_poly_rmap = reverse_map(*inv_poly,*resetg);
//				math::lin_constraint_system<scalar_type>::const_ptr lcsg_f= inv_poly_rmap.get_constraints();
//				lcsg_i->push_back(*lcsg_f);
////      std::cout << "AFTER : \n constraints : \n " << *lcsg_f << std::endl;
//				lcsg_i->expand_equalities();
//				lcsg = lcsg_i;
//			}
//			else{
//				math::lin_constraint_system<scalar_type>::const_ptr lcsg_f= inv_poly->get_constraints();
//				lcsg_i->push_back(*lcsg_f);
//				lcsg_i->expand_equalities();
//				lcsg = lcsg_i;
//			}
//		}
//		else throw std::runtime_error("unable to cast invariant to constr_poly\n");

		lcsg_i->expand_equalities();
		guards.push_back(lcsg_i);
	}

	LOGGER_OS(DEBUG7,"postd_simulation") << "DEBUG Computed guards:"
					<< guards << std::endl;

	//idem for invariant (source) constraints
	const polyhedron<scalar_type>::const_ptr s = boost::dynamic_pointer_cast<
			const polyhedron<scalar_type> >(source_inv);

	if (source_inv && !s) {
		throw std::runtime_error(
				"postd_simulation::post : continuous_set type not supported as invariant");
	}

	lin_cons_sys_ptr inv_s = lin_cons_sys_ptr(new math::lin_constraint_system<scalar_type>( *(s->get_constraints())));
	inv_s->expand_equalities();

	LOGGER_OS(DEBUG7,"postd_simulation") << "DEBUG: Computed Invariant:"
					<< inv_s << std::endl;

	//construct functor for "f"

	math::affine_map<scalar_type>::const_ptr Ab = boost::dynamic_pointer_cast<const ode_affine_dynamics<scalar_type> >( source_dyn );

	const positional_vdomain domAb  = Ab->domain();
	variable_id_set domvars = Ab->domain().get_variable_ids();
	variable_id_set codomvars = Ab->codomain().get_variable_ids();
	variable_id_set A2domvars = domvars;
	set_difference_assign(A2domvars, codomvars);

	vdom_matrix<scalar_type> Mzero = vdom_matrix<scalar_type>( positional_vdomain(A2domvars),Ab->domain());
	ode_affine_dynamics<scalar_type>  Ucst =  ode_affine_dynamics<scalar_type>( math::affine_map<scalar_type>(Mzero) );

	ode_affine_dynamics<scalar_type>  AbUcst = ode_affine_dynamics<scalar_type>(math::compose(*Ab, Ucst));

	// get one state to sync domains from the set of ICs to be explored
	continuous_set_simulation<scalar_type>::const_ptr cset_ICs = boost::dynamic_pointer_cast<const continuous_set_simulation<scalar_type> >(cset);
	if(!cset_ICs) throw std::runtime_error(
			"postd_simulation::post : continuous_set type not supported\n");

	continuous_set_simulation<scalar_type>::const_iterator it=cset_ICs->begin();
	const state& x = it->first;
	positional_vdomain xdom = x.domain();

	vdom_matrix<scalar_type> Anew = AbUcst.get_A();
	Anew.reorder(xdom,xdom);
	vdom_vector<scalar_type> bnew = AbUcst.get_b();
	bnew.reorder(xdom);

//	LOGGER_OS(DEBUG7,"postd_simulation") << "DEBUG: dynamics separation between states and inputs"
//		<< "x:" << x << std::endl
//		<< "Anew:" << Anew << std::endl
//		<< "bnew:" << bnew << std::endl;

	//typed_dynamics<scalar_type>::const_ptr dp = boost::dynamic_pointer_cast<const typed_dynamics<scalar_type> >( source_dyn );

	typed_dynamics<scalar_type>::const_ptr dp = typed_dynamics<scalar_type>::const_ptr(new ode_affine_dynamics<scalar_type>(affine_map<scalar_type>(Anew, bnew)));

	if (!dp) {
		std::stringstream ss;
		ss << typeid(*source_dyn).name();
		throw std::runtime_error(
							"postd_simulation::post : dynamics type not supported : "+ss.str()+"\n");
	}

	//std::cout << "type of dynamics : " << typeid(*dp).name() << std::endl;

	//functor that computes f's Jacobian

	typed_dynamics_state_functor<scalar_type> f(dp);
	typed_dynamics_jacobian_functor<scalar_type> Jac(dp);

	cset_trajs=continuous_set_simulation<scalar_type>::ptr(new continuous_set_simulation<scalar_type>(cset_ICs->get_hbox()));

	//All structures constructed, call solving subfunctions and generating result set.

	// obtain a vector of roots to be explored in the next iteration
	std::vector<continuous_set::ptr> result = std::vector<continuous_set::ptr>();
	result.reserve(trans.size());

	LOGGER_ATTACH(DEBUG2,"postd_simulation::post"," with "+to_string(cset_ICs->size())+" states",logger::get_last_id());
	LOGGER_ATTACH(DEBUG4,"postd_simulation::post"," for "+to_string(trans.size())+" outgoing transitions",logger::get_last_id());

	// @todo GF: Why is result different from interv? Can we not work directly with result and then
	// return result?
	// Right now: interv is a vector of continuous_set_simulation<scalar_type>::ptr,
	// while result is a set of base pointers

	continuous_set_simulation<scalar_type>::ptr inter;
	std::vector<continuous_set_simulation<scalar_type>::ptr> interv;
	interv.reserve(trans.size());
	for(int k=0;k<trans.size();k++){
		inter = continuous_set_simulation<scalar_type>::ptr ( new continuous_set_simulation<scalar_type>(cset_ICs->get_hbox()) );
		interv.push_back(inter);
	}
	continuous_set_simulation<scalar_type>::ptr inters;


	// process each of the roots in cset_roots

	// total count of successor points
	unsigned int successor_count = 0;
	unsigned int transition_count = 0;
	for(continuous_set_simulation<scalar_type>::const_iterator it=cset_ICs->begin();it!=cset_ICs->end();++it){
		// the current root
		const state& x_current = it->first;

		// the trajectory associated with the root
		// right now, we only use its time stamp
		const trajectory& traj = it->second;

		// the current time stamp (global time)
		// the trajectory is just a point (for now), so take the first element
		scalar_type t_current = traj.get_time(0);

		// compute the time horizon for the new trajectory
		scalar_type tmax;
		if (mypostc.get_time_horizon() > 0 && mypostc.get_global_time_horizon() > 0)
			tmax = std::min(t_current + mypostc.get_time_horizon(),
					mypostc.get_global_time_horizon());
		else if (mypostc.get_time_horizon() < 0)
			tmax = mypostc.get_global_time_horizon();
		else
			tmax = t_current + mypostc.get_time_horizon();

		// Compute the trajectory starting from the current root
		simu::state_jump_results res1;

		LOGGER_OS(DEBUG7,"postd_simulation") << "DEBUG postd_simulation::post : exploring init cond : \n"
				<< x_current << std::endl << " at time t = " << t_current
				<< " up to t = " << tmax << ", lth = " << mypostc.get_time_horizon() << ", gth = " << mypostc.get_global_time_horizon() << std::endl
				<< "using simu algo " << my_simu_algo << std::endl;
		{

			LOGGERSWOC(DEBUG5,__FUNCTION__,"Simulating trajectories");

			switch (my_simu_algo) {
			case 0:
				res1 = simulation.compute_urgent_traj(f, Jac, x_current,
						t_current, tmax, guards, inv_s);
				break;

			case 1:
				res1 = simulation.state_jump_intervals_solver(f, Jac, x_current,
						t_current, tmax, guards, inv_s);
				break;

			case 2:
				res1 = simulation.state_jump_intervals_derivative(f, Jac, x_current,
						t_current, tmax, guards, inv_s);
				break;

			default:
				throw std::runtime_error(
						"postd_simulation::post : unsupported internal simulation algorithm\n");
				break;
			}
		}

		// add the trajectory to the set of trajectories for the source location
		cset_trajs->insert(x_current,res1.traj);

		// obtain a set of points of the current trajectory to be mapped by the transitions
		simu::roots_and_trajectory res2 = simulation.select_next_traj_states(res1);

		// Obtain the future roots as a vector of trajectories.
		// Each trajectory corresponds to a set of roots
		//std::cout<< "Postd::post,  trajectory : \n"<<res2.traj <<" and selected future roots : \n" << res2.roots;
		//std::cout<<"DEBUG postd_simulation::post : selected futuure roots : \n" << res2.roots << std::endl;
		res2.roots = simulation.image_next_traj_states(res2.roots,g_assign);

		// add the new roots
		assert(res2.roots.size()==trans.size());
		// for each transition, add the set of roots
		for(int j=0;j<res2.roots.size() ; j++){
			// Each state in the trajectory is a new root.
			// Add them as roots to the new sets
			for (int k=0; k<res2.roots[j].size(); k++){
				const continuous_set_simulation<scalar_type>::state& x = res2.roots[j].get_state(k);
				scalar_type t = res2.roots[j].get_time(k);
				interv[j]->insert(x,trajectory(t,x));
				//	if ( res2.roots[j].size()>0 ) printf("\n Inserting non-empty root at index j=%i \n",j);
				++successor_count;
			}
			if (res2.roots[j].size()>0)
				++transition_count;
		}
	}

	LOGGER(DEBUG3,"postd_simulation::post","found a total of "+to_string(successor_count)+" successor states");
	LOGGER_ATTACH(DEBUG4,"postd_simulation::post"," for "+to_string(transition_count)+" outgoing states",logger::get_last_id());

	for(int k=0;k<trans.size();k++){
			result.push_back(interv[k]);
		}

	return result;
}



}
