#include "math/ode_solving/traj_simu/postc_simulation.h"
#include "core/symbolic_states/symbolic_state_collection_stl_list.h"
#include "core/continuous/polyhedra/hyperbox/bounding_box.h"
#include "core/continuous/polyhedra/hyperbox/hyperbox.h"
#include "math/ode_solving/traj_simu/continuous_set_simulation.h"
#include "math/numeric/approx_comparator.h"


//#include <boost/math/distributions/uniform.hpp>
//#include <boost/random/mersenne_twister.hpp>
#include <boost/random.hpp>
//#include <random>

namespace hybrid_automata {

void postc_simulation::add_post_states(const hybrid_automaton_const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set,
		const symbolic_state_ptr& sstate) const {

	if (boost::dynamic_pointer_cast<continuous::polyhedron<scalar_type> >(
			sstate->get_continuous_set())) {
		waiting_result_set->union_assign(select_states_from_polyhedron(sstate));
	} else {
		waiting_result_set->add(sstate);
	}
}

void postc_simulation::add_post_states(const hybrid_automaton_const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set,
		const symbolic_state_collection_const_ptr& sstate_set) const {

	// @todo GF: Why just check the first one? Careful: sstate_set could be empty.
	if (boost::dynamic_pointer_cast<continuous::polyhedron<scalar_type> >(
			(*(sstate_set->begin()))->get_continuous_set())) {
		for (symbolic_state_collection::const_iterator it = sstate_set->begin(); it
				!= sstate_set->end(); ++it) {
			add_post_states(aut, passed_result_set, waiting_result_set, *it);
		}
	} else {
		waiting_result_set->union_assign(symbolic_state_collection_ptr(sstate_set->clone()));
	}

}

symbolic_state_collection_ptr postc_simulation::select_states_from_polyhedron(const symbolic_state_ptr& poly) const {
	using namespace continuous;

	symbolic_state_collection_ptr resul = symbolic_state_collection_ptr( new symbolic_state_collection_stl_list() );

	support_function_provider::const_ptr sfpoly = boost::dynamic_pointer_cast<support_function_provider>(poly->get_continuous_set() );
	if(!sfpoly){
		throw std::runtime_error("Unable to apply select states from polyhedron : not a sf provider");
	}

	// Get the bounding box of the set

	LOGGER_OS(DEBUG7,"select_states_from_polyhedron") << "here we go, trying to get center of " <<  sfpoly  <<  std::endl;
	math::vdom_vector<scalar_type> v = compute_bounding_box<scalar_type>(*sfpoly).compute_finite_center();

	// chose the domain
	positional_vdomain dom=v.domain();

	// Define a small epsilon-box for the error
	// chosen to be equal to abs_tol
	// TODO support vectors of abs_tol

	math::vector<scalar_type> vg= math::vector<scalar_type>(v.size(), my_params.abs_tol );
	//math::vector<scalar_type> vg= math::vector<scalar_type>(v.size(),(scalar_type)1e-15);

	math::vector<scalar_type> vc= math::vector<scalar_type>(v.size(),scalar_type(0));

	finite_hyperbox<scalar_type> bbox = finite_hyperbox<scalar_type>(vc,vg,dom);


	// Construct the trajectory set
	continuous_set_simulation<scalar_type>::ptr cset(
			new continuous_set_simulation<scalar_type> (bbox));


	if(my_uniform_sampling==0){
		//TODO replace this by myhbox

		//TODO select some points more wisely ?
		cset->insert(v,trajectory(scalar_type(0),v));
	}
	else{
		hyperbox<scalar_type> h =  compute_bounding_box<scalar_type>(*sfpoly);

		polyhedron<scalar_type>::const_ptr pol =  boost::dynamic_pointer_cast<const polyhedron<scalar_type> >(sfpoly);
		math::lin_constraint_system<scalar_type>::const_ptr lcs = pol->get_constraints();

		state s;
		int i = 0;
		unsigned int pass = 0;
		while (i<my_uniform_sampling){
			pass++;
			s = pick_uniform_state_from_hyperbox(h);
			//std::cout << "picking: " << s << std::endl;
			if( math::maybe(lcs->is_satisfied(s)) ){
				cset->insert(s,trajectory(scalar_type(0),s)	);
				i++;
			}
		}
	}

	symbolic_state_ptr sstate(new symbolic_state(poly->get_discrete_set(),cset));

	resul->add(sstate);

	return resul;
}

postc_simulation::state postc_simulation::pick_uniform_state_from_hyperbox(const continuous::hyperbox<scalar_type>& hyp) const{
	const state & l(hyp.get_finite_l());
	const state & u(hyp.get_finite_u());

	state res(l.domain());
	assert (l.size()==u.size());
	for(int i=0; i<l.size();i++){
		// for the dimensions with zero width, pick the value without calling random
		if(math::numeric::is_MEQ(l[i],u[i]) ){
			res[i] = l[i];
		}
		else{
			scalar_type r = gen_random_scalar(l[i],u[i]);
			res[i] =r;
		}

	}

	return res;

}

postc_simulation::scalar_type postc_simulation::gen_random_scalar(scalar_type min, scalar_type max) const {
	// why float???

    boost::uniform_real<float> u((float)min, (float)max);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > gen(my_rng, u);
    float ff=gen();
    //printf("gen() donne = %e\n",ff);
    return scalar_type(ff);
}


}
