/*
 * spacetime_flowpipe_utility.h
 *
 *  Created on: Mar 21, 2013
 *      Author: notroot
 */

#ifndef SPACETIME_FLOWPIPE_UTILITY_H_
#define SPACETIME_FLOWPIPE_UTILITY_H_

#include "utility/logger.h"

namespace continuous {

template<typename scalar_type>
finite_hyperbox<scalar_type> compute_finite_bounding_box(const continuous::spacetime_flowpipe<global_types::float_type>& flowp) {

	typedef spacetime_flowpipe<scalar_type> flowpipe;

	const positional_vdomain& dom = flowp.domain();
	typename flowpipe::direction d(dom);
	size_t N = dom.size();
	typename finite_hyperbox<scalar_type>::vdom_vector_type l(dom),u(dom),c,g;
	if (true) {
		// compute bounding box directions
		// This is just to make sure; if already present no computational overhead.

		// fix const cast
		flowpipe& nonconst_flowp = const_cast<flowpipe&>(flowp);

		for (size_t i = 0; i < N; ++i) {
			// don't compute flowpipes, just obtain what's already computed
			typename flowpipe::error_type err(1e10, 1e10);
			LOGGER(DEBUG5, "compute_finite_bounding_box",
					"computing evolution of variable "+dom.get_variable(i).get_name()+" for bounding box");
			d = typename flowpipe::direction(dom);
			d[i] = 1.0;
			const typename flowpipe::spacetime_plif::annotated_plif& evo_ann_pos =
					nonconst_flowp.get_or_compute_evolution(d, err);
			u[i] = plif::supremum(evo_ann_pos.first.get_upper());

//				plif_graph(evo_ann_pos.first, "X", "/tmp/test_gen");
			d[i] = -1.0;
			const typename flowpipe::spacetime_plif::annotated_plif& evo_ann_neg =
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

}


#endif /* SPACETIME_FLOWPIPE_UTILITY_H_ */
