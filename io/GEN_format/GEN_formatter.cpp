/*
 * GEN_formatter.cpp
 *
 *  Created on: Oct 12, 2010
 *      Author: frehse
 */

#include "GEN_formatter.h"

#include "core/symbolic_states/symbolic_state.h"

namespace io {

void GEN_formatter::output(
		const hybrid_automata::symbolic_state_collection& sstates) {
	if (sstates.is_empty()) {
		// generate a warning that the set is empty
		basic_warning("GEN format output",
				"cannot represent empty set in GEN format",
				basic_warning::AMBIGUOUS_OUTPUT);
	} else {
		base_class::output(sstates);
	}
}

void GEN_formatter::output(const hybrid_automata::symbolic_state& sstate) {
	if (sstate.get_discrete_set()->is_empty()) {
		// generate a warning that the set is empty
		basic_warning("GEN format output",
				"cannot represent empty set in GEN format",
				basic_warning::AMBIGUOUS_OUTPUT);
	} else {
		base_class::output(sstate);
	}
}

void GEN_continuous_set_formatter<continuous::support_function_provider>::output(
		output_formatter& of, const continuous::support_function_provider& c) {
	typedef double scalar_type;
	/*
	 if (of.get_output_variables().empty()) {
	 c.print_double_generators(of.get_os());
	 } else {
	 typename continuous::polyhedron<scalar_type>::ptr cset(c.clone());
	 cset->project_to_variables(of.get_output_variables());
	 cset->print_double_generators(of.get_os());
	 } */
	variable_id x_id = 0;
	variable_id y_id = 0;
	bool found = false;
	if (of.get_output_variables().empty()) {
		variable_id_set vis = c.get_variable_ids();
		if (vis.size() == 2) {
			x_id = *vis.begin();
			y_id = *vis.rbegin();
			found = true;
		}
	} else if (of.get_output_variables().size() == 2) {
		x_id = *of.get_output_variables().begin();
		y_id = *of.get_output_variables().rbegin();
		found = true;
	}
	if (!found) {
		throw std::runtime_error(
				"GEN_formatter:can't handle these dimensions for support_function_provider");
	} else {
		//std::cout << "approx:";
		double eps = static_cast<GEN_formatter&>(of).my_eps;

		try {

			if (c.computes_support_vector()) {
				continuous::support_function::vertice_approx_2D_sv<scalar_type>::sample_list pts;
				pts = continuous::support_function::vertice_approx_2D_sv<
						scalar_type>::approx(c, x_id, y_id, scalar_type(eps));
				//std::cout << pts.size() << std::endl;
				for (continuous::support_function::vertice_approx_2D_sv<
						scalar_type>::sample_list::const_iterator it =
						pts.begin(); it != pts.end(); ++it) {
					of.get_os() << convert_element<double>(it->a) << " "
							<< convert_element<double>(it->b) << std::endl;
				}
				// repeat the first point to close the curve
				of.get_os() << convert_element<double>(pts.begin()->a) << " "
						<< convert_element<double>(pts.begin()->b) << std::endl;
			} else {
				continuous::support_function::vertice_approx_2D<scalar_type>::sample_list pts;
				pts = continuous::support_function::vertice_approx_2D<
						scalar_type>::approx(c, x_id, y_id, scalar_type(eps));
				//std::cout << pts.size() << std::endl;
				for (continuous::support_function::vertice_approx_2D<scalar_type>::sample_list::const_iterator it =
						pts.begin(); it != pts.end(); ++it) {
					of.get_os() << convert_element<double>(it->a) << " "
							<< convert_element<double>(it->b) << std::endl;
				}
				// repeat the first point to close the curve
				of.get_os() << convert_element<double>(pts.begin()->a) << " "
						<< convert_element<double>(pts.begin()->b) << std::endl;
			}

		} catch (continuous::unbounded_set_exception& e) {
			// prepare a user-readable output of the variable bounds
			std::stringstream ss;
			INTV_formatter intvf(ss);
			intvf.set_output_variables(of.get_output_variables());
			intvf.set_context(of.get_context());
			intvf.output(c);
			// generate a warning that the set is unbounded
			basic_warning("GEN format output",
					"cannot represent unbounded set in GEN format. Bounds:\n"
							+ ss.str(), basic_warning::AMBIGUOUS_OUTPUT);

		} catch (continuous::empty_set_exception& e) {
			// generate a warning that the set is empty
			basic_warning("GEN format output",
					"cannot represent empty set in GEN format",
					basic_warning::AMBIGUOUS_OUTPUT);
		}
	}
	// terminate poly with an empty line
	of.get_os() << std::endl;
}

void GEN_formatter::output(
		const math::trajectory<global_types::float_type>& traj) {
	variable_id_list vars;
	if (get_output_variables().empty()) {
		const variable_id_set& vis = traj.get_variable_ids();
		for (variable_id_set::const_iterator it = vis.begin(); it != vis.end();
				++it) {
			vars.push_back(*it);
		}
	} else {
		vars = get_output_variables();
	}
	// Construct a vector of indices
	//std::vector<unsigned int> indices(vars.size());
	std::vector<unsigned int> indices;
	unsigned int max_index = std::numeric_limits<unsigned int>::max();
	unsigned int i = 0;
	for (variable_id_list::const_iterator it = vars.begin(); it != vars.end();
			++it) {
		unsigned int j;
		if (traj.domain().in_domain(variable(*it), j)) {
			indices.push_back(j);
		} else {
			indices.push_back(max_index);
		}
		++i;
	}

	// only output something if there is at least one variable to show
	if (i > 0) {
		double val;
		for (i = 0; i < traj.size(); ++i) {
			for (unsigned int j = 0; j < indices.size(); ++j) {
				if (indices[j] != max_index) {
					val = traj.get_states()(i, indices[j]);
				} else {
					val = 0;
				}
				get_os() << val << " ";
			}
			get_os() << std::endl;
		}
		get_os() << std::endl;
	}
}

void GEN_formatter::output(
		const continuous::continuous_set_simulation<global_types::float_type>& cset_simu) {
	for (continuous::continuous_set_simulation<global_types::float_type>::const_iterator it =
			cset_simu.begin(); it != cset_simu.end(); ++it) {
		output(it->second);
	}
}

void GEN_formatter::output(
		const continuous::polyhedron_collection<global_types::float_type>& polys) {
	for (continuous::polyhedron_collection<global_types::float_type>::const_iterator it =
			polys.begin(); it != polys.end(); ++it) {
		output(**it);
	}
}

void GEN_formatter::output(
		const continuous::spacetime_flowpipe<global_types::float_type>& flowp) {
	typedef continuous::spacetime_flowpipe<global_types::float_type> flowpipe;

	if (true) {

		// compute bounding box directions
		// This is just to make sure; if already present no computational overhead.
		const positional_vdomain& dom = flowp.domain();
		flowpipe::direction d(dom);
		size_t N = dom.size();
		flowpipe::spacetime_plif::annotated_plif evo_ann_pos,
				evo_ann_neg;

		// fix const cast
		flowpipe& nonconst_flowp = const_cast<flowpipe&>(flowp);

		// @todo If we do this, it should only be over the output variables
		std::string temp_fname = continuous::spacetime_flowpipe<double>::temp_filename;
		if (temp_fname != "") {
			LOGGER(MEDIUM, "continuous_post_stc::post",
									"saving "+to_string(N)+" flowpipe projections to "+temp_fname);
			static int plot_counter = 0;
			++plot_counter;
			for (size_t i = 0; i < N; ++i) {
				// don't compute flowpipes, just obtain what's already computed
				flowpipe::error_type err(1e10, 1e10);
				LOGGER(DEBUG5, "continuous_post_stc::post",
						"plotting evolution of variable "+dom.get_variable(i).get_name());
				d = flowpipe::direction(dom);
				d[i] = 1.0;
				evo_ann_pos = nonconst_flowp.get_or_compute_evolution(d, err);
//				plif_graph(evo_ann_pos.first, "X", "/tmp/test_gen");
				d[i] = -1.0;
				evo_ann_neg = nonconst_flowp.get_or_compute_evolution(d, err);
//				plif_graph(-(evo_ann_neg.first), "X", "/tmp/test_gen");

				// get inner and outer bounds
				plif::piecewise_linear_interval_function h_outer(
						-evo_ann_neg.first.get_upper(),
						evo_ann_pos.first.get_upper());
				const plif::piecewise_linear_function& inner_low =
						-evo_ann_neg.first.get_lower();
				const plif::piecewise_linear_function& inner_upp =
						evo_ann_pos.first.get_lower();
				// the inner can't go past the outer
//				inner_low = plif::pointwise_maximum(inner_low,
//						h_outer.get_lower());
//				inner_upp = plif::pointwise_minimum(inner_upp,
//						h_outer.get_upper());
				plif::piecewise_linear_interval_function h_inner(inner_low,
						inner_upp);

				// the inner could be empty, so split
				std::vector<plif::piecewise_linear_interval_function> h_inner_split = plif::split(h_inner);

				// write lower to file
				std::string filename = temp_fname + "_"// "/tmp/GEN_formatter_supp_evo_t_"
						+ to_string(dom.get_variable(i).get_name()) + "_nb_"
						+ to_string(plot_counter);
				std::string filename_inner = filename + "_inner.txt";
				std::ofstream ofile;
				ofile.open(filename_inner.c_str());
				for (std::vector<plif::piecewise_linear_interval_function>::const_iterator it=h_inner_split.begin();it!=h_inner_split.end();++it) {
					plif_vertex_output_closed(ofile, *it);
					ofile << std::endl; // separate each piece
				}
				ofile.close();
				// write upper to file
				std::string filename_outer = filename + "_outer.txt";
				;
				ofile.open(filename_outer.c_str());
				plif_vertex_output_closed(ofile, h_outer);
				ofile.close();

				// plot both
				std::string opt_string = "";
				std::string format = "png";
				std::string epilog;
				if (format != "X") {
					epilog = "> " + filename + "." + format;
				}
				//  -S 4
				std::string command_line = "graph -T" + format + " -C -B -q 0.5"
						+ opt_string + " -m1 " + filename_outer + " -s -m2 "
						+ filename_inner + " " + epilog;
				int res = std::system(command_line.c_str());
			}
		}
	}

	flowpipe::cut_point_method m;
	if (my_eps > 1e-6) {
		m.type = flowpipe::cut_point_method::MIN_CONCAVE_PIECES;
		m.approx_error = flowpipe::error_type(0, my_eps);
		m.lower_is_error_reference = false;
	} else {
		m.type = flowpipe::cut_point_method::ALL_PIECES;
	}
	continuous::polyhedron_collection<global_types::float_type> polys =
			flowp.compute_outer_polyhedra(m);

	double old_my_eps = my_eps;
	// set tolerance to zero before outputting, otherwise error will accumulate
	my_eps = 0;
	output(polys);
	my_eps = old_my_eps;
}

}
