/*
 * choose_directions.cpp
 *
 *  Created on: Dec 7, 2009
 *      Author: frehse
 */

#include "core/continuous/support_function/template_directions/choose_directions.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_constructors.h"

namespace continuous {
namespace support_function {

void direction_chooser::add_directions(const std::string& opt, const std::string& context, const parser::parse_policy& ppol) {
	if (opt == "box") {
		my_types.insert(BOX);
	} else if (opt == "oct") {
		my_types.insert(OCTAGONAL);
	} else if (opt.substr(0, 3) == "uni") {
		my_types.insert(UNIFORM);
		my_number += from_string<unsigned int> (opt.substr(3));
	} else if (opt.substr(0, 1) == "{" && opt.substr(opt.length() - 1, 1)
			== "}") {
		// convert the string to a polyhedron
		constr_polyhedron<scalar_type>::ptr dir_poly =
				construct_constr_polyhedron<scalar_type> (opt.substr(1,
						opt.length() - 2),context, ppol);

		// extract directions from the polyhedron
		math::lin_constraint_system<scalar_type>::ptr cons = math::lin_constraint_system<scalar_type>::ptr(
				new math::lin_constraint_system<scalar_type>(*dir_poly->get_constraints()));
		// add the directions already in the store so that we can unify their domains
		for (direction_store::const_iterator dit = my_store.begin(); dit
				!= my_store.end(); ++dit) {
			math::vdom_vector<scalar_type> vec = dit->first;
			math::lin_expression<scalar_type> l(vec, scalar_type(1));
			math::lin_constraint<scalar_type> con(l,LE);
			cons->push_back(con);
		}
		cons->unify_domains();
		for (math::lin_constraint_system<scalar_type>::const_iterator it = cons->begin(); it
				!= cons->end(); ++it) {
			direction_vector v = it->get_normal();
			my_store.insert(v);
			if (it->is_equality()) {
				my_store.insert(-v);
			}
		}
//std::cout << "in context " + context + " defined directions ";
//my_store.print(std::cout); std::cout << std::endl;
		IFLOGGER(DEBUG3) {
			std::stringstream ss;
			logger::copyfmt_to(ss);
			for (direction_store::const_iterator dit = my_store.begin(); dit
					!= my_store.end(); ++dit) {
				if (dit != my_store.begin())
					ss << ", ";
				math::vdom_vector<scalar_type> vec = dit->first;
				ss << math::lin_expression<scalar_type>(vec, scalar_type(0));
			}
			LOGGER(DEBUG3,"choose_directions","set user directions to " +ss.str());
		}
		my_types.insert(POLYH);
	} else {
		throw std::runtime_error("unknown direction type " + opt);
	}
}

direction_chooser::type_store direction_chooser::my_types;
direction_chooser::direction_store direction_chooser::my_store;
unsigned int direction_chooser::my_number = 0;

}
}
