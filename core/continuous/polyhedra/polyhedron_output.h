#ifndef POLYHEDRON_OUTPUT_H_
#define POLYHEDRON_OUTPUT_H_

#include "core/continuous/polyhedra/vertex_enumerator_bf.h"
#include "core/continuous/polyhedra/vertex_enumerator_3D.h"

#include "hyperbox/bounding_box.h"
#include "math/vdom/vdom_matrix_utility.h"
#include "math/vdom/vdom_vector_utility.h"

#include "polyhedron_operators.h"

namespace continuous {

template<typename scalar_type>
polyhedron<scalar_type>* compute_image(
		const polyhedron<scalar_type>& c, const math::affine_map<
				scalar_type>& t);

template<typename scalar_type> void polyhedron<scalar_type>::print_textual(
		std::ostream& os) const {
	if (is_empty()) {
		os << "false";
	} else if (is_universe()) {
		os << "true";
	} else {
		typename math::lin_constraint_system<scalar_type>::const_ptr
				cons = get_constraints();
		for (typename math::lin_constraint_system<scalar_type>::const_iterator
				it = cons->begin(); it != cons->end(); ++it) {
			if (it != cons->begin()) {
				os << " & ";
			}
			os << *it;
		}
	}
}
;

template<typename scalar_type> void polyhedron<scalar_type>::print_double_constraints(
		std::ostream& os) const {
	typename math::lin_constraint_system<scalar_type>::const_ptr cons =
			get_constraints();

	typename math::lin_expression<scalar_type>::output_format
			old_format =
					math::lin_expression<scalar_type>::get_output_format();
	math::lin_expression<scalar_type>::set_output_format(
			math::lin_expression<scalar_type>::output_format::space_separated());

	for (typename math::lin_constraint_system<scalar_type>::const_iterator
			it = cons->begin(); it != cons->end(); ++it) {
		if (it != cons->begin()) {
			os << std::endl;
		}
		os << it->get_l();
	}

	math::lin_expression<scalar_type>::set_output_format(old_format);
}
;

template<typename scalar_type> void polyhedron<scalar_type>::print_double_generators(
		std::ostream& os) const {

	typename math::lin_constraint_system<scalar_type>::const_ptr cons =
			get_constraints();

	unsigned int dim=cons->get_variable_ids().size();
	//testing
	//cout << "The constraints after projection:" << std::endl << cons << std::endl;
	// @todo should use some kind of general get_vertice function
	typename math::vertex_enumerator_bf<scalar_type>::vertice_vector
			vertices;
	typename math::vertex_enumerator_bf<scalar_type>::face_vector faces;
	if (dim == 3)
		math::vertex_enumerator_3D<scalar_type>(*cons, vertices, faces);
	else {
		math::vertex_enumerator_bf<scalar_type>(*cons, vertices, faces);
		//std::cout << std::endl << "vertex_enumerator_bf" << std::endl << vertices << std::endl << faces << std::endl;
	}

	if (vertices.begin() != vertices.end()) {
		/* in 2D, sort vertices ccw */
		if (dim == 2) {
			math::counter_clockwise_sorter<scalar_type>::angle_sort_2D(
					vertices);
		}

		typename math::vdom_vector<scalar_type>::output_format
				old_format =
						math::vdom_vector<scalar_type>::get_output_format();
		math::vdom_vector<scalar_type>::set_output_format(
				math::vdom_vector<scalar_type>::output_format::space_separated());
		for (typename math::vertex_enumerator_bf<scalar_type>::vertice_vector::const_iterator
				it = vertices.begin(); it != vertices.end(); ++it) {
			if (it != vertices.begin()) {
				os << std::endl;
			}
			os << *it;
		}

		// for 2D, repeat the last vertice to close the polygon
		if (cons->get_variable_ids().size() == 2) {
			os << std::endl;
			os << *vertices.begin();
		}

		math::vdom_vector<scalar_type>::set_output_format(old_format);
	}

}
;

template<typename scalar_type> void polyhedron<scalar_type>::print_JVX(
		std::ostream& os, variable_id_list vil) const {
	unsigned int nb_vars = vil.size();
	bool scale_to_unit_bbox = false;

	// GF20110302: turned off scaling because it causes bad constraints
	//             like 1e-15*z <= 1e-15, which are interpreted as 0*z <= 0.
	// 			   needs some more work (increase minimum diameter, for example)
	// Note that turning scaling off significantly increased speed.

	if (nb_vars > 3 || nb_vars < 2) {
		throw basic_exception("JVX output only possible for 2 or 3 dimensions, not for "+to_string(nb_vars)+'.');
	} else {
		typename math::lin_constraint_system<scalar_type>::const_ptr
				cons = get_constraints();

		// reorder if variable list not empty
		positional_vdomain vil_dom;
		if (!vil.empty()) {
			vil_dom=positional_vdomain(vil);
			// create a copy of the constraints and reorder them
			typename math::lin_constraint_system<scalar_type>::ptr newcons
					= typename math::lin_constraint_system<scalar_type>::ptr(
							new math::lin_constraint_system<scalar_type>(*cons));
			newcons->reorder(vil_dom);
			cons = newcons;
		}

		// scaling Matrix and vector
		math::vdom_vector<scalar_type> scale_inv;
		math::vdom_vector<scalar_type> scale_c;

		if (scale_to_unit_bbox) {
			hyperbox<scalar_type> bbox = compute_bounding_box<scalar_type>(*this);
			// remap to the domain of vil if applicable
			if (!vil.empty()) {
				bbox.reorder(vil_dom);
			}
			const positional_vdomain& dom=bbox.domain();
			// scaling vector, default is 1
			math::vdom_vector<scalar_type> s(dom,scalar_type(1));
			scale_inv = math::vdom_vector<scalar_type>(dom,scalar_type(1));
			scale_c = bbox.compute_finite_center();

			math::vdom_vector<scalar_type> u=bbox.get_finite_u();
			math::vdom_vector<scalar_type> l=bbox.get_finite_l();

			// compute scaling factors
			for (unsigned int i = 0; i<s.size(); ++i) {
				scalar_type d=u[i]-l[i]; // always nonnegative
				if (math::numeric::is_GT(d,scalar_type(0))) {
					s[i]=scalar_type(2)/d;
					scale_inv[i]=d/scalar_type(2);
				}
			}
			// create an affine map
			math::vdom_matrix<scalar_type> M=diagonal_matrix(dom,dom,s);
			math::vdom_vector<scalar_type> cnew=M*scale_c;
			math::affine_map<scalar_type> map(M,-cnew);

			typename polyhedron<scalar_type>::ptr scale_poly;
			scale_poly=typename polyhedron<scalar_type>::ptr(compute_image<scalar_type>(*this,map));
			cons = scale_poly->get_constraints();

			//std::cout << "scaling vector:" << s << " unscale:" << scale_inv << std::endl;
		}


		typename math::vertex_enumerator_3D<scalar_type>::vertice_vector
				vertices;
		typename math::vertex_enumerator_3D<scalar_type>::face_vector
				faces;
		//math::vertex_enumerator_3D<scalar_type>(*cons, vertices, faces);
		math::vertex_enumerator_3D<global_types::precise_float_type>::compute<scalar_type>(*cons, vertices, faces);
		//using namespace math;
		//std::cout << vertices << faces << std::endl << std::flush;
		if (vertices.size() > 0 && faces.size() > 0) {
			os << "<geometry name=\"polyhedron\">" << std::endl;

			os << "<pointSet dim=\"";
			os << get_variable_ids().size();
			os << "\" point=\"hide\">" << std::endl;
			os << "   <points> " << std::endl;

			typename math::vdom_vector<scalar_type>::output_format
					old_format =
							math::vdom_vector<scalar_type>::get_output_format();
			math::vdom_vector<scalar_type>::set_output_format(
					math::vdom_vector<scalar_type>::output_format::space_separated());

			for (typename math::vertex_enumerator_3D<scalar_type>::vertice_vector::const_iterator
					it = vertices.begin(); it != vertices.end(); ++it) {
				if (it != vertices.begin()) {
					os << std::endl;
				}
				os << "      <p> ";
				if (scale_to_unit_bbox) {
					os << element_product(scale_inv,(*it)) + scale_c;
				} else {
					os << *it;
				}
				os << " </p>";
			}
			math::vdom_vector<scalar_type>::set_output_format(
					old_format);

			os << std::endl << "   </points>" << std::endl;
			os << "</pointSet>" << std::endl;
			os << "<faceSet face=\"show\" edge=\"hide\">" << std::endl;
			os << "   <faces>" << std::endl;

			/*
			typename math::lin_constraint_system<scalar_type>::const_iterator
					con_it = cons->begin();
			for (typename math::vertex_enumerator_3D<scalar_type>::face_vector::iterator
					it = faces.begin(); it != faces.end(); ++it) {
				assert(con_it!=cons->end());
				if (it != faces.begin()) {
					os << std::endl;
				}
				if (!it->empty()) {
					// sort the vertices counter-clockwise
					math::counter_clockwise_sorter<scalar_type>::angle_sort_3D(
							vertices, *it, con_it->get_normal());
					os << "      <f>";
					for (std::list<unsigned int>::const_iterator jt =
							it->begin(); jt != it->end(); ++jt) {
						os << " " << *jt;
					}
					os << " </f>";
				}
				++con_it;
			}
			*/
			for (typename math::vertex_enumerator_3D<scalar_type>::face_vector::iterator
					it = faces.begin(); it != faces.end(); ++it) {
				if (it != faces.begin()) {
					os << std::endl;
				}
				if (!it->empty()) {
					os << "      <f>";
					for (std::list<unsigned int>::const_iterator jt =
							it->begin(); jt != it->end(); ++jt) {
						os << " " << *jt;
					}
					os << " </f>";
				}
			}

			os << std::endl << "   </faces>" << std::endl;
			os << "</faceSet>" << std::endl;

			os << "</geometry>" << std::endl;
		}
	}
}
;

}

#endif /*POLYHEDRON_OUTPUT_H_*/
