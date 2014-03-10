#ifndef vdom_VECTOR_UTILITY_H_
#define vdom_VECTOR_UTILITY_H_

#include <algorithm>

namespace math {

/** Compute the average over an STL container. If the container is empty,
 * value_type() is returned.
 * value_type must have the operator/(const value_type& v, const int& c)
 * or equivalent. 
 * */
template<typename container_type> typename container_type::value_type compute_average(
		const container_type& container) {
	if (container.begin() != container.end()) {
		int count = 1;
		typename container_type::const_iterator it = container.begin();
		typename container_type::value_type res(*it);
		for (++it; it != container.end(); ++it) {
			res += *it;
			++count;
		}
		res /= count;
		return res;
	} else {
		typename container_type::value_type def;
		return def;
	}
}

/** Compute the average over selected elements in an STL random access container.
 * The elements are given in a container indices. 
 * If the container is empty,
 * value_type() is returned.
 * value_type must have the operator/(const value_type& v, const int& c)
 * or equivalent. 
 * */
template<typename container_type, typename index_container_type> typename container_type::value_type compute_average(
		const container_type& container, const index_container_type& indices) {
	if (indices.begin() != indices.end()) {
		int count = 1;
		typename index_container_type::const_iterator it = indices.begin();
		typename container_type::value_type res(container[*it]);
		for (++it; it != indices.end(); ++it) {
			res += container[*it];
			++count;
		}
		res /= count;
		return res;
	} else {
		typename container_type::value_type def;
		return def;
	}
}

template<typename scalar_type> class counter_clockwise_sorter {
public:
	typedef vdom_vector<scalar_type> vertice_type;
	typedef std::vector<vertice_type> vertice_vector;
	typedef std::list<unsigned int> index_list;

	/** Comparator for sorting 2D or 3D points in counter clockwise order 
	 * around a center point center, as viewed in the direction opposite
	 * to the normal vector normal_vec.
	 * In 2D, normal_vec is not used. */
	class ccw_comparator {
	public:
		typedef unsigned int index_type;
		ccw_comparator(const vertice_vector& vertices,
				const vertice_type& center, const vertice_type& normal_vec) :
			my_vertices(vertices), my_center(center) {
			if (my_center.size() == 3) {
				my_viewpoint = center + normal_vec;
			}
		}
		;
		bool operator()(index_type x_index, index_type y_index) {
			return operator()(my_vertices[x_index], my_vertices[y_index]);
		}
		;
		bool operator()(const vertice_type& x, const vertice_type& y) {
			assert(x.get_index_to_variable_id_map()
					== y.get_index_to_variable_id_map());
			assert(x.size() == 2 || x.size() == 3);

			bool res = false;
			if (x.size() == 2) {
				res = is_counter_clockwise_2D(x.get_vector(), y.get_vector(),
						my_center.get_vector());
			} else {
				res = is_counter_clockwise_3D(x.get_vector(), y.get_vector(),
						my_center.get_vector(), my_viewpoint.get_vector());
			};
			return res;
		}
		;
	private:
		const vertice_vector& my_vertices;
		const vertice_type& my_center;
		vertice_type my_viewpoint;
	};

	class vector_comparator {
	public:
		bool operator()(const vertice_type& x, const vertice_type& y) {
			return x.get_vector() < y.get_vector();
		}
	};

	/** Sort indices such that they point to the vertices in counter clockwise order with respect to the
	 * normal vector normal_vec. */
	static void sort(const vertice_vector& vertices, index_list& indices,
			const vertice_type& normal_vec) {
		if (!vertices.empty()) {
			/* compute the center */
			vertice_type center = compute_average(vertices, indices);
			indices.sort(ccw_comparator(vertices, center, normal_vec));
		}
	}
	;

	/** Sort vertices such that they point in counter clockwise order with respect to the
	 * normal vector normal_vec. */
	static void sort(vertice_vector& vertices, const vertice_type& normal_vec) {
		if (!vertices.empty()) {
			/* compute the center */
			vertice_type center = compute_average(vertices);

			//			std::sort(vertices.begin(), vertices.end(),vector_comparator());
			std::sort(vertices.begin(), vertices.end(), ccw_comparator(
					vertices, center, normal_vec));

			//			vertice_vector old_v;
			//			while (old_v != vertices) {
			//				old_v = vertices;
			//				//std::cout << center << std::endl;
			//				//std::cout << "before " << vertices << std::endl;
			//				std::stable_sort(vertices.begin(), vertices.end(),
			//						ccw_comparator(vertices, center, normal_vec));
			//			}
			//			//			std::cout << std::endl << "again " << vertices;
		}
	}
	;

	class polar_comparator {
	public:
		typedef unsigned int index_type;
		polar_comparator(const vertice_vector& vertices,
				const vertice_type& center) :
			my_vertices(vertices), my_center(center) {
		}
		;
		polar_comparator(const vertice_vector& vertices,
				const vertice_type& center, const vertice_type& start_vec,
				const vertice_type& normal_vec) :
			my_vertices(vertices), my_center(center) {
			my_startc=start_vec - my_center;
			my_normal=&normal_vec;
		}
		;

		bool operator()(index_type x_index, index_type y_index) {
			return operator()(my_vertices[x_index], my_vertices[y_index]);
		}
		;
		bool operator()(const vertice_type& x, const vertice_type& y) {
			assert(x.get_index_to_variable_id_map()
					== y.get_index_to_variable_id_map());
			//assert(x.size() == 2 && y.size() == 2);

			vertice_type xc = x - my_center;
			vertice_type yc = y - my_center;

			double angle_x;
			double angle_y;
			if (x.size() == 2) {
				angle_x = std::atan2(convert_element<double> (xc[1]),
						convert_element<double> (xc[0]));
				angle_y = std::atan2(convert_element<double> (yc[1]),
						convert_element<double> (yc[0]));
			} else { //if (x.size() == 3) {
				angle_x = angle_3D(xc.get_vector(), my_startc.get_vector(),
						my_normal->get_vector());
				angle_y = angle_3D(yc.get_vector(), my_startc.get_vector(),
						my_normal->get_vector());
			}
			return angle_x < angle_y;
		}
		;
	private:
		const vertice_vector& my_vertices;
		const vertice_type& my_center;
		const vertice_type* my_normal;
		vertice_type my_startc;
	};

	/** Sort 2D vertices such that they point in counter clockwise order
	 * based on their angle w.r.t. the point (1,0). */
	static void angle_sort_2D(vertice_vector& vertices) {
		/* compute the center */
		vertice_type center = compute_average(vertices);

		index_list indices;
		for (unsigned int i = 0; i < vertices.size(); ++i) {
			indices.push_back(i);
		}

		indices.sort(polar_comparator(vertices, center));

		vertice_vector new_vertices(vertices.size());
		index_list::const_iterator it = indices.begin();
		for (unsigned int i = 0; i < vertices.size(); ++i) {
			new_vertices[i] = vertices[*it];
			++it;
		}
		vertices = new_vertices;
	}
	;

	/** Sort 2D vertices such that they point in counter clockwise order
	 * based on their angle w.r.t. to the first point in the set
	 * and defining normal_vec as outside. */
	static void angle_sort_3D(vertice_vector& vertices, index_list& indices,
			const vertice_type& normal_vec) {
		/* compute the center */
		vertice_type center = compute_average(vertices);

		indices.sort(polar_comparator(vertices, center, *vertices.begin(),
				normal_vec));
	}
	;

	/** Sort 2D vertices such that they point in counter clockwise order
	 * based on their angle w.r.t. to the first point in the set
	 * and defining normal_vec as outside. */
	static void angle_sort_3D(vertice_vector& vertices,
			const vertice_type& normal_vec) {
		typename math::vdom_vector<scalar_type> dummy;

		index_list indices;
		for (unsigned int i = 0; i < vertices.size(); ++i) {
			indices.push_back(i);
		}

		angle_sort_3D(vertices, indices, normal_vec);

		vertice_vector new_vertices(vertices.size());
		index_list::const_iterator it = indices.begin();
		for (unsigned int i = 0; i < vertices.size(); ++i) {
			new_vertices[i] = vertices[*it];
			++it;
		}
		vertices = new_vertices;
	}
	;
};

/** Returns the elementwise min of two vectors.
 *
 * Throws if the vectors are not of the same dimension,
 * and if the two domains do not contain the same variables.
 */
template<typename scalar_type> vdom_vector<scalar_type> min(const vdom_vector<scalar_type>& v1, const vdom_vector<scalar_type>& v2) {
	if (v1.get_index_to_variable_id_map() == v2.get_index_to_variable_id_map()) {
		return vdom_vector<scalar_type> (min(v1.get_vector(),v2.get_vector()),
				v1.get_index_to_variable_id_map());
	} else {
		index_to_variable_id_map_ptr new_iimap;
		index_type new_dim;
		index_to_index_bimap remap2;
		get_common_map(v1.get_index_to_variable_id_map(),
				v2.get_index_to_variable_id_map(), new_iimap, new_dim, remap2);
		if (new_dim>v1.size() || new_dim>v2.size())
			throw std::runtime_error("Cannot compute min of vectors with different variables");

		typename vdom_vector<scalar_type>::vector_type new_v(new_dim);
		vdom_vector<scalar_type> new_l(new_iimap);
		for (unsigned int i = 0; i < v1.size(); ++i) {
			new_l[i] = v1[i];
		}
		for (unsigned int i = 0; i < v2.size(); ++i) {
			if (v2[i] < new_l[remap2.get_map(i)])
				new_l[remap2.get_map(i)] = v2[i];
		}
		return new_l;
	}
}

}

#endif /*vdom_VECTOR_UTILITY_H_*/
