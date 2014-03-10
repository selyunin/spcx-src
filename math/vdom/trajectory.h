/*
 * trajectory.h
 *
 *  Created on: Oct 7, 2010
 *      Author: frehse
 */

#ifndef TRAJECTORY_H_
#define TRAJECTORY_H_

#include "utility/tree_node.h"
#include "math/vdom/lin_expression.h"
#include "math/matrix.h"
//#include "math/vdom/index_to_variable_id_map_provider.h"
#include "math/vdom/positional_vdomain.h"
#include "math/vdom/vdom_matrix.h"
#include "math/vdom/vdom_matrix_utility.h"
#include "math/vdom/vdom_matrix_operators.h"
#include "math/vdom/vdom_vector.h"
#include "math/vdom/vdom_vector_utility.h"
#include "math/vdom/vdom_vector_operators.h"
#include "math/numeric/approx_comparator.h"

namespace math {

/** A class for representing a trajectory.
 *
 * A trajectory is a set of points in n-dimensional domain with
 * a time stamp for each point.
 *  */
template<typename scalar_type> class trajectory: public index_to_variable_id_map_provider,
		public printable {
public:
	typedef boost::shared_ptr<trajectory<scalar_type> > ptr;
	typedef boost::shared_ptr<const trajectory<scalar_type> > const_ptr;
	typedef matrix<scalar_type> matrix_type;
	typedef vector<scalar_type> vector_type;
	typedef vdom_matrix<scalar_type> vdom_matrix_type;
	typedef vdom_vector<scalar_type> vdom_vector_type;
	typedef typename matrix_type::size_type size_type;

	/** Define an empty trajectory on an empty domain. */
	trajectory() {
	}
	;

	/** Define an empty trajectory the domain dom. */
	trajectory(const positional_vdomain& dom) :
		index_to_variable_id_map_provider(dom.get_index_to_variable_id_map()) {
	}
	;

	/** Define a single state trajectory with state x at time t. */
	trajectory(scalar_type t, const vdom_vector_type& x) :
		index_to_variable_id_map_provider(
				x.domain().get_index_to_variable_id_map()), my_states(1,
				x.domain().size()), my_times(1, t) {
		my_states.assign_row(0, x.get_vector());
	}
	;

	/** Define a single state trajectory in domain dom with state x at time t. */
	trajectory(const positional_vdomain& dom, scalar_type t,
			const vector_type& x) :
		index_to_variable_id_map_provider(dom.get_index_to_variable_id_map()),
				my_states(1, dom.size()), my_times(1, t) {
		assert(dom.size()==x.size());
		my_states.assign_row(0, x);
	}
	;

	/** Define a trajectory, where X is a matrix whose rows are points
	 * in domain dom and tstamps is a vector with the same number of rows.
	 */
	trajectory(const positional_vdomain& dom, const vector_type& tstamps,
			const matrix_type& X) :
		index_to_variable_id_map_provider(dom.get_index_to_variable_id_map()) {
		assert(X.size1()==tstamps.size());
		assert(X.size2()==dom.size());
		my_states = X;
		my_times = tstamps;
	}
	;

	virtual ~trajectory() {
	}
	;

	virtual bool is_empty() const {
		return my_times.size() == 0;
	}

	/** Returns the number of time points in the trajectory. */
	virtual size_type size() const {
		return my_times.size();
	}

	/** Insert the state x with timestamp t.
	 *
	 * If the domain is empty, the domain of x is adopted.
	 * Otherwise, a differing domain triggers an exception.
	 *
	 * @todo for now this is attached at the back.
	 * maybe this should be sorted according to time.
	 * */
	virtual void insert(scalar_type t, const vector_type& x) {
		assert(domain().size()==x.size());
		unsigned int new_row = size();
		if (is_empty()) {
			my_times = vector_type(1);
			my_states = matrix_type(1, domain().size());
		} else {
			unsigned int new_size = size() + 1;
			my_times.resize(new_size);
			my_states.resize(new_size, domain().size());
		}
		my_times[new_row] = t;
		my_states.assign_row(new_row, x);
	}


	/** Insert the state x with timestamp t.
	 *
	 * If the domain is empty, the domain of x is adopted.
	 * Otherwise, a differing domain triggers an exception.
	 *
	 * @todo for now this is attached at the back.
	 * maybe this should be sorted according to time.
	 * */
	virtual void insert(scalar_type t, const vdom_vector_type& x) {
		if (domain().size()==0)
			set_domain(x.domain());
		else if (x.domain() != domain()) {
			throw std::runtime_error(
					"trying to insert state into trajectory of another domain");
		}
		insert(t,x.get_vector());
	}

	/** Insert the trajectory traj.
	 *
	 * @todo for now this is attached at the back.
	 * maybe this should be sorted according to time. */
	virtual void insert(const trajectory<scalar_type>& traj) {
		if (traj.size() > 0) {
			if (traj.domain() != domain()) {
				throw std::runtime_error(
						"trying to insert traj into trajectory of another domain");
			}
			size_type n = this->domain().size();
			unsigned int old_size = size();
			unsigned int new_size = size() + traj.size();
			my_times.resize(new_size);
			my_states.resize(new_size, domain().size());
			my_times.subvector_assign(traj.my_times, old_size, new_size);
			my_states.submatrix_assign(traj.my_states, old_size, new_size, 0, n);
		}
	}

	/** Insert the k-th state of traj.
	 *
	 * @todo for now this is attached at the back.
	 * maybe this should be sorted according to time. */
	virtual void insert(const trajectory<scalar_type>& traj, unsigned int k) {
		assert(k<traj.size());
		if (traj.size() > 0) {
			if (traj.domain() != domain()) {
				throw std::runtime_error(
						"trying to insert traj into trajectory of another domain");
			}
			size_type n = this->domain().size();
			unsigned int old_size = size();
			unsigned int new_size = size() + 1;
			unsigned int pos = old_size;
			my_times.resize(new_size);
			my_states.resize(new_size, domain().size());
			my_times[pos] = traj.my_times[k];
			for (unsigned int i=0;i<n;++i) {
				my_states(pos,i)=traj.my_states(k,i);
			}
		}
	}

	/**
	 * Converts the elements of the vector to the type result_type.
	 * Uses the template function convert_element, which can be specialized to
	 * implement specific conversions between types.
	 */
	template<typename result_type> trajectory<result_type> convert_to() const {
		typename trajectory<result_type>::matrix_type A =
				my_states.template convert_to<result_type> ();
		typename trajectory<result_type>::vector_type b =
				my_times.template convert_to<result_type> ();
		return trajectory<result_type> (domain(), b, A);
	}
	;

	/** Output the map in the form
	 * TIME name1 ... namen
	 * t1 x1(1) ... xn(1)
	 * ...
	 * tk x1(k) ... xn(k)
	 */
	virtual void print(std::ostream& os) const {
		assert(my_states.size1()==my_times.size());
		assert( my_states.size2()==0 || my_states.size2()==this->domain().size());

		size_type k = size();
		size_type n = this->domain().size();

		// Header
		os << "TIME ";
		for (unsigned int i = 0; i < n; ++i) {
			if (i > 0) {
				os << " ";
			}
			os << get_variable(i).get_name();
		}
		os << std::endl;
		// Time points
		for (unsigned int i = 0; i < k; ++i) {
			os << my_times[i];
			for (unsigned int j = 0; j < n; j++) {
				os << " ";
				os << my_states(i, j);
			}
			os << std::endl;
		}
	}
	;

	/** Returns the states as a matrix where each
	 * row corresponds to a state vector. */
	const matrix_type& get_states() const {
		return my_states;
	}
	;

	/** Returns the timestamps as a vector. */
	const vector_type& get_times() const {
		return my_times;
	}
	;

	/** Returns the state with index k. */
	vdom_vector_type get_state(unsigned int k) const {
		assert(k<size());
		vdom_vector_type v =
				vdom_vector_type(domain(), my_states.vector_from_row(k));
		return v;
	}
	;

	scalar_type get_time(unsigned int k) const{
		assert(k<size());
		return my_times[k];
	}

	/**
	 *  Vector Swapping
	 *  */
	void swap(trajectory<scalar_type>& v) {
		if (this != &v) {
			index_to_variable_id_map_provider::swap(v);
			my_states.swap(v.my_states);
			my_times.swap(v.my_times);
		}
	}
	;

	friend void swap(trajectory<scalar_type> &v1, trajectory<scalar_type> &v2) {
		v1.swap(v2);
	}
	;

protected:
	matrix_type my_states;
	vector_type my_times;
};

}

#endif /* TRAJECTORY_H_ */
