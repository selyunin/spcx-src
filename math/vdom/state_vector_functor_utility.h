/*
 * state_vector_functor_utility.h
 *
 *  Created on: Jun 15, 2011
 *      Author: frehse
 */

#ifndef STATE_VECTOR_FUNCTOR_UTILITY_H_
#define STATE_VECTOR_FUNCTOR_UTILITY_H_

#include "state_vector_functor.h"
#include "state_functor.h"
#include "math/ode_solving/traj_simu/lin_constraint_evaluation.h"
#include "math/vdom/lin_constraint_system.h"
#include "math/vdom/vdom_vector_operators.h"

namespace math {

/** A functor that computes the scalar product of a state with a vector of fixed, given other states
 *
 * The result is returned as a vector of scalars. */
template<typename scalar_type>
class scalar_product_functor: public state_vector_functor<scalar_type> {
public:
	/** Member types */
	typedef typename state_vector_functor<scalar_type>::state state;
	typedef typename state_vector_functor<scalar_type>::vector vector;
	typedef math::state_functor<scalar_type> state_functor;
	typedef std::vector<state> state_vector;

	/** Construct a scalar product functor from a given state functor and
	 * a vector of states that will be used for computing the scalar products.
	 */
	scalar_product_functor(const state_functor& f, const state_vector& W) :
		my_f(f), my_W(W) {
	}
	;

	/** Virtual desctructor for possible derived classes */
	virtual ~scalar_product_functor() {
	}
	;

	/** Map x to w^T f(x) for each element w of the state vector W given at construction time.
	 *
	 * The i-th element of the output vector corresponds to the
	 * scalar product of x with the i-th element of W.
	 */
	virtual vector map(const state& x) const {
		vector res(my_W.size());
		for (typename vector::size_type i = 0; i < my_W.size(); ++i) {
			state y = my_f.map(x);
			res[i] = scalar_product(my_W[i], y);
		}
		return res;
	}
	;
private:
	const state_functor& my_f;
	state_vector my_W;
};

/** An state vector functor implementation for affine maps. */
template<typename scalar_type>
class affine_state_vector_functor: public state_vector_functor<scalar_type> {
public:
	/** Member types */
	typedef typename state_vector_functor<scalar_type>::state state;
	typedef typename state_vector_functor<scalar_type>::vector vector;
	/** Construct a state functor corresponding to an affine map */
	affine_state_vector_functor(const affine_map<scalar_type>& M) :
		my_map(M) {
	}
	;

	/** Virtual destructor for possible derived classes */
	virtual ~affine_state_vector_functor() {
	}
	;

	/** Map the state x to another state according to the function
	 * defined by *this. */
	virtual vector map(const state& x) const {
		return my_map.map(x).get_vector();
	}
	;
private:
	const affine_map<scalar_type>& my_map;
};

/** An state vector functor implementation for vector of affine maps. */
template<typename scalar_type>
class affine_map_vector_state_vector_functor: public state_vector_functor<
		scalar_type> {
public:
	/** Member types */
	typedef typename state_vector_functor<scalar_type>::state state;
	typedef typename state_vector_functor<scalar_type>::vector vector;

	/** Construct a state functor corresponding to an affine map vector
	 *  the i-th int of the select_row parameter decides how the map is done
	 *  	- if it is <0, the whole map is computed.
	 *  	- if it is >=0, only the correspondig row is mapped
	 *  approx_size is an estimation of the size of the output (mapped) vector
	 * */
	affine_map_vector_state_vector_functor(
			const std::vector<typename affine_map<scalar_type>::const_ptr> & M,
			const std::vector<int> & select_row, unsigned int approx_size) :
		my_map(M), my_select_row(select_row), my_approx_size(approx_size) {
	}
	;

	/** Virtual destructor for possible derived classes */
	virtual ~affine_map_vector_state_vector_functor() {
	}
	;

	/** Map the state x to another state according to the function
	 * defined by *this. */
	virtual vector map(const state& x) const {
		vector res(my_approx_size);
		unsigned int pos = 0;

		for (int i = 0; i < my_map.size(); i++) {
			if (my_select_row[i] >= 0) {
				state row_i = my_map[i]->get_A().vdom_vector_from_row(
						my_select_row[i]);
				scalar_type temp = scalar_product(row_i, x)
						+ my_map[i]->get_b()[i];
				if (pos + 1 > res.size())
					res.resize(pos + 1);
				res[pos] = temp;
				pos++;
			} else {
				vector temp = my_map[i]->map(x).get_vector();
				if (pos + temp.size() > res.size())
					res.resize(pos + temp.size());
				res.subvector_assign(temp, pos, pos + temp.size());
				pos += temp.size();
			}
		}
		if (pos != res.size())
			res.resize(pos);
		return res;
	}
	;

	virtual void update_approx_size(unsigned int size) {
		this->my_approx_size = size;
	}
	;

private:
	const std::vector<typename affine_map<scalar_type>::const_ptr> & my_map;
	const std::vector<int> & my_select_row;
	unsigned int my_approx_size;
};

/**Class for defining linear constraint state functor.
 * Assumes that the linear constraints are given in normal form
 *WARNING : does not copy the constraint system, but
 *WARNING : uses reference (be careful with stack)
 */
template<typename scalar_type>
class lin_constraint_state_functor: public state_vector_functor<scalar_type> {
public:

	typedef typename state_vector_functor<scalar_type>::state state;
	typedef typename state_vector_functor<scalar_type>::vector vector;

	typedef struct lin_constraints {
		typename math::lin_constraint_system<scalar_type>::const_iterator begin;
		typename math::lin_constraint_system<scalar_type>::const_iterator end;
	} lin_constraints;

	typedef std::vector<lin_constraints> lin_constraints_vector;

	lin_constraint_state_functor(const lin_constraints_vector &M,
			unsigned int approx_size = 0) :
		my_map(M) {
		this->my_approx_size = approx_size;
		my_premap = NULL;
	}
	;

	lin_constraint_state_functor(const lin_constraints_vector &M,
			const state_functor<scalar_type> * premap,
			unsigned int approx_size = 0) :
		my_map(M) {
		this->my_approx_size = approx_size;
		my_premap = premap;
	}
	;

	virtual ~lin_constraint_state_functor() {
	}
	;

	virtual vector map(const state& x) const {
		if (my_premap)
			return duplicated_with_premap_map(x);
		else
			return simple_map(x);
	}

	virtual vector simple_map(const state& x) const {
		/*to implement*/

		vector temp(my_approx_size);

		typename lin_constraints_vector::const_iterator lit;
		typename math::lin_constraint_system<scalar_type>::const_iterator sit;
		int pos = 0;

		for (lit = my_map.begin(); lit != my_map.end(); ++lit) {
			for (sit = lit->begin; sit != lit->end; ++sit) {

				if (pos + 1 > temp.size())
					temp.resize(2 * temp.size());
				temp[pos] = math::lin_constraint_evaluation<scalar_type>(*sit,
						x);
				pos++;

			}
		}

		temp.resize(pos);

		//				resul.set_vector(temp);
		/*
		 printf("IN ODE_SOLVER::MAP ! %i values mapped : res=(",pos);
		 for(int h=0;h<pos;h++)
		 {
		 printf("%f,",temp[h]);
		 }
		 printf(")\n");
		 */
		return temp;
	}
	;

	virtual vector duplicated_with_premap_map(const state& x) const {
		/*to implement*/

		vector temp(my_approx_size);

		state xprim(my_premap->map(x));
		typename lin_constraints_vector::const_iterator lit;
		typename math::lin_constraint_system<scalar_type>::const_iterator sit;
		int pos = 0;

		for (lit = my_map.begin(); lit != my_map.end(); ++lit) {
			for (sit = lit->begin; sit != lit->end; ++sit) {

				if (pos + 1 > temp.size())
					temp.resize(2 * temp.size());
				temp[pos] = math::lin_constraint_evaluation<scalar_type>(*sit,
						x);
				temp[pos + 1] = math::lin_constraint_evaluation<scalar_type>(
						*sit, xprim);
				pos += 2;

			}
		}

		temp.resize(pos);

		//				resul.set_vector(temp);
		/*
		 printf("IN ODE_SOLVER::MAP ! %i values mapped : res=(",pos);
		 for(int h=0;h<pos;h++)
		 {
		 printf("%f,",temp[h]);
		 }
		 printf(")\n");
		 */
		return temp;
	}
	;

	virtual void update_approx_size(unsigned int size) {
		this->my_approx_size = size;
	}
	;

private:
	const lin_constraints_vector & my_map;
	unsigned int my_approx_size;
	const state_functor<scalar_type> * my_premap;
};

}

#endif /* STATE_VECTOR_FUNCTOR_UTILITY_H_ */
