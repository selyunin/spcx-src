#ifndef VERTEX_ENUMERTOR_BF_H_
#define VERTEX_ENUMERTOR_BF_H_

/** This is a slow vertex enumeration algorithm that
 * creates a list of vertices and faces (defining incidence with vertices).
 * It explores each constraint by brute force and tries to find sets of
 * constraints that have feasible point that lies on all constraints
 * (and therefore is the only feasible point modulo numerical erros).
 */

#include "math/numeric/approx_comparator.h"
#include "math/numeric/container_comp.h"
#include "math/vdom/vdom_vector.h"
#include "math/vdom/lin_constraint_system.h"
#include "math/lp_solving/lp_solver_user.h"
#include "math/lp_solving/lp_solver.h"
#include "math/matrix.h"
#include "math/vdom/lin_constraint_system_utility.h"

namespace math {

template <typename scalar_type> class vertex_enumerator_bf : public lp_solver_user<scalar_type> {
public:
	typedef vdom_vector<scalar_type> vertice_type;
	typedef unsigned int con_id;
	typedef std::set<con_id> con_id_set;
	typedef std::pair<const vertice_type, con_id_set>* vertice_id_type;
	typedef std::map<vertice_type, con_id_set, numeric::lex_comp_less<scalar_type,vdom_vector> > vertice_incidence_list;
	typedef std::vector<vertice_type> vertice_vector;
	typedef std::vector<std::list<unsigned int> > face_vector;
	typedef lin_constraint<scalar_type>* constraint_ref;
	typedef typename lin_constraint<scalar_type>::sign sign_type;
	typedef std::list<constraint_ref> constraint_ref_list;

	vertex_enumerator_bf(const lin_constraint_system<scalar_type>& cons, vertice_vector& vertices,
			face_vector& faces) :
		my_faces(faces) {
		if (cons.size() > 0) {
			my_cons.push_back(cons);
			my_cons.unify_domains();
			my_domain = my_cons.begin()->get_l().domain();
			my_dim = my_domain.size();
			//my_dim=cons.get_variable_ids().size();
					my_faces = face_vector(cons.size());
			if (my_dim > 0) {
				enumerate_vertices_bf(my_cons.begin(), my_cons.end(), 0);
			}
			/* convert the result back to a vertice list */
			vertices = vertice_vector(my_vertices.size());
			unsigned int pos_counter = 0;
			for (typename vertice_incidence_list::const_iterator it =
					my_vertices.begin(); it != my_vertices.end(); ++it) {
				vertices[pos_counter] = it->first;
				for (con_id_set::const_iterator jt = it->second.begin(); jt
						!= it->second.end(); ++jt) {
					my_faces[*jt].push_back(pos_counter);
				}
				++pos_counter;
			}
		} else {
			vertices=vertice_vector();
			faces=face_vector();
		}
	}
	;

	virtual ~vertex_enumerator_bf() {
	}
	;

	/** Check whether every constraint is tight on the point l. */
	bool is_on_all_constraints(const constraint_ref_list& cons,
			const vertice_type& l) {
//std::cout << "constraints: " << cons << std::endl;
		for (typename constraint_ref_list::const_iterator it=cons.begin(); it!=cons.end(); ++it) {
//std::cout << l << " on " << **it << " : " << 	(*it)->is_active(l) << std::endl;
			if (!(*it)->is_active(l))
				return false;
		}
		return true;
	}
	;

	/** Return the active constraints as A*x<=b */
	std::pair<matrix<scalar_type>,vector<scalar_type> > active_matrix() const {
		matrix<scalar_type> A(active_cons.size(),my_dim);
		vector<scalar_type> b(active_cons.size());
		unsigned int i = 0;

		for (typename constraint_ref_list::const_iterator it=active_cons.begin(); it!=active_cons.end(); ++it) {
			lin_expression<scalar_type> l=(*it)->get_canonic_l();
			for (unsigned int j=0;j<my_dim;++j) {
				A(i,j)=l[j];
			}
			b[i]=-(*it)->get_canonic_inh_coeff();
			++i;
		}
		return std::make_pair(A,b);
	};

	void test_for_vertex(vertice_type& v, bool& is_unsat, bool& is_vertex) {
		assert(active_cons.size()>0);

		is_vertex = false;
		if (active_cons.size() == my_dim) {
//std::cout << std::endl << "testing active cons " << active_cons << std::endl;
			std::pair<matrix<scalar_type> , vector<scalar_type> > Ab =
					active_matrix();
			bool singular;

			//check for singular matrix before calling inverse.

			if(math::numeric::is_MEQ<scalar_type>(math::matrix_determinant(Ab.first),scalar_type(0))){
				is_vertex = false;
				singular = true;
				return;
			}
			//

			matrix<scalar_type> Ainv = Ab.first.inverse(singular);
			if (!singular) {
				vector<scalar_type> vec = Ainv * Ab.second;
				v=vertice_type(my_domain,vec);
//std::cout << "potential vertice "<< v << "?";
				// test if it satisfies all constraints
				is_vertex = true;
				for (typename lin_constraint_system<scalar_type>::const_iterator it =
						my_cons.begin(); is_vertex && it != my_cons.end(); ++it) {
					if (!math::maybe(it->is_satisfied(v))) {
						is_vertex = false;
					}
				}
//std::cout << is_vertex << std::endl;
			} else {
				is_vertex = false;
			}
		}
//std::cout << "vertex" << is_vertex << std::endl;
	}
	;

	void test_for_vertex_lp(vertice_type& v, bool& is_unsat, bool& is_vertex) {
		assert(active_cons.size()>0);
		/* Choose an arbitrary nonzero cost function;
		 * here: first variable. */
		is_vertex=false;
		lin_expression<scalar_type> l=lin_expression<scalar_type>((*active_cons.begin())->get_l().get_index_to_variable_id_map());
		l[0]=scalar_type(1);
		typename lp_solver<scalar_type>::lp_result res;
		this->get_lp_solver()->maximize(l, my_cons, res);
//std::cout << std::endl << "testing active cons " << active_cons << std::endl;
//std::cout << std::endl << "testing cons " << my_cons << std::endl;
//std::cout << "obj " << l << std::endl;
//print<scalar_type>(std::cout,res); std::cout << std::endl;
		is_unsat=res.is_unsat;
		if (!is_unsat && res.is_bounded) {
			if (active_cons.size()>=my_dim) {
				is_vertex=is_on_all_constraints(active_cons, res.support_vec);
				if (is_vertex) v=res.support_vec;
			}
		}
//std::cout << "unsat" << is_unsat << "vertex" << is_vertex << std::endl;
	}
	;

	vertice_id_type add_and_get_vertice_id(const vertice_type& l) {
		// @todo: check for duplicates
		/* note that for some reason the examples don't show any duplicates
		 * (at least for non-redundant constraints) */
		typename vertice_incidence_list::iterator it=my_vertices.lower_bound(l); // first element GE
		typename vertice_incidence_list::iterator jt=my_vertices.upper_bound(l); // first element GT
//std::cout << std::setprecision(16) << "searching "<< l << std::endl;
		if (it==my_vertices.end() || it==jt) {
//std::cout << std::setprecision(16) << "inserting new "<< l << std::endl;
			it=my_vertices.insert(it,make_pair(l,con_id_set()));
		}
		return &(*it);
	};

	void add_vertice(const vertice_type& l) {
		vertice_id_type vid=add_and_get_vertice_id(l);
		/* add the active faces to the vertice */
		for (std::list<unsigned int>::const_iterator cit=active_con_ids.begin(); cit
				!=active_con_ids.end(); ++cit) {
			vid->second.insert(*cit);
		}
	};

	void enumerate_vertices_bf(const typename lin_constraint_system<scalar_type>::iterator& beg_it,
			const typename lin_constraint_system<scalar_type>::iterator& end_it, unsigned int con_nr) {
		for (typename lin_constraint_system<scalar_type>::iterator it=beg_it; it!=end_it;) {
			// @todo: don't add linear dependent constraints (replace redundant ones)
//cout << "active constraint : " << *it << std::endl;
			sign_type old_sign=it->get_sign();
			it->set_sign(EQ);
			active_cons.push_back(&(*it));
			active_con_ids.push_back(con_nr);
			bool is_unsat;
			bool is_vertex;
			vertice_type v;
			test_for_vertex(v,is_unsat, is_vertex);
			/* iterate to the next constraint because it needs to be passed to
			 * the next call anyway */
			if (!is_unsat) {
				if (is_vertex) {
					add_vertice(v);
				} else if (it!=end_it) {
					/* fix another constraint */
					//std::cout << "-------------branching : " <<*it << std::endl;
					typename lin_constraint_system<scalar_type>::iterator jt=it;
					++jt;
					enumerate_vertices_bf(jt,end_it,con_nr+1);
					//std::cout << "-------------returning : " <<*it << std::endl;
				}
			}
			active_con_ids.pop_back();
			active_cons.pop_back();
			it->set_sign(old_sign);
			++it;
			++con_nr;
		}
	}
	;

private:
	lin_constraint_system<scalar_type> my_cons;
	constraint_ref_list active_cons;
	std::list<unsigned int> active_con_ids;
	vertice_incidence_list my_vertices;
	face_vector& my_faces;
	unsigned int my_dim;
	positional_vdomain my_domain;
};

}

#endif /*VERTEX_ENUMERTOR_BF_H_*/
