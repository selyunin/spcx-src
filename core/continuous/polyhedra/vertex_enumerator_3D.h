#ifndef VERTEX_ENUMERTOR_3D_H_
#define VERTEX_ENUMERTOR_3D_H_

/** This is a 3D vertex enumeration algorithm that
 * creates a list of vertices and faces (defining incidence with vertices).
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

template<typename scalar_type> class vertex_enumerator_3D {
public:
	typedef unsigned int con_id;
	typedef std::set<con_id> con_id_set;
	typedef unsigned int vert_id;
	typedef std::list<vert_id> vert_id_list;
	typedef std::vector<con_id_set> con_id_set_vector;
	typedef vector<scalar_type> vert_type;
	typedef std::map<vert_type, vert_id, numeric::lex_comp_less<scalar_type,
			vector> > vert_id_map;

	typedef vdom_vector<scalar_type> vertice_type;
	typedef std::vector<vertice_type> vertice_vector;
	typedef std::vector<vert_id_list> face_vector;

	vertex_enumerator_3D(const lin_constraint_system<scalar_type>& cons,
			vertice_vector& vertices, face_vector& faces) :
		my_faces(faces) {
		if (cons.size() > 0) {
			canonic_matrix_form(my_A, my_b, my_domain, cons);

			my_dim = my_domain.size();
			if (my_dim != 3)
				throw basic_exception("vertex_enumerator_3D: can not handle "
						+ to_string(my_dim)
						+ " variables, only 3 variables are allowed.");
		}

		if (cons.size() > 0 && my_dim == 3) {
//std::cout << "A:" << my_A << std::endl;
//std::cout << "b:" << my_b << std::endl;

			my_faces = face_vector(my_A.size1());
			Edeg = con_id_set_vector(my_A.size1());
			try {
				enumerate_vertices();
			} catch ( std::exception& e ) {
				std::stringstream ss;
				ss << cons;
				throw basic_exception("Can not enumerate vertices for polyhedron:"+ss.str(),e);
			}

			/* convert the result back to a vertice list */
			vertices = vertice_vector(my_vertmap.size());
			unsigned int pos_counter = 0;
			for (typename vert_id_map::const_iterator it = my_vertmap.begin(); it
					!= my_vertmap.end(); ++it) {
				vertices[it->second] = vertice_type(my_domain, it->first);
				++pos_counter;
			}
			/* the face incidence lists are already in my_faces */
			// if vertices are not in ccw order, reverse
			for (unsigned int r = 0; r < my_A.size1(); ++r) {
				if (my_faces[r].size() >= 3) {
					vert_id_list::const_iterator it = my_faces[r].begin();
					const vert_type& v1 = vertices[*it].get_vector();
					++it;
					const vert_type& v2 = vertices[*it].get_vector();
					++it;
					const vert_type& v3 = vertices[*it].get_vector();
					vert_type norm_v=my_A.vector_from_row(r);
				if (!is_counter_clockwise_3D(v1,v2,v3,v2+norm_v)) {
					my_faces[r].reverse();
				}
			}
		}
		} else {
			vertices = vertice_vector();
			faces = face_vector();
		}
	}
	;

	virtual ~vertex_enumerator_3D() {
	}
	;

	/** Compute vertices and faces using different number types for computation and final result.
	 */
	template<typename output_type>
	static void compute(
			const lin_constraint_system<output_type>& cons,
			typename vertex_enumerator_3D<output_type>::vertice_vector& vertices,
			typename vertex_enumerator_3D<output_type>::face_vector& faces) {
		lin_constraint_system<scalar_type> scons = cons.template convert_to<
				scalar_type> ();
		typename vertex_enumerator_3D<scalar_type>::vertice_vector svertices;
		vertex_enumerator_3D<scalar_type>(scons, svertices, faces);

		// convert the result vertices to output type
		vertices.resize(svertices.size());
		for (unsigned int i=0;i<svertices.size();++i)
			vertices[i]=svertices[i].template convert_to<output_type> ();
	}
	;

	vert_id add_and_get_vertice_id(const vert_type& v) {
		typename vert_id_map::iterator it = my_vertmap.lower_bound(v); // first element GE
		typename vert_id_map::iterator jt = my_vertmap.upper_bound(v); // first element GT
		//std::cout << std::setprecision(16) << "searching " << v << std::endl;
		if (it == my_vertmap.end() || it == jt) {
			vert_id new_id = my_vertmap.size();
			it = my_vertmap.insert(it, std::make_pair(v, new_id));
			//std::cout << std::setprecision(16) << "inserting new " << v << " with id " << new_id << std::endl;
		}
		return it->second;
	}
	;

	/** Add vertice v to face r */
	void add_vertice(const vert_type& v, unsigned int r) {
		vert_id vid = add_and_get_vertice_id(v);
		my_faces[r].push_back(vid);
	}
	;

	unsigned int largest_nonzero_coeff(const matrix<scalar_type>& A, unsigned int r) {
		using numeric::is_GT;
		unsigned int i = A.size2();
		scalar_type max_coeff = scalar_type(0);
		for (unsigned int ii = 0; ii < 3; ++ii) {
			if (is_GT(A(r, ii),max_coeff)) {
				max_coeff = A(r, ii);
				i = ii;
			} else if (is_GT(-A(r, ii),max_coeff)) {
				max_coeff = -A(r, ii);
				i = ii;
			}
		}
		return i;
	}

	void enumerate_vertices() {
		assert(my_dim==3);
		unsigned int m = my_A.size1();
		unsigned int n = my_A.size2();
		unsigned int i, j, k; // indices of pivot columns
		unsigned int r, s, t; // indices of faces
		unsigned int tl, tu; // indices of the active faces
		scalar_type cj, ck, d; // coefficients of the first pivot column
		unsigned int j1, k1; // indices of pivot columns of the first pivot
		scalar_type ek, f, g, h; // coefficients for the third pivot
		scalar_type lb, ub;
		using numeric::is_MEQ;
		using numeric::is_LT;
		using numeric::is_LE;
		using numeric::is_GE;
		using numeric::is_GT;
		scalar_type max_coeff;
		my_empty = false;
		my_unbounded = false;
		const scalar_type zero(0);
		vert_type vl(3);
		vert_type vu(3);
		vert_type vfirst(3); // first vertice added to face r
		vert_type vlast(3); // last vertice added to face r

		// store one nondegenerate Edges for each face
		Enondeg=vector<unsigned int>(m,m); // default value m:unknown

		// Iterate over faces
		for (r = 0; r < m && !my_empty && !my_unbounded; ++r) {
			unsigned int r_edgecount = 0;
			unsigned int r_verticecount = 0;
			// last chosen edge
			unsigned int last_s = m;

			matrix<scalar_type> A=my_A;
			vector<scalar_type> b=my_b;
			// find a nonzero coefficient, use the largest one
			// that is GT zero (otherwise it's considered null)
			i = largest_nonzero_coeff(A,r);
			if (i >= n) {
				// didn't find a nonzero coeff
				// test if satisfiable
				if (is_LT(b[r], zero)) {
					my_empty = true;
				}
			} else {
//std::cout << std::endl << "Pivoting on " << r << "," << i << std::endl;
				// found a valid face
				// define j,k as the other two indices
				j = 0;
				while (j==i && j < n) {
					++j;
				}
				k = 0;
				while ((k==i || k==j) && k < n) {
					++k;
				}
				// compute pivot coeffs
				cj = -A(r, j) / A(r, i);
				ck = -A(r, k) / A(r, i);
				d = b[r] / A(r, i);

				// we need to remember the columns
				// because we'll reassign j and k,
				// but use cj and ck later
				j1=j;
				k1=k;
//std::cout << "cj " << cj << " ck " << ck << " d " << d << std::endl;

				// apply pivot
				for (s=0;s<m;++s) {
					if (s!=r) {
						A(s, j) += A(s, i) * cj;
						A(s, k) += A(s, i) * ck;
						b[s] -= A(s, i) * d;
						A(s, i) = zero;
					}
				}
//std::cout << A << std::endl;
//std::cout << b << std::endl;

				// find a first edge and pivot index
				// find s!=r with nonzero coeff
				if (Enondeg[r]<m) {
					s = Enondeg[r];
				} else {
					// start looking from zero
					s = 0;
				}
				// find pivot j
				if (s < m)
					j = largest_nonzero_coeff(A,s);
				while (s < m && (s == r || j >= n || Edeg[r].find(s)!=Edeg[r].end())) {
					// if nothing found, go to the next one
					++s;
					if (s < m)
						j = largest_nonzero_coeff(A,s);
				}
				while (s < m && !my_unbounded) {
					++r_edgecount;
					// find pivot j
					j = largest_nonzero_coeff(A,s);
					// k is the remaining pivot
					k = 0;
					while ((k==i || k==j) && k < n) {
						++k;
					}
					//std::cout << std::endl << "Second pivot on " << s << "," << j << std::endl;
					// found a valid edge s
					// compute pivot coeffs

					ek = -A(s,k) / A(s,j);
					f = b[s] / A(s,j);
//std::cout << "ek " << ek << " f " << f << std::endl;
					bool has_lb = false;
					bool has_ub = false;
					bool infeas = false;
					// mark lower and upper bound indices as nonexisting
					tl = m;
					tu = m;
					for (t = 0; t < m && !infeas; ++t) {
						if (t != r && t != s) { // skip r and s
							g = A(t,k) + A(t,j) * ek;
							h = b[t] - A(t,j) * f;
////std::cout << "bound " << t << " is " << h << "/" << g << std::endl;
							if (is_LT(g, zero)) {
								// lower bound
								if (!has_lb || is_GT(h / g, lb)) {
									has_lb = true;
//std::cout << "tightening lb = "<<h / g << " with face " << t << std::endl;
									lb = h / g;
									tl = m;
								}
								if (tl == m && t != last_s && is_MEQ(h / g, lb)
										&& Edeg[r].find(t) == Edeg[r].end()) {
									tl = t;
//std::cout << "lb candidate edge "<< t << std::endl;
								}
							} else if (is_GT(g, zero)) {
								// upper bound
								if (!has_ub || is_LT(h / g, ub)) {
									has_ub = true;
//std::cout << "tightening ub = "<<h / g << " with face " << t << std::endl;
									ub = h / g;
									tu = m;
								}
								if (tu == m && t != last_s && is_MEQ(h / g, ub)
										&& Edeg[r].find(t) == Edeg[r].end()) {
									tu = t;
//std::cout << "ub candidate edge "<< t << std::endl;
								}
							} else {
								// g is zero, so the constraint is infeasible
								// if h<0
								if (is_LT(h, zero))
									infeas=true;
							}
							if (has_lb && has_ub && is_GT(lb,ub)) {
								infeas=true;
							}
						}
					} // for t
					if (!infeas) {
						if (!has_lb || !has_ub) {
							my_unbounded = true;
							throw std::runtime_error(
									"vertex_enumerator_3D called with unbounded polyhedron");
						} else {
							if (is_MEQ(lb, ub)) {
							// the edge (r,s) is degenerate
//std::cout << "found degen edge " << "(" << r << "," << s << ")" << std::endl;
							Edeg[r].insert(s);
							Edeg[s].insert(r);
						}
						// update vertices
//td::cout << "lb=" << lb << ",k=" << k <<",j="<<j<<",i="<<i<<",f="<<f<<",ek="<<ek<<",d="<<d<<",cj="<<cj<<",ck="<<ck<<std::endl;
						vl[k] = lb;
							vl[j] = f + ek * vl[k];
							vl[i] = d + cj * vl[j1] + ck * vl[k1];
							vu[k] = ub;
							vu[j] = f + ek * vu[k];
							vu[i] = d + cj * vu[j1] + ck * vu[k1];
//std::cout << "vl:" << vl << ", vu:" << vu << std::endl << std::flush;
							// mark new s as end of constraints
							last_s = s;
							s = m;
							if (r_verticecount == 0) {
								// try to pick one that's not degenerate
//std::cout << "adding FIRST vl:" << vl << std::endl << std::flush;
								add_vertice(vl, r);
								++r_verticecount;
								if (!is_MEQ(lb,ub)) {
//std::cout << "adding FIRST vu:" << vu << std::endl << std::flush;
									add_vertice(vu, r);
									++r_verticecount;
								}
								if (tl < m) {
//std::cout << " next s from vl:" << tl << std::endl << std::flush;
									s = tl;
									vlast = vl; // this will corrupt the comparison vlast,vu if we don't use take_vu
									vfirst = vu;
									my_faces[r].reverse(); // should have added vl first
								} else {
//std::cout << "  next s from vu:" << tu << std::endl << std::flush;
									s = tu;
									vlast = vu;
									vfirst = vl;
								} // a
							} else {
								bool take_vl = !is_MEQ(vlast, vl);
								bool take_vu = !is_MEQ(vlast, vu);
								// if due to numerical problems both are different,
								// choose the farthest one
								if (take_vl && take_vu) {
									scalar_type dl = norm2_square(vlast-vl);
									scalar_type du = norm2_square(vlast-vu);
//std::cout << "progress in both directions, dl=" << dl << ", du=" << du << std::endl;
									if (du > dl) {
										take_vl=false;
									} else {
										take_vu=false;
									}
								}
								if (take_vl) {
									vlast = vl; // this will corrupt the comparison vlast,vu if we don't use take_vu
									s = tl;
//std::cout << "progress on lb edge " << s << std::endl << std::flush;
								}
								if (take_vu) {
									vlast = vu; // this will corrupt the comparison vlast,vu if we don't use take_vu
									s = tu;
//std::cout << "progress on ub edge " << s << std::endl << std::flush;
								}
								bool same_as_first = is_MEQ(vfirst, vlast);
								if (s < m && !same_as_first) {
									// found a vertex, so add it
//std::cout << "adding vertice:" << vlast << std::endl << std::flush;
									add_vertice(vlast, r);
									++r_verticecount;

									// add the edge as nondeg to s
									if (Enondeg[s] >= m)
										Enondeg[s] = r;
								} else {
									// didn't find a vertex, so try tl or tu if <m
									if (tl < m)
										s = tl;
									else
										s = tu;
//std::cout << "no progress, trying pivot on " << s << std::endl << std::flush;
								}
								// if more than one vertice has been found
								// and we're back to the first one, stop
								if (r_verticecount > 1 && same_as_first) {
									s = m;
								}
							} // if r_verticecount == 0
//std::cout << " next s:" << s << std::endl << std::flush;
						}
					} else {
						// the edge (r,s) is infeasible
						Edeg[r].insert(s);
						Edeg[s].insert(r);
						// look for another s
//std::cout << "infeasible: (" << r << "," << s << ")" << std::flush;
						// I think this can only happen for the first one
						++s;
						// find pivot j
						if (s < m)
							j = largest_nonzero_coeff(A,s);
						while (s < m && (s == r || j >= n || Edeg[r].find(s)!=Edeg[r].end())) {
							// if nothing found, go to the next one
							++s;
							if (s < m)
								j = largest_nonzero_coeff(A,s);
						}
//std::cout << "   found " << s << std::endl << std::flush;
					}
////std::cout << "degen/infeas for " << s << ": "<< Edeg[s] << std::endl << std::flush;
					if (r_edgecount >= m) {
						// There has to have been an error, since there cannot be more
						// edges than faces

						// Fix vertice list for this face by stopping when a double was hit
						if (my_faces[r].size() > 1) {
							vert_id_list::iterator it = my_faces[r].begin();
							++it;
							bool found=false;
							while (it != my_faces[r].end() && !found) {
								// check to see if *it exists in the range [begin(),it)
								// if yes, remove the rest
								found = (find(my_faces[r].begin(),it,*it)!=my_faces[r].end());
								if (!found) {
									++it;
								}
							}
							if (found) {
								my_faces[r].erase(it,my_faces[r].end());
							}
						}

						//std::cout << "Warning: vertex enumeration stopped by cycle detection." << std::endl;

						// Stop exploration of this edge
						s = m;
					}
				} // while s
			}
		} // for r
	}
	;

private:
	matrix<scalar_type> my_A;
	vector<scalar_type> my_b;
	vert_id_map my_vertmap;
	face_vector& my_faces;
	con_id_set_vector Edeg;
	vector<unsigned int> Enondeg;
	unsigned int my_dim;
	positional_vdomain my_domain;
	bool my_empty;
	bool my_unbounded;
};

}

#endif /*VERTEX_ENUMERTOR_3D_H_*/
