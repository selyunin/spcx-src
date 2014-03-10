#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include "boost/shared_ptr.hpp"
#include <stdexcept>
#include <typeinfo>

#include "math/scalar_types/rational.h"
#include "core/continuous/continuous_set_declarations.h"
#include "core/continuous/polyhedra/polyhedron.h"
#include "core/continuous/polyhedra/polyhedron_upcaster.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/continuous/polyhedra/polyhedron_collection.h"

#include "core/continuous/support_function/spacetime_flowpipe.h"

//for simulation
#include "math/ode_solving/traj_simu/continuous_set_simulation.h"



/** Forward declarations */
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;
}

namespace continuous {

template<typename T1, typename T2> class intersection_operator_after_polycast;
template<typename T1, typename T2> class intersection_upcaster;

/** Declaration of the intersection operator. */
template<typename T1, typename T2> class intersection_operator {
public:
	static continuous_set_ptr implement(const T1* t1, const T2* t2) {
		return intersection_operator_after_polycast<
				typename intersection_upcaster<T1, T2>::result1,
				typename intersection_upcaster<T1, T2>::result2>::implement(t1,
				t2);
	}
	;
};

/** If no specialization was found, try casting up. */
template<typename T1, typename T2> class intersection_operator_after_polycast {
public:
	static continuous_set_ptr implement(const T1* t1, const T2* t2) {
		std::string name1 = typeid(t1).name();
		std::string name2 = typeid(t2).name();
		std::string real_name1 = typeid(*t1).name();
		std::string real_name2 = typeid(*t2).name();
		throw std::runtime_error("intersection_operator not defined for "
				+ name1 + " and " + name2 + ". Actual types are " + real_name1 + " and " + real_name2 + ".");
		return continuous_set_ptr();
	}
	;
};

/** If T1 is a polyhedron, try to cast T2 to a polyhedron. */
template<typename T1, typename T2>
class intersection_upcaster {
public:
	typedef typename polyhedron_upcaster<T1>::result result1;
	typedef typename boost::mpl::if_c<polyhedron_upcaster<T1>::is_poly,
			typename polyhedron_upcaster<T2>::result, T2>::type result2;
};

/** Implementation for polyhedra of the same type */
template<typename T> class intersection_operator_after_polycast<polyhedron<T> , polyhedron<T> > {
public:
	static continuous_set_ptr implement(const polyhedron<T>* p1,
			const polyhedron<T>* p2) {
		typename polyhedron<T>::ptr res = typename polyhedron<T>::ptr(p1->clone());
		res->add_constraints(*p2->get_constraints());
		return res;
	}
	;
};

/** Implementation for polyhedra of different types */
template<typename T1,typename T2> class intersection_operator_after_polycast<polyhedron<T1> , polyhedron<T2> > {
public:
	static continuous_set_ptr implement(const polyhedron<T1>* p1,
			const polyhedron<T2>* p2) {
		T1 test1;
		T2 test2;
		std::string name1 = typeid(test1).name();
		std::string name2 = typeid(test2).name();
		throw std::runtime_error("intersection_operator not defined for polyhedra of scalar types "
				+ name1 + " and " + name2 + " .");
		return continuous_set_ptr();
	}
	;
};

typedef constr_polyhedron<Rational> con_poly;
typedef constr_polyhedron<double> con_poly_double;

//template<> class intersection_operator<con_poly, con_poly> {
//public:
//	static continuous_set_ptr implement(const con_poly* p1, const con_poly* p2);
//};

template<> class intersection_operator<ppl_polyhedron::continuous_set_PPL_NNC,
		ppl_polyhedron::continuous_set_PPL_NNC> {
public:
	static continuous_set_ptr implement(
			const ppl_polyhedron::continuous_set_PPL_NNC* p1,
			const ppl_polyhedron::continuous_set_PPL_NNC* p2);
};

/** Implementation for polyhedra with sfm_cont_set */
template<typename T> class intersection_operator_after_polycast<polyhedron<T> , support_function::sfm_cont_set<T> > {
public:
	static continuous_set_ptr implement(const polyhedron<T>* p1,
			const support_function::sfm_cont_set<T>* p2) {
		//continuous_set_ptr res = p2->intersection_with_poly_improved(*p1); // sfm_cont_set Intersection constraint_poly
		typename polyhedron_collection<T>::ptr res =
				typename polyhedron_collection<T>::ptr(
						new polyhedron_collection<T>(
								p2->get_outer_polytope_collection()));
		res->add_constraints(*p1->get_constraints());
		res->remove_empty();
//		std::cout << res << std::endl;
		return res;
	}
	;
};

/** Implementation for polyhedra with spacetime_flowpipe */
template<typename T> class intersection_operator_after_polycast<polyhedron<T> , spacetime_flowpipe<T> > {
public:
	static continuous_set_ptr implement(const polyhedron<T>* p1,
			const spacetime_flowpipe<T>* p2) {

		typename spacetime_flowpipe<T>::cut_point_method m; // use default: all pieces

		//continuous_set_ptr res = p2->intersection_with_poly_improved(*p1); // sfm_cont_set Intersection constraint_poly
		typename polyhedron_collection<T>::ptr res =
				typename polyhedron_collection<T>::ptr(
						new polyhedron_collection<T>(
								p2->compute_outer_polyhedra(m)));
		res->add_constraints(*p1->get_constraints());
		res->remove_empty();
//		std::cout << res << std::endl;
		return res;
	}
	;
};

/** Implementation for polyhedra with a simulation trajectory */
template<typename T> class intersection_operator_after_polycast<polyhedron<T> , continuous_set_simulation<T> > {
public:
	static continuous_set_ptr implement(const polyhedron<T>* p1,
			const continuous_set_simulation<T>* p2) {
		typename math::lin_constraint_system<T>::const_ptr cons = p1->get_constraints();
		if (cons) {
			return p2->intersection_with(*cons);
		} else {
			return continuous_set_ptr(p2->clone());
		}
	}
	;
};

//template<> class intersection_operator<con_poly,
//		ppl_polyhedron::continuous_set_PPL_NNC> {
//public:
//	static continuous_set_ptr implement(const con_poly* p1,
//			const ppl_polyhedron::continuous_set_PPL_NNC* p2);
//};
//
//template<> class intersection_operator<ppl_polyhedron::continuous_set_PPL_NNC,
//		con_poly> {
//public:
//	static continuous_set_ptr
//	implement(const ppl_polyhedron::continuous_set_PPL_NNC* p1,
//			const con_poly* p2);
//};

template<> class intersection_operator<support_function::sfm_cont_set<Rational> , con_poly> {
public:
	static continuous_set_ptr implement(const support_function::sfm_cont_set<Rational>* p1,
			const con_poly* p2);
};
template<> class intersection_operator<support_function::sfm_cont_set<double> , con_poly_double> {
public:
	static continuous_set_ptr implement(const support_function::sfm_cont_set<double>* p1,
			const con_poly_double* p2);
};

template<> class intersection_operator<predicate_continuous_set,
		predicate_continuous_set> {
public:
	static continuous_set_ptr implement(const predicate_continuous_set* p1,
			const predicate_continuous_set* p2);
};


template<> class intersection_operator< continuous_set_simulation<Rational>,
continuous_set_simulation<Rational> > {
public:
	static continuous_set_ptr implement(const continuous_set_simulation<Rational>* p1,
			const continuous_set_simulation<Rational>* p2);
};

template<> class intersection_operator<continuous_set_simulation<double>,
continuous_set_simulation<double> > {
public:
	static continuous_set_ptr implement(const continuous_set_simulation<double>* p1,
			const continuous_set_simulation<double>* p2);
};

}

#endif /*INTERSECTION_H_*/
