#ifndef CONTAINMENT_H_
#define CONTAINMENT_H_

#include <stdexcept>
#include <typeinfo>
#include <boost/mpl/if.hpp>

#include "core/continuous/polyhedra/polyhedron_upcaster.h"
#include "core/continuous/support_function_provider_upcaster.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"

//for simulation
#include "math/ode_solving/traj_simu/continuous_set_simulation.h"

// For temporary fix until we can pass the containment operator in the scenario
#include "core/continuous/support_function/sf_base/sf_unary_template.h"
#include "core/continuous/polyhedra/polyhedron_collection.h"

/** Forward declarations */
namespace continuous {
class support_function_provider;
template<typename T> class polyhedron;
template<typename scalar_type> math::tribool compute_poly_containment(const polyhedron<
		scalar_type>& t1, const support_function_provider& t2);
}

namespace continuous {

/** Declaration of the containment operator. */
template<typename T1, typename T2> class containment_operator {
public:
	static math::tribool implement(const T1* t1, const T2* t2) {
		// Temporary fix until we can pass the containment operator in the scenario
		if (const support_function::sf_unary_template<double,
				polyhedron_collection<double> >* d2 =
				dynamic_cast<const support_function::sf_unary_template<double,
						polyhedron_collection<double> >*>(t2)) {
			if (const support_function::sf_unary_template<double,
					polyhedron_collection<double> >* d1 =
					dynamic_cast<const support_function::sf_unary_template<
							double, polyhedron_collection<double> >*>(t1)) {
				if (d1->get_map() || d2->get_map()) {
					throw basic_warning("containment_operator",
							"can't yet handle mapped sf_unary_templates.",
							basic_warning::MISSING_IMPLEMENTATION);
				}
				return d1->get_implementor()->element_wise_contains(
						*(d2->get_implementor()));
			} else {
				// call elementwise containment for each element of d2
				const polyhedron_collection<double>* c2 =
						static_cast<const polyhedron_collection<double>*>(d2->get_implementor().get());
				math::tribool contains_res(true);
				for (typename polyhedron_collection<double>::const_iterator it =
						c2->begin();
						it != c2->end() && math::maybe(contains_res); ++it) {
					contains_res = contains_res && t1->contains(*it);
				}
				return contains_res;
			}
		} else if (const support_function::sf_unary_template<double,
				polyhedron_collection<double> >* d1 =
				dynamic_cast<const support_function::sf_unary_template<double,
						polyhedron_collection<double> >*>(t1)) {

			// call elementwise containment for each element of d2
			const polyhedron_collection<double>* c1 =
					static_cast<const polyhedron_collection<double>*>(d1->get_implementor().get());
			math::tribool contains_res(false);
			for (typename polyhedron_collection<double>::const_iterator it =
					c1->begin(); it != c1->end() && math::maybe(!contains_res);
					++it) {
				contains_res = contains_res
						&& containment_test((*it).get(), t2);
			}
			return contains_res;
		}
		// check if it's two mapped sf_sets with the same map
		if (const support_function::sf_unary<double>* d2 =
				dynamic_cast<const support_function::sf_unary<double>*>(t2)) {
			if (const support_function::sf_unary<double>* d1 =
					dynamic_cast<const support_function::sf_unary<double>*>(t1)) {
				// only check if the maps are the same and the bloating is contained
				bool check = true;
				if (!d1->get_map()) {
					if (!d2->get_map()) {
						// check containment between root of d1 and root of d2
						return d1->get_unmapped_set()->contains(d2->get_unmapped_set());
					} else {
						// check containment between root of d1 and d2
						return d1->get_unmapped_set()->contains(d2->get_const_ptr());
					}
				} else {
					if (!d2->get_map()) {
						// check containment between d1 and root of d2
						// @note we could try to invert the map of d1 here and check a mapped d2 against the root of d1
						return d1->contains(d2->get_unmapped_set());
					} else {
						// check if both maps are the same
						if (*d1->get_map() == *d2->get_map()) {
							return d1->get_unmapped_set()->contains(d2->get_unmapped_set());
						}
					}
				}
			}
		}


		std::string name1 = typeid(t1).name();
		std::string name2 = typeid(t2).name();
		std::string real_name1 = typeid(*t1).name();
		std::string real_name2 = typeid(*t2).name();
		basic_warning("containment_operator","containment_operator not defined for "
				+ name1 + " and " + name2 + ". Actual types are " + real_name1 + " and " + real_name2 + ".",basic_warning::MISSING_IMPLEMENTATION);
		IFLOGGER(DEBUG7) {
			LOGGER_OS(DEBUG7, __FUNCTION__) << "sets are " << *t1 << " contains(?) "<< *t2 << std::endl;
			throw std::runtime_error("missing implementation of containment operator");
		}
		return math::indeterminate();
	}
	;
};

/** If T1 is a polyhedron, try to cast T2 to a support_function_provider. */
template<typename T1, typename T2>
class containment_test_upcaster {
public:
	typedef typename polyhedron_upcaster<T1>::result result1;
	typedef typename boost::mpl::if_c<polyhedron_upcaster<T1>::is_poly,
			typename support_function_provider_upcaster<T2>::result,
			typename polyhedron_upcaster<T2>::result>::type
			result2;
};

/** Implementation for polyhedra with support_function_provider */
template<typename T> class containment_operator<polyhedron<T> ,
		support_function_provider> {
public:
	static math::tribool implement(const polyhedron<T>* p1,
			const support_function_provider* p2) {
		return compute_poly_containment<T> (*p1, *p2);
	}
	;
};

///** Implementation for support_function_provider with polyhedra*/
//template<typename T> class containment_operator<support_function_provider,
//					polyhedron<T> > {
//public:
//	static bool implement(const polyhedron<T>* p1,
//			const support_function_provider* p2) {
//		return compute_poly_containment<T> (*p1, *p2);
//	}
//	;
//};
/** Implementation for sfm_cont_set */
template<typename T> class containment_operator<support_function::sfm_cont_set<T> , support_function::sfm_cont_set<
		T> > {
public:
	static math::tribool implement(const support_function::sfm_cont_set<T>* p1,
			const support_function::sfm_cont_set<T>* p2) {
		return p1->contains_initial_set(*p2);
	};
};

/** Implementation for sfm_cont_set */
template<typename T> class containment_operator<support_function::sfm_cont_set<T> , polyhedron<T> > {
public:
	static math::tribool implement(const support_function::sfm_cont_set<T>* p1,
			const polyhedron<T>* p2) {
		return p1->contains_initially(p2->get_const_ptr());
	};
};

/** Implementation for continuous_set_simulation */
template<typename T> class containment_operator<continuous_set_simulation<T> , continuous_set_simulation<
		T> > {
public:
	static math::tribool implement(const continuous_set_simulation<T>* p1,
			const continuous_set_simulation<T>* p2) {
		return p1->contains_roots_of(*p2);
	};
};

template<typename T> class containment_operator<continuous_set_simulation<T> , polyhedron<
		T> > {
public:
	static math::tribool implement(const continuous_set_simulation<T>* p1,
			const polyhedron<T>* p2) {
		//return p1->contains_roots_of(*p2);
		return math::indeterminate();
	};
};


}

#endif /*CONTAINMENT_H_*/
