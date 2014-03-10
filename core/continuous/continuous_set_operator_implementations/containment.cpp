#include "core/continuous/continuous_set_operator_implementations/containment.h"

namespace continuous {

//bool containment_test(const continuous_set_const_ptr& p1, const continuous_set_const_ptr& p2){
//	return dispatching::double_dispatch_tc_cast<bool, containment_operator,
//	continuous_set, continuous_set_typelist,polyhedron_upcaster,support_function_provider_upcaster>(p1.get(), p2.get());
//}
//
//bool containment_operator<polyhedron<double>, support_function_provider>::implement(
//		const polyhedron<double>* p1, const support_function_provider* p2) {
//	return compute_poly_containment<double> (*p1, *p2);
//}
//
//bool containment_operator<polyhedron<Rational>, support_function_provider>::implement(
//		const polyhedron<Rational>* p1, const support_function_provider* p2) {
//	return compute_poly_containment<Rational> (*p1, *p2);
//}

}
