#include "core/continuous/continuous_set_operator_implementations/intersection.h"

#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/continuous/support_function/sfm/sfm_cont_set.h"
#include "core/continuous/predicate_continuous_set.h"

namespace continuous {

//continuous_set_ptr intersection_operator<con_poly, con_poly>::implement(
//		const con_poly* p1, const con_poly* p2) {
//	con_poly::ptr res(new con_poly(*p1));
//	res->intersection_assign(*p2);
//	return res;
//}

continuous_set_ptr intersection_operator<
		ppl_polyhedron::continuous_set_PPL_NNC,
		ppl_polyhedron::continuous_set_PPL_NNC>::implement(
		const ppl_polyhedron::continuous_set_PPL_NNC* p1,
		const ppl_polyhedron::continuous_set_PPL_NNC* p2) {
	ppl_polyhedron::continuous_set_PPL_NNC::ptr pret(
			new ppl_polyhedron::continuous_set_PPL_NNC(*p1));
	pret->intersection_assign(*p2);
	return pret;
}

//continuous_set_ptr intersection_operator<
//		ppl_polyhedron::continuous_set_PPL_NNC, con_poly>::implement(
//		const ppl_polyhedron::continuous_set_PPL_NNC* p1, const con_poly* p2) {
//	ppl_polyhedron::continuous_set_PPL_NNC::ptr res(
//			new ppl_polyhedron::continuous_set_PPL_NNC(*p1));
//	res->add_constraints(p2->get_constraints());
//	return res;
//}
//
//continuous_set_ptr intersection_operator<con_poly,
//		ppl_polyhedron::continuous_set_PPL_NNC>::implement(const con_poly* p1,
//		const ppl_polyhedron::continuous_set_PPL_NNC* p2) {
//	return intersection_operator<ppl_polyhedron::continuous_set_PPL_NNC,
//			con_poly>::implement(p2, p1);
//}

continuous_set_ptr intersection_operator<support_function::sfm_cont_set<Rational> , con_poly>::implement(
		const support_function::sfm_cont_set<Rational>* p1, const con_poly* p2) {
	return p1->intersection_with_poly_improved(*p2); // sfm_cont_set Intersection constraint_poly
//	return p1->intersection_with_poly(*p2); // sfm_cont_set Intersection constraint_poly
}

continuous_set_ptr intersection_operator<support_function::sfm_cont_set<double> , con_poly_double>::implement(
		const support_function::sfm_cont_set<double>* p1, const con_poly_double* p2) {
	continuous_set_ptr res = p1->intersection_with_poly_improved(*p2); // sfm_cont_set Intersection constraint_poly
	return res;
//	return p1->intersection_with_poly(*p2); // sfm_cont_set Intersection constraint_poly
}

continuous_set_ptr intersection_operator<predicate_continuous_set,
		predicate_continuous_set>::implement(
		const predicate_continuous_set* p1, const predicate_continuous_set* p2) {
	predicate_continuous_set::ptr pnew(p1->clone());
	pnew->intersection_assign(*p2);
	return pnew;
}

continuous_set_ptr intersection_operator<continuous_set_simulation<Rational>,
	continuous_set_simulation<Rational> >::implement(
		const continuous_set_simulation<Rational>* p1, const continuous_set_simulation<Rational>* p2) {
	return  p1->intersection_with(*p2);
}

continuous_set_ptr intersection_operator<continuous_set_simulation<double>,
	continuous_set_simulation<double> >::implement(
		const continuous_set_simulation<double>* p1, const continuous_set_simulation<double>* p2) {
	return  p1->intersection_with(*p2);
}

}

