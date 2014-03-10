#ifndef DOUBLE_DISPATCH_H_
#define DOUBLE_DISPATCH_H_

/** Recursive double dispatching, visitor-based implementation of multimethods.
 * It converts two base pointers to the corresponding derived pointers using
 * a visitor (dispatcher).
 * See Shopyrin: Multimethods Multimethods in C++ Using Recursive Deferred Dispatching
 * */

#include "utility/dispatching/impl/cross_dispatcher.h"
#include "utility/dispatching/impl/template_class_implementor_wrapper.h"

namespace dispatching {

/** Double dispatch to implementor class.
 *
 * The class implementor must provide a function
 * return_type implement(d1,d2) for any derived type d1 and d2.
 * A useful way to guarantee this is to let implement be a template
 * member function.
 * However, this cannot be partially specialized.
 * For partial specializations, use double_dispatch_tc.
 */
template<typename return_type, class implementor, class base_class,
		typename types_to_visit> return_type double_dispatch(
		const base_class* base1, const base_class* base2) {
	cross_dispatcher<return_type, implementor, base_class, types_to_visit,
			types_to_visit> d(base2);
	base1->accept(d);
	return d.get();

}
;

/** Double dispatch to templated implementor class.
 *
 * The class implementor<T1,T2> must be templated, and will be instantiated
 * for all combinations of derived classes.
 * The class must provide
 * return_type implement(const T1* t1,const T2* t2) for derived types T1 and T2.
 * Note that the class can be partially specialized.
 */
template<typename return_type,
		template<typename , typename > class template_class_implementor,
		class base_class, typename types_to_visit> return_type double_dispatch_tc(
		const base_class* base1, const base_class* base2) {
	typedef template_class_double_implementor_wrapper<return_type,
			template_class_implementor> implementor;
	return double_dispatch<return_type, implementor, base_class, types_to_visit> (
			base1, base2);
}
;

/** Double dispatch to templated implementor class with casting.
 *
 * The class implementor<T1,T2> must be templated, and will be instantiated
 * for all combinations of derived classes.
 * The class must provide
 * return_type implement(const T1* t1,const T2* t2) for derived types T1 and T2.
 * Note that the class can be partially specialized.
 * The derived classes will be cast according to the provided casters.
 */
template<typename return_type,
		template<typename , typename > class template_class_implementor,
		class base_class, typename types_to_visit,
		template<typename , typename > class caster> return_type double_dispatch_tc_cast(
		const base_class* base1, const base_class* base2) {
	typedef template_class_double_implementor_wrapper<return_type,
			template_class_implementor, caster> implementor;
	return double_dispatch<return_type, implementor, base_class, types_to_visit> (
			base1, base2);
}
;

}

#endif /*DOUBLE_DISPATCH_H_*/
