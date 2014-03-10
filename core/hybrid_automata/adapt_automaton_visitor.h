/*
 * adapt_automaton_visitor.h
 *
 *  Created on: Sep 11, 2009
 *      Author: frehse
 */

#ifndef ADAPT_AUTOMATON_VISITOR_H_
#define ADAPT_AUTOMATON_VISITOR_H_

#include "core/discrete/discrete_set.h"
#include "core/continuous/continuous_set.h"
#include "core/continuous/continuous_dynamics/continuous_dynamics_base.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"
#include "core/hybrid_automata/hybrid_automaton_visitor.h"
#include "global/global_types.h"
//#include "assignment_traitment.h"
//#include "../parser/automaton_creation.h"
//#include "../abstract_framework/continuous/continuous_dynamics/continuous_dynamics.h"
//#include "../abstract_framework/continuous/continuous_set_transforms/continuous_set_transform_composition.h"
//#include "../abstract_framework/continuous/predicate_continuous_set.h"
//#include "../set_implementations/convert_predicate.h"

namespace hybrid_automata {

class adapt_automaton_visitor;
class adapt_continuous_set_visitor;
class adapt_dynamics_visitor;
class adapt_transform_visitor;

typedef boost::shared_ptr<adapt_automaton_visitor> adapt_automaton_visitor_ptr;
typedef boost::shared_ptr<adapt_continuous_set_visitor>
		adapt_continuous_set_visitor_ptr;
typedef boost::shared_ptr<adapt_dynamics_visitor> adapt_dynamics_visitor_ptr;
typedef boost::shared_ptr<adapt_transform_visitor> adapt_transform_visitor_ptr;

class adapt_discrete_set_visitor: public discrete::discrete_set::const_visitor {
public:
	adapt_discrete_set_visitor();
	virtual ~adapt_discrete_set_visitor();
	void reset();
	bool get_success();
	discrete::discrete_set::ptr get_discrete_set();
protected:
	bool success;
	discrete::discrete_set::ptr my_discrete_set;
};

class adapt_continuous_set_visitor: public continuous::continuous_set::const_visitor {
public:
	adapt_continuous_set_visitor();
	virtual ~adapt_continuous_set_visitor();
	void reset();
	bool get_success();
	continuous::continuous_set::ptr get_continuous_set();
protected:
	bool success;
	continuous::continuous_set::ptr my_continuous_set;
};

/** Default implementation for constant_bound_dynamics and relation_dynamics
 * is to convert the underlying continuous_set using the adapt_continuous_set_visitor
 * passed in the init method.
 * Optionally, the invariant can be set so that dynamics can take
 * it into account.
 */
class adapt_dynamics_visitor: public continuous::continuous_dynamics::const_visitor {
public:
	adapt_dynamics_visitor();
	void init(adapt_continuous_set_visitor_ptr con_adn);
	void set_invariant(continuous::continuous_set::ptr inv);
	virtual ~adapt_dynamics_visitor();
	void reset();

	virtual void dispatch(const continuous::constant_bound_dynamics* c);
	virtual void dispatch(const continuous::relation_dynamics* c);

	bool get_success();
	continuous::continuous_dynamics::ptr get_dynamics() const;
	continuous::continuous_set::ptr get_invariant() const;
protected:
	bool success;
	continuous::continuous_dynamics::ptr my_dynamics;
	continuous::continuous_set::ptr my_invariant;
	adapt_continuous_set_visitor_ptr my_con_ad; // in case continuous_sets need adapting (relation_dynamics)
};

/** Default implementation for relation_transform,
 * constant_bound_time_elapse_transform, intersection_transform is to convert the
 * underlying continuous_set using the adapt_continuous_set_visitor
 * passed in the init method.
 */
class adapt_transform_visitor: public continuous::continuous_set_transform::const_visitor {
public:
	adapt_transform_visitor();
	void init(adapt_continuous_set_visitor_ptr con_ad);
	virtual ~adapt_transform_visitor();

	virtual adapt_transform_visitor* clone();
	void reset();

	virtual void dispatch(const continuous::relation_transform* c);
	virtual void dispatch(
			const continuous::constant_bound_time_elapse_transform* c);
	virtual void dispatch(const continuous::intersection_transform* c);
	virtual void dispatch(const continuous::reset_affine_transform<
			global_types::rational_type>* c);
	virtual void dispatch(const continuous::reset_affine_transform<
			global_types::float_type>* c);
	virtual void dispatch(const continuous::reset_function_transform* c);
	virtual void dispatch(const continuous::sequence_transform* c);

	bool get_success();
	continuous::continuous_set_transform::ptr get_transform();
protected:
	bool success;
	continuous::continuous_set_transform::ptr my_transform;
	adapt_continuous_set_visitor_ptr my_con_ad; // in case continuous_sets need adapting (relation_transform)
};

class adapt_automaton_visitor: public hybrid_automaton_visitor {
public:
	typedef global_types::coefficient_type coefficient_type;

	adapt_automaton_visitor();
	virtual ~adapt_automaton_visitor();
	virtual adapt_automaton_visitor* clone() const;
	/** The init call should be overridden by derived classes, which should call define() with the
	 * appropriate visitors.
	 */
	virtual void init();
	void define(adapt_continuous_set_visitor_ptr con_ad,
			adapt_dynamics_visitor_ptr dyn_ad,
			adapt_transform_visitor_ptr trans_ad, coefficient_type bool_t,
			coefficient_type number_t);
	virtual void visit(hybrid_automaton& h, location& l, location_id l_id);
	virtual void visit(hybrid_automaton& h, transition& t, transition_id t_id);
	virtual void epilogue(hybrid_automaton& h);
	virtual void reset();
	virtual bool get_success() const;
	adapt_continuous_set_visitor_ptr get_adapt_continuous_set_visitor();
	adapt_dynamics_visitor_ptr get_adapt_dynamics_visitor();
	adapt_transform_visitor_ptr get_adapt_transform_visitor();
	coefficient_type get_bool_type() const;
	coefficient_type get_number_type() const;
protected:
	void set_bool_type(coefficient_type t);
	void set_number_type(coefficient_type t);
	unsigned int my_loc_count;
	unsigned int my_trans_count;
private:
	void successful_visit();
	void failed_visit();
	bool my_success;
	adapt_continuous_set_visitor_ptr my_con_ad;
	adapt_dynamics_visitor_ptr my_dyn_ad;
	adapt_transform_visitor_ptr my_trans_ad;
	coefficient_type my_bool_type;
	coefficient_type my_number_type;
};

}

#endif /* ADAPT_AUTOMATON_VISITOR_H_ */
