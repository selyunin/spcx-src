#ifndef LP_SOLVER_USER_H_
#define LP_SOLVER_USER_H_

#include "math/lp_solving/lp_solver_chooser.h"

namespace math {

/** A base class that provides access to an lp solver chosen via the
 * lp_solver_chooser. */

template<typename scalar_type> class lp_solver_user {
public:
	// --------------------------------------------
	// Constructors
	// --------------------------------------------
	lp_solver_user() :
			my_lp_solver(0) {
		initialize();
	}
	;

	lp_solver_user(math::lp_solver<scalar_type>* s) :
			my_lp_solver(0) {
		set_lp_solver(s);
	}
	;

	/** Override default copy constructor, to create a new solver */
	lp_solver_user(const lp_solver_user<scalar_type>& s) :
			my_lp_solver(0) {
		initialize();
	}
	;

	/** Override default assignment constructor, to create a new solver */
	lp_solver_user<scalar_type>& operator=(
			const lp_solver_user<scalar_type>& s) {
		if (this != &s) {
			initialize();
		}
		return *this;
	}
	;

	virtual ~lp_solver_user() {
		destruct_solver();
	}
	;

	/** Choose the lp solver s. */
	void set_lp_solver(math::lp_solver<scalar_type>* s) {
		destruct_solver();
		my_lp_solver = s;
		my_changed = true;
	}
	;

	math::lp_solver<scalar_type>& get_lp_solver() {
		if (!my_lp_solver) {
			create_new_solver();
		}
		return *my_lp_solver;
	}
	;

	void swap(lp_solver_user& s) {
		/* the following is wrong, since set_lp_solver deletes the previous solver */
//		math::lp_solver<scalar_type>* tmp=get_lp_solver();
//		set_lp_solver(s.get_lp_solver());
//		s.set_lp_solver(tmp);
		std::swap(my_lp_solver, s.my_lp_solver);
		std::swap(my_changed, s.my_changed);
	}
	;

	/** Returns true after construction and if solver changed since last call to has_changed
	 *
	 * @attention Only reports true once, unless solver changes afterwards. */
	bool lp_solver_has_changed() {
		bool ret_val = my_changed;
		my_changed = false;
		return ret_val;
	}

protected:
	/** Initialize
	 *
	 * Could be either to get a new instance or delay (get zero). */
	void initialize() {
		if (true) {
			my_lp_solver = 0;
			my_changed = true;
		} else {
			create_new_solver();
		}
	}

	/** Destruct the solver if there is one */
	void destruct_solver() {
		if (my_lp_solver) {
//			LOGGER_OS(DEBUG7, __FUNCTION__) << "user destroying solver "
//					<< my_lp_solver << " of user " << this;
			delete my_lp_solver;
			/** @note delete does not necessarily set pointer to null, so we need to do it by hand */
			my_lp_solver = 0;
		}
	}

	/** Obtain a new solver */
	void create_new_solver() {
		set_lp_solver(math::lp_solver_chooser<scalar_type>::get_instance());
//		LOGGER_OS(DEBUG7, __FUNCTION__) << "user created solver "
//				<< my_lp_solver << " of user " << this;
	}

private:
	math::lp_solver<scalar_type>* my_lp_solver;
	bool my_changed; // true after construction and if solver changed since last call to  has_changed
};

}

#endif /*LP_SOLVER_USER_H_*/
