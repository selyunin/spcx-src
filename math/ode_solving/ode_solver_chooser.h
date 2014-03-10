/*
 * ode_solver_chooser.h
 *
 *  Created on: Oct 7, 2010
 *      Author: frehse
 */

#ifndef ODE_SOLVER_CHOOSER_H_
#define ODE_SOLVER_CHOOSER_H_

#include "math/ode_solving/ode_solver_cvode.h"

namespace math {
namespace ode {

template<typename scalar_type> class ode_solver;
template<typename scalar_type> class ode_solver_cvode;

template<typename scalar_type> class ode_solver_chooser {
public:
	static ode_solver<scalar_type>* get_instance() {
		ode_solver<scalar_type>* inst = new ode_solver_cvode<scalar_type> ();
		return inst;
	}
private:
	ode_solver_chooser(const ode_solver_chooser<scalar_type>&);
	ode_solver_chooser<scalar_type>& operator =(const ode_solver_chooser<
			scalar_type>& other);
};

}
}

#endif /* ODE_SOLVER_CHOOSER_H_ */
