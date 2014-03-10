/*
 * ode_solver.cpp
 *
 *  Created on: Oct 22, 2010
 *      Author: frehse
 */

#include "ode_solver.h"

#include "utility/basic_exception.h"

namespace math {
namespace ode {

double ode_defaults::my_abs_tol = 1.0e-6;
double ode_defaults::my_rel_tol = 1.0e-9;

double ode_defaults::get_abs_tol() {
	return my_abs_tol;
}
double ode_defaults::get_rel_tol() {
	return my_rel_tol;
}
void ode_defaults::set_abs_tol(double x) {
	if (x < 0.0)
		throw basic_exception("Cannot set absolute tolerance for ode solver to a negative value.");
	my_abs_tol=x;
}
void ode_defaults::set_rel_tol(double x) {
	if (x < 0.0)
		throw basic_exception("Cannot set relative tolerance for ode solver to a negative value.");
	my_rel_tol=x;
}


}
}
