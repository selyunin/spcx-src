#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"

std::ostream& operator<<(std::ostream& os, const continuous::constr_polyhedron<
		calc_string>::ptr& p) {
	p->get_poly().print(os);
	return os;
}

