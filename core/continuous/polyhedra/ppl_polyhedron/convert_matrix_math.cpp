#include "core/continuous/polyhedra/ppl_polyhedron/convert_matrix_math.h"

namespace ppl_polyhedron {

using namespace std;
using namespace Parma_Polyhedra_Library;

/** Converts a Generator object to a rational vector
 */
math::rational_vector convert_to_rational_vector(
		const Generator& g) {
	math::rational_vector p=
			math::rational_vector(g.space_dimension());
//	mpq_class q(0);
	for (unsigned int i=0; i<g.space_dimension(); ++i) {
		p[i]=Rational(g.coefficient(Variable(i)), g.divisor());
	}
	return p;
}

math::lin_expression<Rational> convert_to_lin_expression(
		const Generator& g, const index_to_variable_id_map_ptr& iimap) {
	return math::lin_expression<Rational>(convert_to_rational_vector(g),Rational(0),iimap);
}

math::vdom_vector<Rational> convert_to_vdom_vector(
		const Generator& g, const index_to_variable_id_map_ptr& iimap) {
	return math::vdom_vector<Rational>(convert_to_rational_vector(g),iimap);
}

math::lin_constraint<Rational> convert_to_lin_constraint(
		const Parma_Polyhedra_Library::Constraint& c, const index_to_variable_id_map_ptr& iimap)
{
	math::rational_vector p=
			math::rational_vector(c.space_dimension());
//	mpq_class q(0);
	for (unsigned int i=0; i<c.space_dimension(); ++i) {
		p[i]=Rational(c.coefficient(Variable(i))); 
	}	
	/* convert from a1 x1 + ... an xn + b >= 0 (PPL) to a1 x1 + ... an xn + b >= 0 (lin_constraint) */ 
	Rational b=Rational(c.inhomogeneous_term());// Note : b has to change sign
	
	math::lin_expression<Rational> l=math::lin_expression<Rational>(p,b,iimap);
	if (c.is_equality())
		return math::lin_constraint<Rational>(l,EQ);
	else if(c.is_nonstrict_inequality())
		return math::lin_constraint<Rational>(l,GE);
	else
		return math::lin_constraint<Rational>(l,GT);
}
/**
 * Converts a Generator object to a double vector
 */

math::double_vector convert_to_double_vector(const Generator& g) {
	math::double_vector p=
			math::double_vector(g.space_dimension());
	Rational temp;

	for (unsigned int i=0; i<g.space_dimension(); ++i) {
		temp = Rational(g.coefficient(Variable(i)), g.divisor());
		p[i] = temp.get_double();
	}
	return p;
}

/** Converts the rational_vector v into a linear expression le and a denominator d such that le[i]/d = v[i]. */
void convert_to_Linear_Expression(const math::rational_vector& v,
		Linear_Expression& le, Integer& d) {
	le = Linear_Expression();
	d = Integer(1);
	if (v.size()>0) {
		// Get common denominator
		unsigned int maxdim=v.size()-1;
		for (unsigned int i = maxdim; i >= 0 && i<=maxdim; i--)
			d *= v[i].get_den();

		for (unsigned int i = maxdim; i >= 0 && i<=maxdim; i--) {
			Integer d_temp=Integer(1);
			for (unsigned int j = maxdim; j >= 0 && j<=maxdim; j--) {
				if (j!=i) {
					d_temp *= v[j].get_den();
				}
			}
			le += v[i].get_num() * d_temp * Variable(i);
		}
	}
}

/** Converts the rational_vector v into a linear expression le and a denominator d such that le[i]/d = v[i]. */
void convert_to_Linear_Expression(const math::double_vector& v,
		Linear_Expression& le, Integer& d) {
	math::rational_vector vr=v.convert_to<Rational>();
	convert_to_Linear_Expression(vr,le,d);
}

bool convert_to_Linear_Expression(
		const math::lin_expression<Rational>& l, Linear_Expression& le,
		Integer& denominator, const index_to_variable_id_map_ptr& iimap) {
	if (iimap==l.get_index_to_variable_id_map()) {
		convert_to_Linear_Expression(l.get_vector(), le, denominator);
		// Add the inhomogenous coefficient
		Integer b=l.get_inh_coeff().get_num()*denominator;
		le*=l.get_inh_coeff().get_den();
		denominator*=l.get_inh_coeff().get_den();
		le+=b;
	} else {
		bool has_id;
		// Start with the inhomogenous term
		le = Linear_Expression(l.get_inh_coeff().get_num());
		denominator = l.get_inh_coeff().get_den();
		if (l.size()>0) {
			// Get common denominator
			unsigned int maxdim=l.size()-1;
			for (unsigned int i = maxdim; i >= 0 && i<=maxdim; i--)
				denominator *= l[i].get_den();

			for (unsigned int i = maxdim; i >= 0 && i<=maxdim; i--) {
				Integer d_temp=Integer(l.get_inh_coeff().get_den());
				for (unsigned int j = maxdim; j >= 0 && j<=maxdim; j--) {
					if (j!=i) {
						d_temp *= l[j].get_den();
					}
				}
				index_type index_in_l=l.get_index_to_variable_id_map()->get_id(i);
				index_type j=iimap->check_for_index(index_in_l, has_id);
				if (!has_id && l[i]!=Rational(0)) {
					return true;
				}
				le += l[i].get_num() * d_temp * Variable(j);
			}
		}
	}
	return false; // everything ok, no other variables
}

}
