#ifndef LIN_EXPRESSION_OUTPUT_H_
#define LIN_EXPRESSION_OUTPUT_H_

#include "global/global_types.h" // for __float128 operator<<
#include "core/predicates/dot_context_lookup.h"

#include "math/matrix.h"

namespace math {

using ::operator<<;

template<typename scalar_type> class lin_expression<scalar_type>::output_format {
public:
	output_format() { // default values
		preamble="";
		epilogue="";
		element_separator="+";
		value_and_name_separator="*";
		write_variable_names=true;
		write_values_as_double=false;
		omit_zero_coeffs=true;
	}
	;
	static output_format matlab() {
		output_format f;
		f.preamble="[";
		f.epilogue="]";
		f.element_separator=",";
		f.write_variable_names=false;
		f.write_values_as_double=true;
		f.omit_zero_coeffs=false;
		return f;
	}
	;
	static output_format space_separated() {
		output_format f;
		f.preamble="";
		f.epilogue="";
		f.element_separator=" ";
		f.write_variable_names=false;
		f.write_values_as_double=true;
		f.omit_zero_coeffs=false;
		return f;
	}
	;
	std::string preamble; // written before matrix
	std::string epilogue; // written after matrix
	std::string element_separator; // written between elements
	std::string value_and_name_separator; // written between name and value

	bool write_variable_names;
	bool write_values_as_double;
	bool omit_zero_coeffs;
};

template<typename scalar_type> void lin_expression<scalar_type>::print(
		std::ostream& os) const {
	bool no_output_so_far=true;
	const lin_expression<scalar_type>& v=*this;
	const typename lin_expression<scalar_type>::output_format& of=v.get_output_format();

	for (unsigned int i=0; i<v.size(); i++) {
		if (!of.omit_zero_coeffs || v[i]!=scalar_type(0)) {
			if (!no_output_so_far && (!of.write_variable_names || v[i]>scalar_type(0))) {
				os << of.element_separator;
			}
			if (of.write_variable_names && v[i]==-scalar_type(1)) {
				os << "-";
			} else if (!of.write_variable_names || v[i]!=scalar_type(1)) {
				if (of.write_values_as_double)
					os << convert_element<double,scalar_type>(v[i]);
				else
					os << v[i];
				if (of.write_variable_names) 
					os << of.value_and_name_separator;
			}
			if (of.write_variable_names) {
				os << variable(v.get_index_to_variable_id_map()->get_id(i));
			}
			no_output_so_far=false;
		}
	}
	if (no_output_so_far || !of.omit_zero_coeffs || v.get_inh_coeff()!=scalar_type(0)) {
		if (!of.write_variable_names || (!no_output_so_far && v.get_inh_coeff()>=scalar_type(0))) {
			os << of.element_separator;
		}
		if (of.write_values_as_double)
			os << convert_element<double,scalar_type>(v.get_inh_coeff());
		else
			os << v.get_inh_coeff();
		no_output_so_far=false;
	}
}
;

template<typename scalar_type> void lin_expression<scalar_type>::print_hom(
		std::ostream& os) const {
	bool no_output_so_far=true;
	const lin_expression<scalar_type>& v=*this;
	const typename lin_expression<scalar_type>::output_format& of=v.get_output_format();

	for (unsigned int i=0; i<v.size(); i++) {
		if (!of.omit_zero_coeffs || v[i]!=scalar_type(0)) {
			if (!no_output_so_far && (!of.write_variable_names || v[i]>scalar_type(0))) {
				os << of.element_separator;
			}
			if (of.write_variable_names && v[i]==-scalar_type(1)) {
				os << "-";
			} else if (!of.write_variable_names || v[i]!=scalar_type(1)) {
				if (of.write_values_as_double)
					os << convert_element<double,scalar_type>(v[i]);
				else
					os << v[i];
				if (of.write_variable_names)
					os << of.value_and_name_separator;
			}
			if (of.write_variable_names) {
				os << variable(v.get_index_to_variable_id_map()->get_id(i));
			}
			no_output_so_far=false;
		}
	}
	if (no_output_so_far) {
		if (of.write_values_as_double)
			os << convert_element<double,scalar_type>(scalar_type(0));
		else
			os << scalar_type(0);
		no_output_so_far=false;
	}
}
;


/** Print a matrix as linear expressions over a domain */
template<typename scalar_type> void print_lin_expressions(std::ostream& os, const matrix<scalar_type>& A, positional_vdomain dom) {
	for (size_t i = 0; i<A.size1(); ++i) {
		vdom_vector<scalar_type> vec(dom,A.vector_from_row(i));
		lin_expression<scalar_type> expr(vec);
		os << expr << std::endl;
	}
}

}

#endif /*LIN_EXPRESSION_OUTPUT_H_*/
