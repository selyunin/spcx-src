#ifndef vdom_VECTOR_OUTPUT_H_
#define vdom_VECTOR_OUTPUT_H_

namespace math {

using ::operator<<;

template<typename scalar_type>
class vdom_vector<scalar_type>::output_format {
public:
	output_format() { // default values
		preamble="[";
		epilogue="]";
		element_separator=",";
		value_and_name_separator="=";
		write_variable_names=true;
		write_values_as_double=false;
	}
	;
	static output_format matlab() {
		output_format f;
		f.preamble="[";
		f.epilogue="]";
		f.element_separator=",";
		f.write_variable_names=false;
		f.write_values_as_double=true;
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
		return f;
	}
	;
	std::string preamble; // written before matrix
	std::string epilogue; // written after matrix
	std::string element_separator; // written between elements
	std::string value_and_name_separator; // written between name and value

	bool write_variable_names;
	bool write_values_as_double;
};

template<typename scalar_type>
void vdom_vector<scalar_type>::print(std::ostream& os) const {
	const vdom_vector<scalar_type>& v=*this;
	const typename vdom_vector<scalar_type>::output_format& of=v.get_output_format(); 
	os << of.preamble;
	for (unsigned int i=0; i<v.size(); i++) {
		//if (v[i]!=scalar_type(0)) {
			if (i>0) {
				os << of.element_separator;
			}
			if (of.write_variable_names) {
				os << variable(v.get_index_to_variable_id_map()->get_id(i));
				os << of.value_and_name_separator;
			}
			if (of.write_values_as_double)
				os << convert_element<double,scalar_type>(v[i]);
			else
				os << v[i];
		//}
	}
	os << of.epilogue;		
}

}
#endif /*vdom_VECTOR_OUTPUT_H_*/
