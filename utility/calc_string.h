#ifndef CALC_STRING_H_
#define CALC_STRING_H_

#include <stdexcept>
#include <string>
#include <typeinfo>
#include "math/type_conversion.h"

class calc_string {
public:
	calc_string();

	explicit calc_string(double b);

	explicit calc_string(int b);

	explicit calc_string(bool b);

	explicit calc_string(const std::string& s);

	calc_string(const calc_string& cs);

	virtual ~calc_string();

	const std::string& get_my_string() const;

	static bool get_parentheses();

	static void set_parentheses(bool b);

private:
	std::string my_string;
	static bool parentheses;
};

inline std::ostream& operator<<(std::ostream& os, const calc_string& s) {
	os << s.get_my_string();
	return os;
}

/** Don't convert from calc_string to anything,
 * but let the code compile. */
template<typename result_type>
class converter<result_type, calc_string> {
public:
	static result_type convert(const calc_string& x) {
		std::string name1=typeid(x).name();
		result_type rtmp;
		std::string name2=typeid(rtmp).name();
		throw std::runtime_error("convert_element from "
				+ name1 + " to "
				+ name2
				+ " not possible.");
		return result_type();
	}
	;
};

/** The above converter is ambiguous with the predefined converter if
 * scalar_type is Rational. The following full specialization resolves the
 * ambiguity. */
template<>
class converter<calc_string, calc_string> {
public:
	static calc_string convert(const calc_string& x) {
		return x;
	}
	;
};

#endif /*CALC_STRING_H_*/
