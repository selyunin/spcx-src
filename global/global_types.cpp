/*
 * global_types.cpp
 *
 *  Created on: Nov 10, 2009
 *      Author: frehse
 */

#include "global/global_types.h"

namespace global_types {

//template<typename T> coefficient_type type_identifier<T>::coeff =
//		UNDEFINED_TYPE;

template<> const coefficient_type
		type_identifier<type_selector<STD_BOOL>::type>::coeff = STD_BOOL;
template<> const coefficient_type
		type_identifier<type_selector<STD_INT>::type>::coeff = STD_INT;
template<> const coefficient_type
		type_identifier<type_selector<STD_DOUBLE>::type>::coeff = STD_DOUBLE;
template<> const coefficient_type
		type_identifier<type_selector<STD_LONG_DOUBLE>::type>::coeff = STD_LONG_DOUBLE;
template<> const coefficient_type
		type_identifier<type_selector<STD_FLOAT128>::type>::coeff = STD_FLOAT128;
template<> const coefficient_type
		type_identifier<type_selector<GMP_RATIONAL>::type>::coeff =
				GMP_RATIONAL;
template<> const coefficient_type
		type_identifier<type_selector<CALC_STR>::type>::coeff = CALC_STR;

//template<typename T> std::string type_identifier<T>::name = "UNDEFINED_TYPE";

template<> const std::string type_identifier<type_selector<STD_BOOL>::type>::name =
		"STD_BOOL";
template<> const std::string type_identifier<type_selector<STD_INT>::type>::name =
		"STD_INT";
template<> const std::string type_identifier<type_selector<STD_DOUBLE>::type>::name =
		"STD_DOUBLE";
template<> const std::string type_identifier<type_selector<STD_LONG_DOUBLE>::type>::name =
		"STD_LONG_DOUBLE";
template<> const std::string type_identifier<type_selector<STD_FLOAT128>::type>::name =
		"STD_FLOAT128";
template<> const std::string
		type_identifier<type_selector<GMP_RATIONAL>::type>::name =
				"GMP_RATIONAL";
template<> const std::string type_identifier<type_selector<CALC_STR>::type>::name =
		"CALC_STR";

std::string coefficient_type_summary()
{
//	  std::cout << "int        : " << measure_precision<int>() << std::endl;
//	  std::cout << "float      : " << measure_precision<float>() << std::endl;
//	  std::cout << "double     : " << measure_precision<double>() << std::endl;
//	  std::cout << "long double: " << measure_precision<long double>() << std::endl;
//	  std::cout << "__float128 : " << measure_precision<__float128>() << std::endl;

    int float_precision = measure_precision<float_type> ();
	int precise_float_precision = measure_precision<precise_float_type> ();

	return to_string(float_precision) + "-bit float, " + to_string(
			precise_float_precision) + "-bit precise float";
}

}

