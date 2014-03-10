/*
 * aaflib_loader.h
 *
 *  Created on: Mar 21, 2013
 *      Author: notroot
 */

#ifndef AAFLIB_LOADER_H_
#define AAFLIB_LOADER_H_

/** This header file includes the aaflib and defines various extensions */

#include "aaflib/aa.h"

/** Exponentiation
 *
 * This is a forwarding call to a function
 * with non-standard name.
 */
inline
AAF pow( const AAF& x, const AAF& y) {
	// try to convert to int if possible
	if (y.rad()==0.0) {
		double y_double = y.getcenter();
		double intpart;
		double fractpart = modf(y_double,&intpart);
		if (fractpart == 0.0 && (int)y_double == y_double) {
			return aaf_pow(x,(int)y_double);
		}
	}
	return aaf_pow(x,y);
};

/** Product
 *
 * This is a forwarding call to the member operator that is nonconst */
inline
AAF operator*(const AAF& x, const AAF& y) {
	AAF z(x);
	return z.operator*(y);
}
;

/** Element Conversion */
template<>
inline
AAF convert_element<AAF,double>(
		const double& x) {
	return AAF(x);
}
;


#endif /* AAFLIB_LOADER_H_ */
