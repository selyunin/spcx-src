/*
 * omega_model_factory.h
 *
 *  Created on: Jan 19, 2011
 *      Author: frehse
 */

#ifndef OMEGA_MODEL_FACTORY_H_
#define OMEGA_MODEL_FACTORY_H_

#include "omega_forward.h"
#include "omega_interpol.h"
#include "omega_interpol_fb.h"

namespace continuous {
namespace support_function {

/** Factory for creating omega models */

class omega_model_factory {
public:
	template<typename scalar>
	static omega_model<scalar>* create(
			const typename omega_model<scalar>::matrix_type& A,
			const scalar& delta, const support_function_provider* X0,
			const support_function_provider* U) {
		omega_model<scalar>* inst;
		switch (active_model) {
		case FORW:
			inst = new omega_forward<scalar> (A, delta, X0, U);
			break;
		case INTERP:
			inst = new omega_interpol<scalar> (A, delta, X0, U);
			break;
		case INTERPFB:
			inst = new omega_interpol_fb<scalar> (A, delta, X0, U);
			break;
		case INTERPFBMIN:
			inst = new omega_interpol_fb<scalar> (A, delta, X0, U, true);
			break;
		default:
			throw std::runtime_error("unknown omega_model");
		}
		return inst;
	}

	static void set_model(std::string s) {
		if (s=="forw") {
			active_model = FORW;
		} else if (s=="interp") {
			active_model = INTERP;
		} else if (s=="interpfb") {
			active_model = INTERPFB;
		} else if (s=="interpfbmin") {
			active_model = INTERPFBMIN;
		} else {
			throw std::runtime_error("unknown choice of error-model '"+s+"'");
		}
	}

private:
	omega_model_factory(const omega_model_factory&);
	omega_model_factory& operator =(const omega_model_factory& other);

	enum model_type {
		FORW, INTERP, INTERPFB,INTERPFBMIN
	};
	static model_type active_model;
};

}
}

#endif /* OMEGA_MODEL_FACTORY_H_ */
