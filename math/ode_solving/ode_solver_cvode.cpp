/*
 * ode_solver_cvode.cpp
 *
 *  Created on: Oct 7, 2010
 *      Author: frehse
 */

#include "ode_solver_cvode.h"

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

namespace math {
namespace ode {

/** Forward declarations of impl functions */
/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y).
 */

static int cvodes_f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int cvodes_g(realtype t, N_Vector y, realtype *gout, void *user_data);
static int cvodes_jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int check_flag(void *flagvalue, std::string funcname, int opt);

struct ode_solver_cvode<sundial::realtype>::impl_data {
	impl_data() :
		NEQ(0), reltol(0), abstol(0), t(0), tout(0), y(0), cvode_mem(0),
				flag(0), flagr(0), iout(0), rootsfound(0), nrtfn(0) {
	}
	;
	~impl_data() {
		/* Free y vector */
		if (y)
			N_VDestroy_Serial(y);
		if (cvode_mem)
			CVodeFree(&cvode_mem);
		if (rootsfound)
			delete[] rootsfound;
	}
	;
	long int NEQ; // dimension of x
	realtype reltol, abstol, t, tout;
	N_Vector y; // current state
	void *cvode_mem;
	int flag, flagr, iout;
	int* rootsfound; // beginning of rootsfound array
	unsigned int nrtfn; // dimension of g(x)
	const state_functor* f;
	const state_vector_functor* g;
	const state_matrix_functor* j;
	positional_vdomain dom;
};

ode_solver_cvode<sundial::realtype>::ode_solver_cvode() :
	my_data(new impl_data) {
}

ode_solver_cvode<sundial::realtype>::~ode_solver_cvode() {
	// Delete all data
	delete my_data;
}

void error_message(std::string msg) {
	throw std::runtime_error("ode_solver_cvode<sundial::realtype>: error "
			+ msg);
}

/** Create a new vector with the contents of the N_Vector y. */
vector<sundial::realtype> N_Vector_to_vector(N_Vector y) {
	unsigned int n = NV_LENGTH_S(y);
	ode_solver_cvode<sundial::realtype>::vector xvec(n);
	for (unsigned int i = 0; i < n; ++i) {
		xvec[i] = NV_Ith_S(y, i);
	}
	return xvec;
}

/** Assign the contents of the vector x to the N_Vector y.
 *
 * Note that y must be previously created elsewhere. */
void assign_to_N_Vector(N_Vector y, const vector<sundial::realtype>& x) {
	unsigned int n = x.size();
	for (unsigned int i = 0; i < n; ++i) {
		NV_Ith_S(y, i) = x[i];
	}
}

void assign_to_DlsMat(DlsMat Jac, const matrix<sundial::realtype>& M){

	//number of colums
	unsigned int m=M.size2();
	unsigned int n=M.size1();
	for(int j=0;j<m;j++){
		sundial::realtype * col = DENSE_COL(Jac,j);
		for(int i=0;i<n;i++){
			col[i] = M(i,j);
		}
	}

}

void ode_solver_cvode<sundial::realtype>::init_ivp_impl() {
	/* Destroy existing data if present.
	 */
	if (my_data) {
		delete my_data;
		my_data = new impl_data;
	}

	/* Define the problem data */
	my_data->t = my_t0;
	my_data->f = my_f;
	my_data->dom = my_x0.domain();

	/* Get the number of states */
	unsigned int n = my_x0.size();
	my_data->NEQ = n;
	/* Get tolerances */
	my_data->reltol = my_odepar.rel_tol;
	my_data->abstol = my_odepar.abs_tol;

	/* Create serial vector of length NEQ for I.C. */
	my_data->y = N_VNew_Serial(my_data->NEQ);
	if (check_flag((void *) my_data->y, "N_VNew_Serial", 0))
		error_message("allocating y");

	/* Initialize y */
	assign_to_N_Vector(my_data->y, my_x0.get_vector());

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	my_data->cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *) my_data->cvode_mem, "CVodeCreate", 0))
		error_message("creating solver memory, BDF and Newton");

	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function in y'=f(t,y), the inital time T0, and
	 * the initial dependent variable vector y. */
	my_data->flag = CVodeInit(my_data->cvode_mem, cvodes_f, my_t0, my_data->y);
	if (check_flag(&my_data->flag, "CVodeInit", 1))
		error_message("initialize the integrator memory");

	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	my_data->flag = CVodeSStolerances(my_data->cvode_mem, my_data->reltol,
			my_data->abstol);
	if (check_flag(&my_data->flag, "CVodeSVtolerances", 1))
		error_message("setting tolerances");

	//	/* Call CVodeRootInit to specify the root function g with 2 components */
	//	flag = CVodeRootInit(cvode_mem, 2, cvodes_g);
	//	if (check_flag(&flag, "CVodeRootInit", 1))
	//		return (1);

	/* Call CVDense to specify the CVDENSE dense linear solver */
	my_data->flag = CVDense(my_data->cvode_mem, my_data->NEQ);
	if (check_flag(&my_data->flag, "CVDense", 1))
		error_message("specify the CVDENSE dense linear solver");

	//	/* Set the Jacobian routine to Jac (user-supplied) */
	//	flag = CVDlsSetDenseJacFn(cvode_mem, cvodes_jac);
	//	if (check_flag(&flag, "CVDlsSetDenseJacFn", 1))
	//		return (1);

	/* Call CVodeSetUserData to specify user data that is passed on to f and g */
	my_data->flag = CVodeSetUserData(my_data->cvode_mem, (void*) my_data);
	if (check_flag(&my_data->flag, "CVodeSetUserData", 0))
		error_message("setting user data");

	/* Call CVodeSetMaxStep to specify max. time step. */
	my_data->flag = CVodeSetMaxStep(my_data->cvode_mem, std::max(0.0,
			my_odepar.max_timestep));
	if (check_flag(&my_data->flag, "CVodeSetMaxStep", 0))
		error_message("setting max time step");

	/* Set number of iterations to 0. */
	my_data->iout = 0;
}

ode_solver_cvode<sundial::realtype>::ivp_result ode_solver_cvode<
		sundial::realtype>::step_ivp(sundial::realtype tfinal) {
	//std::cout << "solving from "<< my_data->t << " to "<< tfinal;
	my_data->flag = CVode(my_data->cvode_mem, tfinal, my_data->y,
			&(my_data->t), CV_NORMAL);
	//std::cout << " ended at "<< my_data->t << std::endl;

	/* Copy new state to result */
	ivp_result res;
	//	PrintOutput(t, Ith(y,1), Ith(y,2), Ith(y,3));

	//	if (my_data->flag == CV_ROOT_RETURN) {
	//		flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
	//		if (check_flag(&flagr, "CVodeGetRootInfo", 1))
	//			return (1);
	//		PrintRootInfo(rootsfound[0], rootsfound[1]);
	//	}

	if (check_flag(&my_data->flag, "CVode", 1))
		error_message("CVode solver call");
	if (my_data->flag == CV_SUCCESS) {
		res.status = ivp_result::COMPLETED;
		my_data->iout++;
	} else if (my_data->flag == CV_ROOT_RETURN) {
		res.status = ivp_result::INTERRUPTED;
		my_data->iout++;
	} else {
		res.status = ivp_result::ERROR;
	}
	res.stop_time = my_data->t;
	unsigned int n = my_data->NEQ;
	vector xvec = N_Vector_to_vector(my_data->y);
	res.stop_state = state(my_data->dom, xvec);

	res.traj
			= trajectory<sundial::realtype> (my_data->dom, res.stop_time, xvec);

	return res;
}

void ode_solver_cvode<sundial::realtype>::init_rootf_impl() {
	assert(my_g);
	init_ivp_impl();

	my_data->g = my_g;

	// call g with x0 to get the codomain
	vector g_of_x0 = my_g->map(my_x0);
	my_data->nrtfn = g_of_x0.size();
	my_data->rootsfound = new int[my_data->nrtfn];

	/* Call CVodeRootInit to specify the root function g with nrtfn components */
	my_data->flag = CVodeRootInit(my_data->cvode_mem, my_data->nrtfn, cvodes_g);
	if (check_flag(&my_data->flag, "CVodeRootInit", 1))
		error_message("defining root functions with CVodeRootInit");
}

ode_solver_cvode<sundial::realtype>::rootf_result ode_solver_cvode<
		sundial::realtype>::step_rootf(sundial::realtype tfinal) {

	ivp_result ivp_res = step_ivp(tfinal);
	rootf_result res(ivp_res);

	// Retrieve infos on roots
	unsigned int gdim = my_data->nrtfn;
	res.root_info = math::matrix<rootf_result::root_status> (1, gdim);

	if (my_data->flag == CV_ROOT_RETURN) {
		my_data->flagr = CVodeGetRootInfo(my_data->cvode_mem,
				my_data->rootsfound);
		if (check_flag(&my_data->flagr, "CVodeGetRootInfo", 1))
			error_message("getting info on roots");
		for (unsigned int i = 0; i < gdim; ++i) {
			if (my_data->rootsfound[i] > 0)
				res.root_info(0, i) = rootf_result::FROM_BELOW;
			else if (my_data->rootsfound[i] < 0)
				res.root_info(0, i) = rootf_result::FROM_ABOVE;
			else
				res.root_info(0, i) = rootf_result::NO_ROOT;
		}
	} else {
		for (unsigned int i = 0; i < gdim; ++i) {
			res.root_info(0, i) = rootf_result::NO_ROOT;
		}
	}

	return res;
}

void ode_solver_cvode<sundial::realtype>::change_rootf(const state_vector_functor& g)
{
	assert(my_g);
	this->my_g=&g;

	//cut-paste of init_rootf_impl
		my_data->g = my_g;

		// call g with x0 to get the codomain
		vector g_of_x0 = my_g->map(my_x0);
		my_data->nrtfn = g_of_x0.size();
		my_data->rootsfound = new int[my_data->nrtfn];

		/* Call CVodeRootInit to specify the root function g with nrtfn components */
		my_data->flag = CVodeRootInit(my_data->cvode_mem, my_data->nrtfn, cvodes_g);
		if (check_flag(&my_data->flag, "CVodeRootInit", 1))
			error_message("defining root functions with CVodeRootInit");

}
void ode_solver_cvode<sundial::realtype>::specify_jacobian(const state_matrix_functor & j){

	this->my_jac = &j;
	my_data->j=my_jac;

	my_data->flag = CVDlsSetDenseJacFn(my_data->cvode_mem, cvodes_jac);
	if (check_flag(&my_data->flag, "CVDlsSetDenseJacFn", 0))
		error_message("specifying Jacobian function with CVDlsSetDenseJacFn");
}

 sundial::realtype ode_solver_cvode<sundial::realtype>::get_solver_current_time() const{

	sundial::realtype tcur;
	my_data->flag=CVodeGetCurrentTime(my_data->cvode_mem,&tcur);
	if (check_flag(&my_data->flag, "CVodeGetCurrentTime", 0))
				error_message("getting current integration time with CVodeGetCurrentTime");
	return tcur;
}

std::string ode_solver_cvode<sundial::realtype>::get_solver_name() const {
	return "CVODE";
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y).
 */

static int cvodes_f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
	ode_solver_cvode<sundial::realtype>::impl_data* my_data =
			static_cast<ode_solver_cvode<sundial::realtype>::impl_data*>(user_data);

	//std::cout << "Entering cvodes_f..." << std::endl;

	// convert y to a state
	ode_solver_cvode<sundial::realtype>::vector xvec = N_Vector_to_vector(y);
	ode_solver_cvode<sundial::realtype>::state x(my_data->dom, xvec);

	// compute xdot=f(x)
	ode_solver_cvode<sundial::realtype>::state xdot = my_data->f->map(x);

	//std::cout << "time t:" << t << std::endl;
	//std::cout << "state x:" << x << std::endl;
	//std::cout << "deriv x':" << xdot << std::endl;


	// convert xdot to N_Vector
	assign_to_N_Vector(ydot, xdot.get_vector());

	return (0);
}

/*
 * g routine. Compute functions g_i(t,y) for i = 0,1.
 */

static int cvodes_g(realtype t, N_Vector y, realtype *gout, void *user_data) {

	ode_solver_cvode<sundial::realtype>::impl_data* my_data =
			static_cast<ode_solver_cvode<sundial::realtype>::impl_data*>(user_data);

	// convert y to a state
	ode_solver_cvode<sundial::realtype>::vector xvec = N_Vector_to_vector(y);
	ode_solver_cvode<sundial::realtype>::state x(my_data->dom, xvec);

	// compute xdot=f(x)
	ode_solver_cvode<sundial::realtype>::vector gres = my_data->g->map(x);

	// convert xdot to N_Vector
	for (unsigned int i = 0; i < my_data->nrtfn; ++i) {
		gout[i] = gres[i];
	}

	return (0);
}


/**
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */
static int cvodes_jac(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
	realtype y1, y2, y3;

	ode_solver_cvode<sundial::realtype>::impl_data* my_data =
				static_cast<ode_solver_cvode<sundial::realtype>::impl_data*> (user_data);

	ode_solver_cvode<sundial::realtype>::vector xvec = N_Vector_to_vector(y);
	ode_solver_cvode<sundial::realtype>::state x(my_data->dom, xvec);
	ode_solver_cvode<sundial::realtype>::matrix M = my_data->j->map(x);
	M.reorder(my_data->dom,my_data->dom);
	assign_to_DlsMat(J,M.get_matrix());

	return (0);
}

static int check_flag(void *flagvalue, std::string funcname, int opt) {
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		fprintf(stderr,
				"\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname.c_str());
		return (1);
	}

	/* Check if flag < 0 */
	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
					funcname.c_str(), *errflag);
			return (1);
		}
	}

	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == NULL) {
		fprintf(stderr,
				"\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
				funcname.c_str());
		return (1);
	}

	return (0);
}

}
}
