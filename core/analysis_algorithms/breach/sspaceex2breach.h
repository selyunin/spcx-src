#ifndef SSPACEEX2BREACH_H
#define SSPACEEX2BREACH_H

#include "whatever_is_needed.h"

class SSpaceExModel {

 private:
  
  /* do whatever you need or not to hide here */

 public:

  /* Constructor */

  SSpaceExModel(string ModelName);  // parse the file ModelName and create an automaton

  /* Methods */ 

  loc_id get_loc_id(string loc_name);
  string  get_loc_name(loc_id lid);
  
  param_id get_param_id(string param_name);
  string  get_param_name(param_id pid);

  int get_dimx();  // dimension of the state vector
  int get_dimp();  // dimension of the parameter; initial state is part of p so dimp >= dimx 
  int get_dimg();  // number of guards

  /* init function: set initial state, location and parameters; if no parameter, p is x0 */

  void init(const double_vector &p, loc_id location);

  /* f function: compute the dynamics xdot at state x and current location */

  void f(const double_vector &x, double_vector &xdot); // dynamics at current state x

  /* g (guard) function: function; gout should contain (some) distances to *all* guards */
  /*                               if a guard g_i cannot be reached in the current   */
  /*                               location, set g_i = inf or some big value         */
  
  void g(const double_vector &x, double_vector &gout);

  /* r (reset) function: if x is in a guard take a transition and update the state and the location */
  /*                     otherwise, do nothing.                                                     */
   
  void r(const double_vector &x, double_vector &new_x, loc_id &new_location);

};

#endif
