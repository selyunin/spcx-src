#ifndef VALUATION_H_
#define VALUATION_H_

//#include "global/global_types.h" // for __float128 operator<<

#include <map>
#include <iostream>
#include "boost/shared_ptr.hpp"
#include "utility/stl_helper_functions.h"
#include "math/vdom/variable.h"
#include "utility/shared_ptr_user.h"

namespace valuation_functions {

//using ::operator<<;

// Note: I tried to declare operator<< as a friend, but couldn't get it to compile.
// It couldn't make the association with the operator<< in stl_helper_functions to print the map.

// forward declaration needed, otherwise friend declaration below gives an error
// see http://www.parashift.com/c++-faq-lite/templates.html#faq-35.15
// and http://bytes.com/forum/thread63289.html
//template<typename T> class valuation;
//template<typename T> std::ostream& operator<<(std::ostream& os, const valuation<T>& v);

/** A valuation attributes a value of type \p T to a number of variables.
 * Variables are considered to be of type variable, i.e., accessible via
 * their name (a string) or an identifier of type \p variable_id.
 *
 * \todo Provide a template member function for casting to different types.
 */
template<typename T, typename codomain_type=variable_id> class valuation {
public:
	typedef typename boost::shared_ptr<valuation<T,codomain_type> > ptr;
	typedef typename boost::shared_ptr<const valuation<T,codomain_type> > const_ptr;

	valuation(const valuation& v){
		my_valuation_map = v.get_valuation_map();
	}
	;

	valuation(){
	}
	;

	virtual ~valuation(){
	}
	;


	/** Attribute the value val to the variable with identifier id. */
	virtual void add(const codomain_type& id, const T& val) {
		my_valuation_map[id]=val;
	}
	;

	/** Attribute the value val to the codomain_type with name n.*/
	//virtual void add(const std::string& n, const T& val) = 0;

	/** Returns the value of the variable with identifier id.
	 * Throws an exception if the variable has no value attributed to it.*/
	virtual const T& get_value(const codomain_type& id) const {
		// need to include typename keyword, see http://www.parashift.com/c++-faq-lite/templates.html#faq-35.18
		typename valuation_map_t::const_iterator pos=my_valuation_map.find(id);
		if (pos==my_valuation_map.end()) {
			print(std::cerr);
			std::cerr << std::endl << std::flush;
			throw std::runtime_error("Unknown codomain_type "+int2string(id)+".");
		}
		return pos->second;
	}
	;

	/** Returns the value of the codomain_type with name n.*/
	//virtual const T& get_value(const std::string& n) const = 0;

	/** Returns true if the valuation has a value attributed to the variable with id \p vid.*/
	virtual bool has_value(codomain_type id) const {
		// need to include typename keyword, see http://www.parashift.com/c++-faq-lite/templates.html#faq-35.18
		typename valuation_map_t::const_iterator pos=my_valuation_map.find(id);
		return !(pos==my_valuation_map.end());
	}
	;

	// \see chapter C.13.2 of Stroustrup, TC++PL3 on why <> is necessary
	// friend std::ostream& operator<< <>(std::ostream& os, const valuation<T>& v);

	/** Output all variables with their values in the form:
	 * \code [x->1,y->2,z->3] \endcode
	 */
	virtual void print (std::ostream& os) const { 	os << "[";
	for (typename valuation_map_t::const_iterator it=my_valuation_map.begin(); it!=my_valuation_map.end(); ++it)
	{
		if (it!=my_valuation_map.begin()) os << ",";
		os << it->first << "->" << it->second;
	}
	os << "]"; };

	virtual const std::map<codomain_type, T>& get_valuation_map() const{
		return my_valuation_map;
	}

protected:
	typedef std::map<codomain_type, T> valuation_map_t;
	valuation_map_t my_valuation_map;

};

template<typename T> class variable_valuation : public valuation<T, variable_id> {
public:
	typedef typename boost::shared_ptr<variable_valuation<T> > ptr;
	typedef typename boost::shared_ptr<const variable_valuation<T> > const_ptr;

	variable_valuation(const valuation<T, variable_id>& v) : valuation<T, variable_id>(v){
	}
	;

	variable_valuation(){
	}
	;

	virtual ~variable_valuation(){
	}
	;

	/** The following using directives are needed to redirect name look up
	 * to the base class. If left out, the base class methods are simply not
	 * seen by the compiler.
	 * (Welcome to the happy of world of template programming)
	 */
	using valuation<T, variable_id>::add;
	using valuation<T, variable_id>::get_value;

	/** Attribute the value val to the variable with name n.*/
	virtual void add(const std::string& n, const T& val) {
		valuation<T, variable_id>::add(variable::get_or_add_variable_id(n), val);
	}
	;

	/** Returns the value of the variable with name n.
	 * Throws an exception if the variable has no value attributed to it. */
	virtual const T& get_value(const std::string& n) const {
		return valuation<T, variable_id>::get_value(variable::get_or_add_variable_id(n));
	}
	;

	/** Output all variables with their values in the form:
	 * \code [x->1,y->2,z->3] \endcode
	 */
	virtual void print (std::ostream& os) const { 	os << "[";
	for (typename valuation<T, variable_id>::valuation_map_t::const_iterator it=this->my_valuation_map.begin(); it!=this->my_valuation_map.end(); ++it)
	{
		if (it!=this->my_valuation_map.begin()) os << ",";
		os << variable(it->first) << "->" << it->second;
	}
	os << "]"; };
};

}

/** Output all variables with their values in the form:
 * \code [x->1,y->2,z->3] \endcode
 */
template<typename T> std::ostream& operator<<(std::ostream& os,
		const valuation_functions::valuation<T>& v) {
//	os << v.my_valuation_map; // I have no idea why this doesn't work
	v.print(os);
	return os;
}

#endif /*VALUATION_H_*/
