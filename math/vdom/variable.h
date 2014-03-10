#ifndef variable_H_
#define variable_H_

#include <set>
#include <list>
#include <vector>
#include <map>
#include <string>
#include "variable_formatter.h"

typedef unsigned int variable_id;
typedef std::set<variable_id> variable_id_set;
typedef std::list<variable_id> variable_id_list;

/** A class for attributing variable names (string) to identifiers of type \p variable_id.
 * The association is cached in a bidirectional map.
 *
 * Variable can be "primed", i.e., one can obtain the id corresponding to the variable x' from x (and vice versa)
 * via a simple calculation. This facilitates dealing with relations etc., which can be seen as continuous sets over
 * the variables and the primed variables, e.g., the identity being \{ (x,x') | x=x' \}.
 * To obtain the image (or preimage) of the relation, one needs to project onto the primed (unprimed) variables.
 * For more complex computations, variables can be primed to an arbitrary degree \p prime_count.
 *
 * Variables can be vectors of a dimension up to variable::max_dim().
 * A scalar variable has dimension 0, a vector variable a dimension>0.
 * A vector variable has a special id to denote the entire vector, called its vector_id.
 * Each element of the vector has a different id that is different from the vector_id.
 * The names of the vector elements are "x(k)", where x is the name of the vector variable
 * and k is the index of the element. The index ranges from 1 to the dimension of the vector.
 *
 * Note that the primes of vector elements come after the parentheses, e.g., x(1)'''.
 */
class variable {
public:
	typedef std::vector<std::string> name_store_type;

	/** Variable without valid id. */
	variable();

	/** Create a variable with name n and dimension dim.
	 *
	 * If dim is omitted, the variable is a scalar.
	 *
	 * @note Variable names are unique, i.e., two variables with the same
	 * name refer the same values.
	 */
	explicit variable(std::string n, unsigned int dim = 0);

	/** Construct a variable using an existing id.
	 *
	 * @attention Does not check whether id exists already in the cache.
	 * For expert use only. */
	explicit variable(variable_id id);

	/** Get the id of the variable. */
	const variable_id& get_id() const;

	/** Get the id of the corresponding vector variable.
	 *
	 * Returns 0 if not a vector. */
	variable_id get_vector_id() const;

	/** Returns true iff the variable is a vector. */
	bool is_vector() const;

	/** Get the name of the variable. */
	std::string get_name() const;

	/** Get the dimension of the variable. */
	unsigned int get_dimension() const;

	/** Get the id of the k-th element of a vector variable.
	 *
	 * k may range from 1 to n, where n is the dimension of the variable.
	 *
	 * Applies if *this is either a vector variable,
	 * otherwise throws.*/
	variable_id get_element(unsigned int k) const;

	/** Get the index of the element of a vector variable.
	 *
	 * Returns 0 if the variable is not element of a vector
	 * (or if it's the a vector variable).
	 * .*/
	unsigned int get_element_pos() const;

	/** Less operator
	 *
	 * The comparison induces a complete order, for use with std::set etc. */
	bool operator<(const variable& x) const;

	/** Equal operator */
	bool operator==(const variable& x) const;

	/** Equal operator */
	bool operator!=(const variable& x) const;

	/** Returns true if a variable with that name is known. */
	static bool has_variable(const std::string name);

	/** Obtain a unique identifier for the variable name \p new_name with dimension dim.
	 * If the identifier is not found in the cache, then a new one is created.
	 */
	static variable_id get_or_add_variable_id(std::string new_name, unsigned int dim = 0);

	/** Obtain a unique identifier for the variable name \p new_name with dimension dim,
	 *  primed to the degree \p prime_count.
	 */
	static variable_id get_or_add_variable_id_primed(std::string new_name,
			unsigned int prime_count, unsigned int dim = 0);

	/**
	 * Returns the id of varible with name var_name
	 */
	static variable_id get_variable_id(std::string var_name);

	/** Obtain a unique identifier for a temporary variable whose id is
	 * not yet in the cache. The variable is attributed a name that is
	 * not yet in the cache.
	 */
	static variable_id get_temp_variable_id();

	/** Returns the name associated with the identifier \p id.
	 */
	static std::string get_name(variable_id id);

	/** Returns all known names. */
	static const name_store_type& get_names();

	/** Returns the id primed to the degree \p prime_count. */
	static variable_id get_primed_id(variable_id id, unsigned int prime_count = 1);

	/** Returns the id primed to its current degree plus \p prime_count.
	 */
	static variable_id get_id_primedness_increased(variable_id id, unsigned int prime_count = 1);

	/** Returns the id primed to its current degree minus \p prime_count.
	 */
	static variable_id get_id_primedness_decreased(variable_id id, unsigned int prime_count = 1);

	/** Returns the degree of primedness of the variable with identifier \p id.
	 */
	static unsigned int get_prime_count(variable_id id);

	/** Output the list of names and associated ids. */
	static void print_variable_cache(std::ostream& o);

	/** Clear the cache of all variables. */
	static void reset();

	/** Remove the variable with id \p id and its name from the cache.
	 * @attention This is very slow as all subsequent elements must be moved.
	 */
	static void remove_variable_id(variable_id id);

	/** Returns the max dimension of vector variables. */
	static unsigned int max_dim();

	/** Helper function to get the name without the primes, and add the primes in the string to the prime count. */
	static void get_unprimed_name_and_add_prime_count(std::string& new_name,
			unsigned int& prime_count);

private:
	variable_id my_id;

	/** Helper function that returns an id that is guaranteed to be new. */
	static variable_id get_new_id();
	/** Helper function to compute the corresponding vector id.
	 * Returns 0 if id is neither a vector variable nor an element of a vector. */
	static variable_id compute_vector_id(variable_id id);
	/** Helper function to compute the id of the kth element of vector variable id.
	 * k may range from 1 to n, where n is the dimension of the vector. */
	static variable_id compute_element_id(variable_id id, unsigned int k);
	/** Get the dimension of the variable id. */
	static unsigned int get_dimension(variable_id id);
	/** Get the index of the element of a vector variable.
	 *
	 * Returns 0 if the variable is not element of a vector, or if it is
	 * a vector variable.
	 * */
	static unsigned int get_element_pos(variable_id id);
	/** Helper function that returns the cache_id associated with a variable id. */
	static variable_id cache_id(variable_id id);
	/** Helper function that returns the vector or scalar id associated with a cache id. */
	static variable_id cache_to_variable_id(variable_id id);
	/** Helper function for obtaining the vector name from an element name.
	 *
	 * For example "x(4)" gives the vector name "x" and the element index 4.
	 * If the name is a vector name, the name is unchanged and the element index returned is 0.
	 */
	static void get_vector_name(std::string& name, unsigned int& element_index);

	static std::map<std::string,variable_id> string_to_variable_id_cache;
	static std::vector<std::string> variable_id_to_string_cache;
	static std::vector<unsigned int> variable_id_to_dim_cache;
	static unsigned int highest_id;
	static unsigned int primed_factor; // primeness is defined as variable_id div primed_factor, the unprimed variable_id is variable_id mod primed_factor
	static unsigned int vector_factor; // space between variable ids to account for vectors of max size vector_factor-1
	static char prime_char; // the character representing primes (')

	friend class index_to_variable_id_map_provider;
	friend class positional_vdomain;
};

/** Return the variables in vis with their unprimed ids. */
variable_id_set get_unprimed_variables(const variable_id_set& vis);

/** Return the primed variables in vis with their primed ids. */
variable_id_set get_primed_variables(const variable_id_set& vis);

/** Throw if vis1 does not contain vis2 and produce a user-readable message
 * of the form "prologuevar1,var2,var3epilogue". */
void throw_if_not_contains(const variable_id_set& v1,
		const variable_id_set& v2, std::string prologue, std::string epilogue);

void print_variable_id_set(std::ostream& os,const variable_id_set& vis);

#endif /*variable_H_*/
