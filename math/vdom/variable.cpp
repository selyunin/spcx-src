#include "math/vdom/variable.h"

#include <limits>
#include <iostream>
#include <iostream>
#include <stdexcept>
#include <assert.h>

#include "utility/basic_exception.h"
#include "utility/stl_helper_functions.h"

using namespace std;
std::map<string, variable_id> variable::string_to_variable_id_cache;
std::vector<string> variable::variable_id_to_string_cache;
std::vector<unsigned int> variable::variable_id_to_dim_cache;
variable_id variable::primed_factor = (1 << 20); // max nb of variables 2^(20-8)-1
variable_id variable::vector_factor = (1 << 8); // max vector size 255
// max. number of variables = 2^(20-8)-1 = 4095
variable_id variable::highest_id = variable::vector_factor-1;
char variable::prime_char = '\'';

variable::variable() : my_id(0) {
}

variable::variable(std::string n, unsigned int dim) {
	my_id = get_or_add_variable_id(n,dim);
}

variable::variable(variable_id id) : my_id(id) {
}

const variable_id& variable::get_id() const {
	return my_id;
}

variable_id variable::get_vector_id() const {
	return compute_vector_id(my_id);
}

bool variable::is_vector() const {
	return compute_vector_id(my_id)==my_id;
}

std::string variable::get_name() const {
	return get_name(my_id);
}

unsigned int variable::get_dimension() const {
	return get_dimension(my_id);
}

variable_id variable::get_element(unsigned int k) const {
	if (is_vector()) {
		if (k > 0 && k <= get_dimension())
			return compute_element_id(my_id,k);
		else
			throw std::runtime_error(
					"variable::get_element: Index out of bounds.");
	} else
		throw std::runtime_error("variable::get_element: not a vector.");
}

unsigned int variable::get_element_pos() const {
	return get_element_pos(my_id);
}

bool variable::operator<(const variable& x) const {
	return my_id<x.my_id;
}

bool variable::operator==(const variable& x) const {
	return my_id==x.my_id;
}

bool variable::operator!=(const variable& x) const {
	return my_id!=x.my_id;
}

variable_id variable::get_or_add_variable_id(std::string new_name, unsigned int dim) {
	return get_or_add_variable_id_primed(new_name, 0, dim);
}

variable_id variable::get_variable_id(std::string var_name)
{
	if(string_to_variable_id_cache.find(var_name)!=string_to_variable_id_cache.end())
		return string_to_variable_id_cache.find(var_name)->second;
	else {
		std::stringstream ss;
		print_variable_cache(ss);
		throw std::runtime_error(" variable::get_variable_id: unknown variable \"" + var_name + "\". Known variables:\n"+ss.str());
	}
}

variable_id variable::get_temp_variable_id() {
	// take as string the end of the map and add "z"
	if (string_to_variable_id_cache.empty())
		return get_or_add_variable_id("zzzzztemp");
	else
		return get_or_add_variable_id(string_to_variable_id_cache.rbegin()->first + "z");
}

void variable::remove_variable_id(variable_id id) {
	string n = variable_id_to_string_cache[cache_id(id)];
	string_to_variable_id_cache.erase(n);
	variable_id_to_string_cache.erase(variable_id_to_string_cache.begin() + cache_id(id));
	variable_id_to_dim_cache.erase(variable_id_to_dim_cache.begin() + cache_id(id));
}

void variable::get_unprimed_name_and_add_prime_count(std::string& new_name, unsigned int& prime_count) {
	string::size_type prime_pos = new_name.find(variable::prime_char, 0);
	if (prime_pos != string::npos) {
		unsigned int prime_nr = 0;
		while (prime_pos != string::npos) {
			prime_pos = new_name.find(variable::prime_char, prime_pos + 1);
			++prime_nr;
		}
		prime_count += prime_nr;
		// Take as name the substring up to the first occurrence of '
		new_name = new_name.substr(0, new_name.find(variable::prime_char, 0));
	}
}

bool variable::has_variable(const std::string name) {
	return string_to_variable_id_cache.find(name)!=string_to_variable_id_cache.end();
}

const variable::name_store_type& variable::get_names() {
	return variable_id_to_string_cache;
}

void variable::get_vector_name(std::string& name, unsigned int& element_index) {
	element_index=0;
	string::size_type open_pos = name.find('(', 0);
	string::size_type close_pos = name.find(')', 0);

	//std::cout << "getting:" << name <<"." << std::endl;
	if (open_pos != string::npos && close_pos != string::npos) {
		element_index = from_string<unsigned int>(name.substr(open_pos+1,close_pos-open_pos-1));
		name=name.substr(0,open_pos);
	}
	//std::cout << "name:"<<name << ". idx:"<< element_index << std::endl;
}

variable_id variable::get_or_add_variable_id_primed(std::string name,
		unsigned int prime_count, unsigned int dim) {
	std::string new_name=name;
	variable_id res;
	// Get the primedness of new_name and add it to prime_count
	get_unprimed_name_and_add_prime_count(new_name, prime_count);

	unsigned int element_index=0;
	get_vector_name(new_name,element_index);

	std::map<string, variable_id>::const_iterator pos =
			string_to_variable_id_cache.find(new_name);

	if (pos == string_to_variable_id_cache.end()) {
		// can't create a vector element without before creating the vector
		if (element_index>0) {
			throw std::runtime_error("Can't find vector element variable "
					+ name + " because vector variable " + new_name
					+ " is unknown.");
		}
		// create a new id
		variable_id new_id = get_new_id();
		string_to_variable_id_cache.insert(make_pair(new_name, new_id));
		variable_id cid=cache_id(new_id);
//std::cerr << "cid:" << cid << "bef:" << variable_id_to_string_cache.size() << std::endl;
		if (cid>=variable_id_to_string_cache.size()) {
			variable_id_to_string_cache.resize(cid+1);
//			std::cerr << "resize to " << cid+1 << " result " << variable_id_to_string_cache.size() << std::endl;
		}
		variable_id_to_string_cache[cid] = new_name;
		if (dim + 1 >= variable::vector_factor)
					throw std::out_of_range("vector dimension " + to_string(dim)
					+ " higher than limit (" + to_string(
					variable::vector_factor - 1) + ")");
		if (cid>=variable_id_to_dim_cache.size())
			variable_id_to_dim_cache.resize(cid+1);
		variable_id_to_dim_cache[cid] = dim;
		res = get_primed_id(new_id, prime_count);
	} else // already in the cache
	{
		variable_id vid=pos->second; // vector/scalar id
		variable_id cid=cache_id(vid);
		// If it's an element id, then the index has to be within the vector dimension
		if (element_index > 0) {
			if (dim > 1) {
				throw std::runtime_error("vector element "+name+" needs to be a scalar");
			}
			if (element_index > get_dimension(vid)) {
				throw std::out_of_range(
						"requested vector element index "+name+" exceeds dimension");
			}
			res = get_primed_id(compute_element_id(vid,element_index), prime_count);
		} else {
			// It's a vector or scalar. Check dimension with what's already in the cache.
			if (get_dimension(vid) != dim) {
				throw std::out_of_range(
						"requested adding variable "+name+" with dimension "+to_string(dim)+" to cache but it already exists with different dimension "+to_string(get_dimension(vid)));
			}
			res = get_primed_id(vid, prime_count);
		}
	}
	return res;
}

unsigned int variable::max_dim() {
	return vector_factor-1;
}

string variable::get_name(variable_id id) {
	// we use the fact that ids are given contiguously:
	// id is in the cache iff id<get_variable_count()
	variable_id uid = get_primed_id(id, 0);
	if (uid > highest_id)
		throw std::out_of_range("id not in cache of variable");
	unsigned int pos = get_element_pos(id);
	if (pos > 0) {
		return variable_id_to_string_cache[cache_id(uid)] + '(' + int2string(pos) + ')' + string(get_prime_count(id),
				variable::prime_char);
	} else {
		return variable_id_to_string_cache[cache_id(uid)] + string(get_prime_count(id),
				variable::prime_char);
	}
}

variable_id variable::get_id_primedness_increased(variable_id id, unsigned int prime_count) {
	return get_primed_id(id, get_prime_count(id) + prime_count);
}

variable_id variable::get_id_primedness_decreased(variable_id id, unsigned int prime_count) {
	unsigned int new_prime_count=get_prime_count(id);
	if (new_prime_count<prime_count) new_prime_count=0;
	else
		new_prime_count-=prime_count;
	return get_primed_id(id, new_prime_count);
}

variable_id variable::get_primed_id(variable_id id, unsigned int prime_count) {
	if (prime_count > std::numeric_limits<variable_id>::max()
			/ variable::primed_factor)
		throw basic_exception("prime count exceeded allocatable limit");
	return (id % variable::primed_factor) + variable::primed_factor * prime_count;
}

unsigned int variable::get_prime_count(variable_id id) {
	return id / variable::primed_factor;
}

variable_id variable::cache_id(variable_id id) {
	return (id % variable::primed_factor) / variable::vector_factor;
}

variable_id variable::cache_to_variable_id(variable_id id) {
	return id * variable::vector_factor;
}

variable_id variable::compute_vector_id(variable_id id) {
	variable_id cid=cache_id(id);
	// strip vector element bits without removing the primedness bits
	if (variable_id_to_dim_cache[cid] > 0)
		return (id / variable::vector_factor) * variable::vector_factor;
	else
		return 0;
}

variable_id variable::compute_element_id(variable_id id, unsigned int k) {
	assert(compute_vector_id(id)>0); // it has to be a vector
	return compute_vector_id(id)+k;
}

unsigned int variable::get_dimension(variable_id id) {
	variable_id cid=cache_id(id);
	// strip vector element bits
	return variable_id_to_dim_cache[cid];
}

unsigned int variable::get_element_pos(variable_id id) {
	variable_id vid=compute_vector_id(id);
	if (vid>0) {
		return id-compute_vector_id(id);
	} else
		return 0;
}


void variable::print_variable_cache(std::ostream& o) {
	for (std::map<std::string, variable_id>::const_iterator it =
			string_to_variable_id_cache.begin(); it != string_to_variable_id_cache.end(); ++it) {
		o << it->first << " <-> " << it->second << endl;
	}
}

variable_id variable::get_new_id() {
	highest_id+=vector_factor;
	if (highest_id+vector_factor>=primed_factor)
		throw basic_exception("max. number of allocatable variables exceeded");
	return (highest_id / variable::vector_factor) * variable::vector_factor;
}

void variable::reset() {
	string_to_variable_id_cache.clear();
	variable_id_to_string_cache.clear();
	variable_id_to_dim_cache.clear();
	highest_id = vector_factor-1;
}

variable_id_set get_unprimed_variables(const variable_id_set& vis) {
	variable_id_set ret;
	variable_id_set::iterator jt = ret.begin();
	for (variable_id_set::const_iterator it = vis.begin(); it != vis.end(); ++it) {
		ret.insert(jt, variable::get_primed_id(*it, 0));
	}
	return ret;
}

variable_id_set get_primed_variables(const variable_id_set& vis) {
	variable_id_set ret;
	variable_id_set::iterator jt = ret.begin();
	for (variable_id_set::const_iterator it = vis.begin(); it != vis.end(); ++it) {
		if (variable::get_prime_count(*it) > 0) {
			ret.insert(jt, *it);
		}
	}
	return ret;
}

void print_variable_id_set(std::ostream& os,const variable_id_set& v) {
	for (variable_id_set::const_iterator it = v.begin(); it != v.end(); ++it) {
		if (it != v.begin())
			os << ",";
		os << variable(*it);
	}
}

void throw_if_not_contains(const variable_id_set& v1,
		const variable_id_set& v2, std::string prologue, std::string epilogue) {
//	std::cout << "testing: ";
//	print_variable_id_set(std::cout,v1);
//	std::cout << " in ";
//	print_variable_id_set(std::cout,v2);
//	std::cout << "?" << std::endl;

	if (!set_contains(v1, v2)) {
		variable_id_set vis(v2);
		set_difference_assign(vis, v1);
		std::stringstream ss("");
		print_variable_id_set(ss,vis);
		throw basic_exception(prologue + ss.str() + epilogue);
	}
}

