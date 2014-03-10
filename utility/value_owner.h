#ifndef VALUE_OWNER_H_
#define VALUE_OWNER_H_

/** A simple class for storing a value of type value_type.
 * It is useful as a base class, whose pointer can be passed
 * between classes that know nothing else about each other. 
 * @todo Which design pattern does this correspond to? Traits?
 * */

template<typename value_type> class value_owner {
public:
	void set(value_type new_value) {
		ret_value=new_value;
	}
	;
	value_type get() {
		return ret_value;
	}
	;
private:
	value_type ret_value;
};

#endif /*VALUE_OWNER_H_*/
