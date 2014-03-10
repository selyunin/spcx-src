#ifndef ID_TO_SHARED_PTR_CACHE_H_
#define ID_TO_SHARED_PTR_CACHE_H_

#include <map>
#include <boost/shared_ptr.hpp>

/** A class for associating "ids" with some object of type \p T.
 * ids effectively act like pointers, but have the advantage that they
 * can remain unchanged when cloning containers of ids.
 */
template <typename T> class id_to_object_cache {
public:
	typedef unsigned int id_type;
	typedef std::map<id_type,T> id_to_object_map_type;

	id_to_object_cache() :
		my_highest_id(0) {
	}
	;
	virtual ~id_to_object_cache() {
	}
	;

	/** Insert a copy of object \p p into the cache, and return a new id that is associated to it. */
	id_type insert(const T& p) {
		id_type new_id=get_new_id();
		my_map.insert(std::make_pair(new_id, p));
		return new_id;
	}
	;

	/** Return the object associated with \p id. If none is found, return a T(). */
	T get(id_type id) const {
		typename id_to_object_map_type::const_iterator it= my_map.find(id);
		if (it!=my_map.end())
			return it->second;
		else
			return T();
	}
	;
private:
	id_type get_new_id() {
		return ++my_highest_id;
	}
	id_to_object_map_type my_map;
	id_type my_highest_id;
};

/** A class for associating "ids" with a shared pointer to some object.
 * ids effectively act like pointers, but have the advantage that they
 * can remain unchanged when cloning containers of ids.
 */
template <typename T> class id_to_shared_ptr_cache {
public:
	typedef unsigned int id_type;
	typedef boost::shared_ptr<T> ptr;
	typedef std::map<id_type,ptr> id_to_ptr_map_type;

	id_to_shared_ptr_cache() :
		my_highest_id(0) {
	}
	;
	virtual ~id_to_shared_ptr_cache() {
	}
	;

	/** Insert the new pointer \p p into the cache, and return a new id that is associated to it. */
	id_type insert(const ptr& p) {
		id_type new_id=get_new_id();
		my_map.insert(std::make_pair(new_id, p));
		return new_id;
	}
	;

	/** Return the pointer associated with \p id. If none is found, return a shared_ptr<T>() (equivalent to a null pointer). */
	ptr get_ptr(id_type id) const {
		typename id_to_ptr_map_type::const_iterator it= my_map.find(id);
		if (it!=my_map.end())
			return it->second;
		else
			return ptr();
	}
	;
	
	typedef typename id_to_ptr_map_type::const_iterator const_iterator;
	virtual const_iterator begin() const { return my_map.begin(); };
	virtual const_iterator end() const { return my_map.end(); };	
	virtual typename id_to_ptr_map_type::size_type size() const { return my_map.size(); };

private:
	id_type get_new_id() {
		return ++my_highest_id;
	}
	id_to_ptr_map_type my_map;
	id_type my_highest_id;
};

#endif /*ID_TO_SHARED_PTR_CACHE_H_*/
