/*
 * collection_base.h
 *
 *  Created on: Jun 20, 2009
 *      Author: frehse
 */

#ifndef COLLECTION_BASE_H_
#define COLLECTION_BASE_H_

#include <list>
#include <boost/shared_ptr.hpp>
#include "utility/simple_iterators/proxy_forward_iterator.h"

namespace simple_iterators {

/** A base class defining the interface of a "collection", i.e.,
 * an enumerable set whose elements are accessible via a forward const_iterator.
 *
 * In order to avoid constructing proxy iterators in begin() and end(),
 * they are stored in private members. Template functions
 * are provided for easy access.
 */
template<typename T> class collection_const_base {
public:
	virtual ~collection_const_base() {
	}
	;
	typedef T object_type;
	typedef proxy_forward_const_iterator<T> const_iterator;
	virtual const const_iterator& begin() const = 0;
	virtual const const_iterator& end() const = 0;
protected:
	template<typename implementation_iterator_type> const const_iterator& get_const_begin(
			const implementation_iterator_type& it) const {
		// test if my_begin already points to the same object
		// gf: correction : the test is illegal since either iterator might not be dereferenceable
		//if (my_const_begin.operator->() != it.operator->()) {

		// assign it to my_begin and return my_begin
			const_iterator& b=const_cast<const_iterator&>(my_const_begin);
			b.set_implementation_iterator(it);
		//}
		return my_const_begin;
	}
	;
	template<typename implementation_iterator_type> const const_iterator& get_const_end(
			const implementation_iterator_type& it) const {
		// test if my_begin already points to the same object
		// gf: correction : the test is illegal since either iterator might not be dereferenceable
		//if (my_const_end.operator->() != it.operator->()) {
			// assign it to my_begin and return my_begin
			const_iterator& b=const_cast<const_iterator&>(my_const_end);
			b.set_implementation_iterator(it);
		//}
		return my_const_end;
	}
	;
private:
	const_iterator my_const_begin;
	const_iterator my_const_end;
};

/** A base class defining the interface of a "collection", i.e.,
 * an enumerable set whose elements are accessible via a forward iterator.
 *
 * In order to avoid constructing proxy iterators in begin() and end(),
 * they are stored in private members. Template functions
 * are provided for easy access.
 */
template<typename T> class collection_nonconst_base {
public:
	virtual ~collection_nonconst_base() {
	}
	;
	//typedef T object_type;
	typedef proxy_forward_iterator<T> iterator;
	virtual const iterator& begin() = 0;
	virtual const iterator& end() = 0;
protected:
	template<typename implementation_iterator_type> const iterator& get_begin(
			const implementation_iterator_type& it) {
		// test if my_begin already points to the same object
		// gf: correction : the test is illegal since either iterator might not be dereferenceable
		//if (my_begin.operator->() != it.operator->()) {
			// assign it to my_begin and return my_begin
			my_begin.set_implementation_iterator(it);
		//}
		return my_begin;
	}
	;
	template<typename implementation_iterator_type> const iterator& get_end(
			const implementation_iterator_type& it) {
		// test if my_begin already points to the same object
		// gf: correction : the test is illegal since either iterator might not be dereferenceable
		//if (my_end.operator->() != it.operator->()) {
			// assign it to my_begin and return my_begin
			my_end.set_implementation_iterator(it);
		//}
		return my_end;
	}
	;
private:
	iterator my_begin;
	iterator my_end;
};

/** A utility function for converting iterators. */

/*
 template<typename object_type, typename implementation_const_iterator_type >
 polymorphic_forward_iterator<const object_type> get_const_iterator(
 const implementation_const_iterator_type& it) {
 boost::shared_ptr<forward_iterator<const object_type> > p=boost::shared_ptr<
 forward_iterator<const object_type> >(new forward_iterator_wrapper<const object_type,implementation_const_iterator_type>(it));
 return polymorphic_forward_iterator<object_type>(p);
 };
 */

/** This is an illustration on how to use the polymorphic_forward_iterator class. */
namespace simple_iterator_demo_class {

//class shape;
typedef int shape;
typedef boost::shared_ptr<shape> shape_ptr;

class shape_container;
class shape_container: public collection_const_base<shape_ptr> {
public:
	typedef boost::shared_ptr<shape_container> ptr;
	typedef boost::shared_ptr<const shape_container> const_ptr;
	virtual ~shape_container() {
	}
	;
	virtual void init() = 0;
};

class nonconst_shape_container: public virtual shape_container, public collection_nonconst_base<
		shape_ptr> {
public:
	typedef boost::shared_ptr<nonconst_shape_container> ptr;
	typedef boost::shared_ptr<const nonconst_shape_container> const_ptr;
	//	typedef collection_nonconst_base<shape_ptr>::iterator iterator;
	virtual ~nonconst_shape_container() {
	}
	;
};

class stl_container: public virtual shape_container {
public:
	typedef boost::shared_ptr<stl_container> ptr;
	typedef boost::shared_ptr<const stl_container> const_ptr;
	typedef shape_container::const_iterator const_iterator;
	virtual const const_iterator& begin() const {
		return get_const_begin(my_set.begin());
	}
	;
	virtual const const_iterator& end() const {
		return get_const_end(my_set.end());
	}
	;
	void init() {
		shape_ptr p;
		p = shape_ptr(new shape(1));
		my_set.push_back(p);
		p = shape_ptr(new shape(2));
		my_set.push_back(p);
		p = shape_ptr(new shape(3));
		my_set.push_back(p);
	}
	;
protected:
	std::list<shape_ptr> my_set;
};

class acc_stl_container: public stl_container, public collection_nonconst_base<shape_ptr> {
public:
	typedef boost::shared_ptr<acc_stl_container> ptr;
	typedef boost::shared_ptr<const acc_stl_container> const_ptr;

	typedef collection_nonconst_base<shape_ptr>::iterator iterator;

	/** The const versions of begin() and end() need to be defined inside
	 * the derived class. If they are inherited, the compiler will erroneously
	 * choose the nonconst version instead.
	 */
	virtual const const_iterator& begin() const {
		return stl_container::begin();
	}
	;
	virtual const const_iterator& end() const {
		return stl_container::end();
	}
	;

	virtual const iterator& begin() {
		return get_begin(this->my_set.begin());
	}
	;
	virtual const iterator& end() {
		return get_end(this->my_set.end());
	}
	;
};

}

}

#endif /* COLLECTION_BASE_H_ */
