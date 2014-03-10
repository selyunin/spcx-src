/*
 * forward_iterator_wrapper.h
 *
 *  Created on: Jun 20, 2009
 *      Author: frehse
 */

#ifndef FORWARD_ITERATOR_WRAPPER_H_
#define FORWARD_ITERATOR_WRAPPER_H_

#include "utility/simple_iterators/forward_iterator_interface.h"

namespace simple_iterators {

/** This wrapper hides the implementation type of the actual iterator from the users, which
 * have access to the iterator via the forward_iterator interface.
 * It may be used to turn any STL-like iterator into a forward_iterator.
 * */

/** This structure defines how the iterator is dereferenced. Normal use is
 * it-> and *it, but alternative implementations can insert additional code,
 * such as accessing only the first item of a pair etc. */
template<typename arg_iterator_type>
struct normal_dereferencing {
	typedef arg_iterator_type iterator_type;
typedef	typename iterator_type::value_type object_type;
	static object_type& get_reference(
			iterator_type& my_it) {
		return my_it.operator*();
	}
	;
	static object_type* get_pointer(
			iterator_type& my_it) {
		return my_it.operator->();
	}
	;
	static const object_type& get_const_reference(const
			iterator_type& my_it) {
		return my_it.operator*();
	}
	;
	static const object_type* get_const_pointer(const
			iterator_type& my_it) {
		return my_it.operator->();
	}
	;
};

/** Dereference the first out of a pair. */
template<typename arg_iterator_type,typename arg_obj_type=typename arg_iterator_type::value_type::first_type>
struct first_in_pair_dereferencing {
	typedef arg_iterator_type iterator_type;
	typedef arg_obj_type object_type;
	static object_type& get_reference(
			iterator_type& my_it) {
		return (my_it.operator*()).first;
	}
	;
	static object_type* get_pointer(
			iterator_type& my_it) {
		return &((my_it.operator*()).first);
	}
	;
	static const object_type& get_const_reference(const
			iterator_type& my_it) {
		return (my_it.operator*()).first;
	}
	;
	static const object_type* get_const_pointer(const
			iterator_type& my_it) {
		return &((my_it.operator*()).first);
	}
	;
};

/** Dereference the second out of a pair. */
template<typename arg_iterator_type>
struct second_in_pair_dereferencing {
	typedef arg_iterator_type iterator_type;
	typedef typename iterator_type::value_type::second_type object_type;
	static object_type& get_reference(
			iterator_type& my_it) {
		return (my_it.operator*()).second;
	}
	;
	static object_type* get_pointer(
			iterator_type& my_it) {
		return &((my_it.operator*()).second);
	}
	;
	static const object_type& get_const_reference(const
			iterator_type& my_it) {
		return (my_it.operator*()).second;
	}
	;
	static const object_type* get_const_pointer(const
			iterator_type& my_it) {
		return &((my_it.operator*()).second);
	}
	;
};

template<typename iterator_type, typename dereferencing_policy> class forward_const_iterator_wrapper;

template<typename iterator_type,
typename dereferencing_policy = normal_dereferencing<iterator_type> > class forward_iterator_wrapper: public forward_iterator_interface<
typename dereferencing_policy::object_type> {
public:
	typedef forward_iterator_wrapper<iterator_type,dereferencing_policy> mytype;
	typedef typename iterator_type::difference_type difference_type;
	typedef typename dereferencing_policy::object_type object_type;

	explicit forward_iterator_wrapper(iterator_type it) :
	my_it(it) {
		//object_type* test=operator->(); // for testing if the types are ok
	}
	;
	virtual ~forward_iterator_wrapper() {
	}
	;
	/* Virtual copy constructor. */
	virtual forward_iterator_interface<object_type>* clone() const {
		return new forward_iterator_wrapper<iterator_type,dereferencing_policy>(my_it);
	}
	;
	virtual forward_const_iterator_interface<object_type>* const_clone() const {
		return new forward_const_iterator_wrapper<iterator_type,dereferencing_policy>(my_it);
	}
	;
	virtual object_type* operator->() {
		return dereferencing_policy::get_pointer(my_it);
	}
	;
	virtual object_type& operator*() {
		return dereferencing_policy::get_reference(my_it);
	}
	;
	virtual const object_type* operator->() const {
		return dereferencing_policy::get_const_pointer(my_it);
	}
	;
	virtual const object_type& operator*() const {
		return dereferencing_policy::get_const_reference(my_it);
	}
	;
	virtual forward_iterator_wrapper<iterator_type,dereferencing_policy>& operator=(const forward_iterator_wrapper<iterator_type,dereferencing_policy>& it) {
		if (this!=&it) {
			my_it=it.get_implementation_iterator();
		}
		return *this;
	}
	;
	virtual forward_iterator_interface<object_type>& operator++() {
		++my_it;
		return *this;
	}
	;
	virtual bool operator==(const forward_iterator_interface<object_type>& it) const {
		if (this==&it)
		return true;
		const mytype& jt=static_cast<const mytype&>(it);
		return (my_it==jt.my_it);
		//return this->operator->()==it.operator->();
	}
	;
	virtual bool operator!=(const forward_iterator_interface<object_type>& it) const {
		return !operator==(it);
	}
	;
	virtual bool operator==(const forward_const_iterator_interface<object_type>& it) const {
		const forward_const_iterator_wrapper<iterator_type,dereferencing_policy>& jt=static_cast<const forward_const_iterator_wrapper<iterator_type,dereferencing_policy>&>(it);
		return (my_it==jt.my_it);
		//return this->operator->()==it.operator->();
	}
	;
	virtual bool operator!=(const forward_const_iterator_interface<object_type>& it) const {
		return !operator==(it);
	}
	;
	virtual iterator_type& get_implementation_iterator() {
		return my_it;
	}
	;
	virtual const iterator_type& get_implementation_iterator() const {
		return my_it;
	}
	;
	template<typename iterator_typeT,typename dereferencing_policyT> friend class forward_const_iterator_wrapper;
protected:
	iterator_type my_it;
private:
	forward_iterator_wrapper();
};

template<typename iterator_type,
typename dereferencing_policy = normal_dereferencing<iterator_type> > class forward_const_iterator_wrapper: public forward_const_iterator_interface<
typename dereferencing_policy::object_type> {
public:
	typedef forward_const_iterator_wrapper<iterator_type,dereferencing_policy> mytype;
	typedef typename iterator_type::difference_type difference_type;
	typedef typename dereferencing_policy::object_type object_type;

	explicit forward_const_iterator_wrapper(iterator_type it) :
	my_it(it) { //assert(typename iterator_type::value_type* tst(new object_type()));
	}
	;
	virtual ~forward_const_iterator_wrapper() {
	}
	;
	/* Virtual copy constructor. */
	virtual forward_const_iterator_interface<object_type>* const_clone() const {
		return new forward_const_iterator_wrapper<iterator_type,dereferencing_policy>(my_it);
	}
	;
	virtual const object_type* operator->() const {
		return dereferencing_policy::get_const_pointer(my_it);
	}
	;
	virtual const object_type& operator*() const {
		return dereferencing_policy::get_const_reference(my_it);
	}
	;
	virtual forward_const_iterator_wrapper<iterator_type,dereferencing_policy>& operator=(const forward_iterator_wrapper<iterator_type,dereferencing_policy>& it) {
		my_it=it.get_implementation_iterator();
		return *this;
	}
	;
	virtual forward_const_iterator_wrapper<iterator_type,dereferencing_policy>& operator=(const forward_const_iterator_wrapper<iterator_type,dereferencing_policy>& it) {
		if (this!=&it) {
			my_it=it.get_implementation_iterator();
		}
		return *this;
	}
	;
	virtual forward_const_iterator_interface<object_type>& operator++() {
		++my_it;
		return *this;
	}
	;
	virtual bool operator==(const forward_iterator_interface<object_type>& it) const {
		const forward_iterator_wrapper<iterator_type,dereferencing_policy>& jt=static_cast<const forward_iterator_wrapper<iterator_type,dereferencing_policy>&>(it);
		return (my_it==jt.my_it);
		// return this->operator->()==it.operator->();
	}
	;
	virtual bool operator!=(const forward_iterator_interface<object_type>& it) const {
		return !operator==(it);
	}
	;
	virtual bool operator==(const forward_const_iterator_interface<object_type>& it) const {
		if (this==&it)
		return true;
		const mytype& jt=static_cast<const mytype&>(it);
		return (my_it==jt.my_it);
		// return this->operator->()==it.operator->();
	}
	;
	virtual bool operator!=(const forward_const_iterator_interface<object_type>& it) const {
		return !operator==(it);
	}
	;
	virtual iterator_type& get_implementation_iterator() {
		return my_it;
	}
	;
	virtual const iterator_type& get_implementation_iterator() const {
		return my_it;
	}
	;
	template<typename iterator_typeT,typename dereferencing_policyT> friend class forward_iterator_wrapper;
protected:
	iterator_type my_it;
private:
	forward_const_iterator_wrapper();
};

}

#endif /* FORWARD_ITERATOR_WRAPPER_H_ */
