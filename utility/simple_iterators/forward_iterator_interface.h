/*
 * forward_iterator_interface.h
 *
 *  Created on: Jun 20, 2009
 *      Author: frehse
 */

#ifndef FORWARD_ITERATOR_INTERFACE_H_
#define FORWARD_ITERATOR_INTERFACE_H_

#include <iterator>
#include <iostream>
#include <stdexcept>
//#include <boost/shared_ptr.hpp>

namespace simple_iterators {

/** An abstract base class, defining the interface for a forward iterator (an "input iterator" in STL speak).
 * The iterator can be incremented, and dereferences to a const reference to an object of type \p object_type.
 */

template<typename object_type> class forward_const_iterator_interface;

template<typename object_type> class forward_iterator_interface {
public:
	typedef object_type* pointer;
	typedef object_type& reference;
	typedef object_type value_type;
	typedef std::forward_iterator_tag iterator_category;

	virtual ~forward_iterator_interface() {
	}
	;
	/* Virtual copy constructor. */
	virtual forward_iterator_interface<object_type>* clone() const = 0;
	virtual forward_const_iterator_interface<object_type>* const_clone() const = 0;
	virtual object_type* operator->() = 0;
	virtual object_type& operator*() = 0;
	virtual const object_type* operator->() const = 0;
	virtual const object_type& operator*() const = 0;
	virtual forward_iterator_interface<object_type>& operator++() = 0; // incrementing
	virtual bool operator==(const forward_iterator_interface<object_type>& it) const = 0;
	virtual bool operator!=(const forward_iterator_interface<object_type>& it) const = 0;
	virtual bool operator==(const forward_const_iterator_interface<object_type>& it) const = 0;
	virtual bool operator!=(const forward_const_iterator_interface<object_type>& it) const = 0;
};

template<typename object_type> class forward_const_iterator_interface {
public:
	typedef const object_type* pointer;
	typedef const object_type& reference;
	typedef object_type value_type;
	typedef std::forward_iterator_tag iterator_category;

	virtual ~forward_const_iterator_interface() {
	}
	;
	/* Virtual copy constructor. */
	virtual forward_const_iterator_interface<object_type>* const_clone() const = 0;
	virtual pointer operator->() const = 0;
	virtual reference operator*() const = 0;
	virtual forward_const_iterator_interface<object_type>& operator++() = 0; // incrementing
	virtual bool operator==(const forward_iterator_interface<object_type>& it) const = 0;
	virtual bool operator!=(const forward_iterator_interface<object_type>& it) const = 0;
	virtual bool operator==(const forward_const_iterator_interface<object_type>& it) const = 0;
	virtual bool operator!=(const forward_const_iterator_interface<object_type>& it) const = 0;
};

}

#endif /* FORWARD_ITERATOR_INTERFACE_H_ */
