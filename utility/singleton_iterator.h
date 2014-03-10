/*
 * singleton_iterator.h
 *
 *  Created on: Jun 22, 2009
 *      Author: frehse
 */

#ifndef SINGLETON_ITERATOR_H_
#define SINGLETON_ITERATOR_H_

/** An iterator class for a singleton_set.
 * An iterator into a singleton set has only two states:
 * True (at the beginning) or False (at the end).
 */

template<typename object_type> class singleton_const_iterator;

template<typename object_type> class singleton_iterator {
public:
	typedef object_type* pointer;
	typedef object_type& reference;
	typedef object_type value_type;
	typedef std::forward_iterator_tag iterator_category;
	typedef ptrdiff_t difference_type;

	singleton_iterator(object_type* obj) :
		my_object_ptr(obj), my_pos(true) {
	}
	;

	singleton_iterator(object_type* obj, bool pos) :
		my_object_ptr(obj), my_pos(pos) {
	}
	;

	virtual ~singleton_iterator() {
	}
	;

	virtual object_type* operator->() {
		if (my_pos)
			return my_object_ptr;
		else
			throw std::runtime_error("error dereferencing singleton_iterator at end()");
	}
	;
	virtual object_type& operator*() {
		if (my_pos)
			return *my_object_ptr;
		else
			throw std::runtime_error("error dereferencing singleton_iterator at end()");
	}
	;
	virtual const object_type* operator->() const {
		if (my_pos)
			return my_object_ptr;
		else
			throw std::runtime_error("error dereferencing singleton_iterator at end()");
	}
	;
	virtual const object_type& operator*() const {
		if (my_pos)
			return *my_object_ptr;
		else
			throw std::runtime_error("error dereferencing singleton_iterator at end()");
	}
	;
	virtual singleton_iterator<object_type>& operator++() {
		my_pos = false;
		return *this;
	}
	;

	virtual bool operator==(const singleton_iterator<object_type>& it) const {
		return my_pos == it.my_pos && (my_object_ptr == it.my_object_ptr);
	}
	;
	virtual bool operator!=(const singleton_iterator<object_type>& it) const {
		return !operator==(it);
	}
	;

	friend class singleton_const_iterator<object_type>;
private:
	object_type* my_object_ptr;
	bool my_pos;
};

template<typename object_type> class singleton_const_iterator {
public:
	typedef const object_type* pointer;
	typedef const object_type& reference;
	typedef object_type value_type;
	typedef std::forward_iterator_tag iterator_category;
	typedef ptrdiff_t difference_type;

	singleton_const_iterator(const object_type* obj) :
		my_object_ptr(obj), my_pos(true) {
	}
	;

		singleton_const_iterator(const object_type* obj, bool pos) :
		my_object_ptr(obj), my_pos(pos) {
	}
	;

	virtual ~singleton_const_iterator() {
	}
	;

	virtual const object_type* operator->() const {
		if (my_pos)
			return my_object_ptr;
		else
			throw std::runtime_error("error dereferencing singleton_iterator at end()");
	}
	;
	virtual const object_type& operator*() const {
		if (my_pos)
			return *my_object_ptr;
		else
			throw std::runtime_error("error dereferencing singleton_iterator at end()");
	}
	;
	virtual singleton_const_iterator<object_type>& operator++() {
		my_pos = false;
		return *this;
	}
	;

	virtual bool operator==(const singleton_const_iterator<object_type>& it) const {
		return my_pos == it.my_pos && (my_object_ptr == it.my_object_ptr);
	}
	;
	virtual bool operator!=(const singleton_const_iterator<object_type>& it) const {
		return !operator==(it);
	}
	;

	virtual bool operator==(const singleton_iterator<object_type>& it) const {
		return (my_object_ptr == it.my_object_ptr) && my_pos == it.my_pos;
	}
	;
	virtual bool operator!=(const singleton_iterator<object_type>& it) const {
		return !operator==(it);
	}
	;

private:
	const object_type* my_object_ptr;
	bool my_pos;
};

#endif /* SINGLETON_ITERATOR_H_ */
