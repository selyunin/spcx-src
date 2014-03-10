/*
 * proxy_forward_const_iterator.h
 *
 *  Created on: Jun 20, 2009
 *      Author: frehse
 */

#ifndef PROXY_FORWARD_CONST_ITERATOR_H_
#define PROXY_FORWARD_CONST_ITERATOR_H_

#include "utility/simple_iterators/proxy_forward_iterator.h"

namespace simple_iterators {

template<typename object_type> class proxy_forward_const_iterator: public forward_const_iterator_interface<
		object_type> {
public:
	typedef ptrdiff_t difference_type;

	proxy_forward_const_iterator(forward_const_iterator_interface<object_type>* it_p) :
		my_it(it_p) {
	}
	;
	/*	proxy_forward_iterator(const forward_iterator_interface<object_type>& it_p) {
	 my_it=it_p.clone();
	 }
	 ;
	 */
	proxy_forward_const_iterator() :
		my_it(0) {
	}
	;

	proxy_forward_const_iterator(const proxy_forward_const_iterator<object_type>& it) {
		if (this != &it) {
			// make a copy of the object pointed to by my_it
			//std::cout << "new const proxy iter " << this << " copied from " << &it << std::endl<< std::flush;
			my_it = it.const_clone();
			//std::cout << "csetting (new) ptr in const proxy iter " << this << " with ptr " << my_it<< std::endl << std::flush;
		}
	}
	;

	template<typename implementation_iterator_type> proxy_forward_const_iterator(
			const implementation_iterator_type& it);

	template<typename dereferencing_policy>
	static proxy_forward_const_iterator<object_type> //
	create(const typename dereferencing_policy::iterator_type& it);

	virtual ~proxy_forward_const_iterator() {
		/* std::cout << "deleting proxy const_iter " << this << " with ptr " << my_it << std::endl << std::flush; */
		delete my_it;
		my_it = 0;
	}
	;
	/* Virtual copy constructor. */
	virtual forward_const_iterator_interface<object_type>* const_clone() const {
		if (my_it)
		return my_it->const_clone();
		else
		return 0;
	}
	;
	virtual const object_type* operator->() const {
		if (my_it)
		return my_it->operator->();
		else
		return 0;
	}
	;
	virtual const object_type& operator*() const {
		assert(my_it);
		return my_it->operator*();
	}
	;
	/*
	 virtual forward_iterator_interface<object_type>& operator=(const forward_iterator_interface<
	 object_type>& it) {
	 // make a copy of the object pointed to by my_it
	 delete my_it;
	 my_it = it.clone();
	 return *this;
	 }
	 ;
	 */
	virtual proxy_forward_const_iterator<object_type>& operator=(
			const proxy_forward_const_iterator<object_type>& it) {
		if (this != &it) {
			// make a copy of the object pointed to by my_it
			delete my_it;
			my_it = it.const_clone();
		}
		return *this;
	}
	;
	virtual proxy_forward_const_iterator<object_type>& operator=(const
			proxy_forward_iterator<
			object_type>& it) {
		// make a copy of the object pointed to by my_it
		delete my_it;
		my_it = it.const_clone();
		return *this;
	}
	;
	virtual forward_const_iterator_interface<object_type>& operator++() {
		assert(my_it);
		my_it->operator++();
		assert(my_it);
		return *this;
	}
	;
	virtual bool operator==(const forward_iterator_interface<object_type>& it) const {
		//		if ((forward_iterator_interface<object_type>*) this == &it)
		//			return true;

		// gf: the following is invalid if either pointer is at end()
		//return this->operator->() == it.operator->();

		// assume they are both of the same type and directly test iterator
		assert(my_it);
		const proxy_forward_iterator<object_type>& jt=static_cast<const proxy_forward_iterator<object_type>& >(it);
		return my_it->operator==(*jt.my_it);
	}
	;
	virtual bool operator!=(const forward_iterator_interface<object_type>& it) const {
		return !operator==(it);
	}
	;
	virtual bool operator==(const forward_const_iterator_interface<object_type>& it) const {
		//		if ((forward_iterator_interface<object_type>*) this == &it)
		//			return true;
		// gf: the following is invalid if either pointer is at end()
		//return this->operator->() == it.operator->();

		// assume they are both of the same type and directly test iterator
		assert(my_it);
		const proxy_forward_const_iterator<object_type>& jt=static_cast<const proxy_forward_const_iterator<object_type>& >(it);
		return my_it->operator==(*jt.my_it);
	}
	;
	virtual bool operator!=(const forward_const_iterator_interface<object_type>& it) const {
		return !operator==(it);
	}
	;
	/** Attempt to get back at the original iterator by dynamic casting. */
	template<typename implementation_type> implementation_type&retrieve_implementation_iterator() {
		if (forward_iterator_wrapper<implementation_type> * jt
				= dynamic_cast<forward_iterator_wrapper<implementation_type>*> (my_it)) {
			return jt->get_implementation_iterator();
		} else {
			throw std::runtime_error("cannot retrieve implementation iterator");
			/* do something to make the compiler warnings go away */
			implementation_type* dummy = new implementation_type;
			return *dummy;
		}
	}
	;

	/** Set a new target implementation iterator. */
	template<typename implementation_type> void set_implementation_iterator(implementation_type it);

	template<typename object_typeT> friend class proxy_forward_iterator;
protected:
	forward_const_iterator_interface<object_type>* my_it;
};

	template<typename implementation_iterator_type> struct proxy_forward_const_iterator_copier {
		static forward_const_iterator_interface<typename implementation_iterator_type::value_type>* implement(
				const implementation_iterator_type& it) {
			//std::cout << "wrapping const iter " << &it << " obj " << it.operator->() << std::endl<< std::flush;
			forward_const_iterator_interface<typename implementation_iterator_type::value_type>* p =
			new forward_const_iterator_wrapper<implementation_iterator_type> (it);
			assert(p);
			//std::cout << " w--> " << p << std::endl << std::flush;
			return p;
		}
		;
	};

	///** Partially specialize copying for derived classes of forward_iterator_interface. */
	template<typename object_type> struct proxy_forward_const_iterator_copier<
	proxy_forward_iterator<object_type> > {
		static forward_const_iterator_interface<object_type>* implement(const proxy_forward_iterator<
				object_type>& it) {
			//std::cout << "cloning " << it.operator->() << std::endl << std::flush;
			forward_const_iterator_interface<object_type>* p = it.const_clone();
			//std::cout << " c--> " << p << std::endl << std::flush;
			return p;
		}
		;
	};

	template<typename object_type> struct proxy_forward_const_iterator_copier<
	proxy_forward_const_iterator<object_type> > {
		static forward_const_iterator_interface<object_type>* implement(
				const proxy_forward_const_iterator<object_type>& it) {
			//std::cout << "cloning " << it.operator->() << std::endl << std::flush;
			forward_const_iterator_interface<object_type>* p = it.const_clone();
			//std::cout << " c--> " << p << std::endl << std::flush;
			return p;
		}
		;
	};

	template<typename implementation_iterator_type, typename dereferencing_policy> struct proxy_forward_const_iterator_copier_der {
		static forward_const_iterator_interface<typename dereferencing_policy::object_type>* implement(
				const implementation_iterator_type& it) {
			//std::cout << "wrapping_d " << it.operator->() << std::endl << std::flush;
			forward_const_iterator_interface<typename dereferencing_policy::object_type> * p =
			new forward_const_iterator_wrapper<implementation_iterator_type,
			dereferencing_policy> (it);
			assert(p);
			//std::cout << " w--> " << p << std::endl << std::flush;
			return p;
		}
		;
	};

	/** Partially specialize copying for derived classes of forward_iterator_interface. */
	template<typename dereferencing_policy> struct proxy_forward_const_iterator_copier_der<
	proxy_forward_iterator<typename dereferencing_policy::object_type> , dereferencing_policy> {
		static forward_const_iterator_interface<typename dereferencing_policy::object_type>* implement(
				const proxy_forward_iterator<typename dereferencing_policy::object_type>& it) {
			//std::cout << "cloning_d " << it.operator->() << std::endl << std::flush;
			forward_const_iterator_interface<typename dereferencing_policy::object_type>* p =
			it.const_clone();
			//std::cout << " c--> " << p << std::endl << std::flush;
			return p;
		}
		;
	};

	/** Partially specialize copying for derived classes of forward_iterator_interface. */
	template<typename dereferencing_policy> struct proxy_forward_const_iterator_copier_der<
	proxy_forward_const_iterator<typename dereferencing_policy::object_type> ,
	dereferencing_policy> {
		static forward_const_iterator_interface<typename dereferencing_policy::object_type>* implement(
				const proxy_forward_const_iterator<typename dereferencing_policy::object_type>& it) {
			//std::cout << "cloning_d " << it.operator->() << std::endl << std::flush;
			forward_const_iterator_interface<typename dereferencing_policy::object_type>* p =
			it.const_clone();
			//std::cout << " c--> " << p << std::endl << std::flush;
			return p;
		}
		;
	};

	template<typename object_type>
	template<typename implementation_iterator_type> proxy_forward_const_iterator<object_type>::proxy_forward_const_iterator(
			const implementation_iterator_type& it) {
		//std::cout << "new proxy_const_iter start" << this << std::endl << std::flush;
		my_it
		= proxy_forward_const_iterator_copier<implementation_iterator_type>::implement(
				it);
		//std::cout << "new proxy_const_iter " << this << std::endl << std::flush;
	}
	;

	template<typename object_type>
	template<typename dereferencing_policy> proxy_forward_const_iterator<object_type> proxy_forward_const_iterator<
	object_type>::create(const typename dereferencing_policy::iterator_type& it) {
		//std::cout << "create proxy_const_iter start"<< std::endl << std::flush;
		proxy_forward_const_iterator<object_type> jt;
		jt.my_it=(proxy_forward_const_iterator_copier_der<
				typename dereferencing_policy::iterator_type, dereferencing_policy>::implement( it));
		//std::cout << "created proxy_const_iter " << &jt << std::endl << std::flush;
		return jt;
	}
	;

	template<typename object_type>
	template<typename implementation_type> void proxy_forward_const_iterator<object_type>::set_implementation_iterator(
			implementation_type it) {
		if (my_it) {
			if (forward_const_iterator_wrapper<implementation_type> * jt
					= dynamic_cast<forward_const_iterator_wrapper<implementation_type>*> (my_it)) {
				jt->get_implementation_iterator() = it;
			} else {
				// must create a new one
				delete my_it;
				my_it = proxy_forward_const_iterator_copier<implementation_type>::implement(it);
			}
		} else {
			// must create a new one
			my_it = proxy_forward_const_iterator_copier<implementation_type>::implement(it);
		}
	}
	;

}

#endif /* PROXY_FORWARD_CONST_ITERATOR_H_ */
