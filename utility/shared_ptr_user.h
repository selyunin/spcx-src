#ifndef SHARED_PTR_USER_H_
#define SHARED_PTR_USER_H_

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

/** An interface class to provide useful functions with boost::shared_ptr.
 * The instantiation takes the name of the derived class:
 * \code class foo : public shared_ptr_user<foo> {...} \endcode
 *
 */
template <class T> class shared_ptr_user :
	public boost::enable_shared_from_this<T> {
public:
	shared_ptr_user() {
	}
	;
	virtual ~shared_ptr_user() {
	}
	;

	/** Some typedefs to make the code prettier. Users are highly encouraged to use thes instead of boost::shared_ptr<foo>. */
	typedef boost::shared_ptr<T> ptr;
	typedef boost::shared_ptr<const T> const_ptr;

	/** Return a shared_ptr to *this.
	 * It is necessary because \code return shared_ptr<this>; \endcode is not allowed.
	 * \see http://www.boost.org/doc/libs/1_35_0/libs/smart_ptr/sp_techniques.html#from_this
	 */
	ptr myself() {
		return boost::enable_shared_from_this<T>::shared_from_this();
	}

	/** Return a shared_ptr to const *this.
	 */
	const_ptr const_myself() const {
		return boost::enable_shared_from_this<T>::shared_from_this();
	}
};

/** An interface class to provide useful functions with boost::shared_ptr.
 * The instantiation takes the name of the derived class:
 * \code class foo : public shared_ptr_user<foo> {...} \endcode
 *
 */
template <class T> class ptr_interface :
	public boost::enable_shared_from_this<T> {
public:
	virtual ~ptr_interface() {
	}
	;

	/** Some typedefs to make the code prettier. Users are highly encouraged to use thes instead of boost::shared_ptr<foo>. */
	typedef boost::shared_ptr<T> ptr;
	typedef boost::shared_ptr<const T> const_ptr;

	/** Return a shared_ptr to *this.
	 * It is necessary because \code return shared_ptr<this>; \endcode is not allowed.
	 * \see http://www.boost.org/doc/libs/1_35_0/libs/smart_ptr/sp_techniques.html#from_this
	 */
	ptr get_ptr() {
		ptr p=boost::enable_shared_from_this<T>::shared_from_this();
		return p;
		//return boost::enable_shared_from_this<T>::shared_from_this();
	};

	/** Return a shared_ptr to const *this.
	 */
	const_ptr get_const_ptr() const {
		const_ptr p=boost::enable_shared_from_this<T>::shared_from_this();
		return p;
		//return boost::enable_shared_from_this<T>::shared_from_this();
	};

	/** Return a shared_ptr to b.
	 */
	template <class target_class>
	boost::shared_ptr<target_class> cast_ptr() {
		typename boost::shared_ptr<target_class> p=boost::static_pointer_cast<target_class, T>(get_ptr());
		return p;
	};

	/** Return a shared_ptr to b.
	 */
	template <class base_class>
	static ptr cast_ptr(const boost::shared_ptr<base_class>& b) {
		ptr p=boost::static_pointer_cast<T, base_class>(b);
		return p;
	};

	/** Return a const shared_ptr to b.
	 */
	template <class base_class>
	static const_ptr cast_ptr(const boost::shared_ptr<const base_class>& b) {
		return boost::static_pointer_cast<const T, const base_class>(b);
	};

	/** Return true if *this can be cast to an object of class T (test is done using const).
	 */
	template <class base_class>
	static bool is_ptr(const boost::shared_ptr<base_class>& b) {
		return dynamic_cast<const T>(b.get());
	};
};

/** An interface class to provide useful functions with boost::shared_ptr,
 * aimed at a derived class that wishes to return pointers to a base class.
 * The instantiation takes the name of the derived class:
 * \code class foo : public shared_ptr_user<foo> {...} \endcode
 *
 */
template <class T, class base_class> class DEPRECATED_derived_ptr_interface : public
		ptr_interface<T> {
private:
	/** Some typedefs to make the code prettier. Users are highly encouraged to use thes instead of boost::shared_ptr<foo>. */
	typedef boost::shared_ptr<base_class> base_ptr;
	typedef boost::shared_ptr<const base_class> const_base_ptr;

public:
	virtual ~DEPRECATED_derived_ptr_interface() {
	}
	;

	/** Return a shared_ptr to *this.
	 * It is necessary because \code return shared_ptr<this>; \endcode is not allowed.
	 * \see http://www.boost.org/doc/libs/1_35_0/libs/smart_ptr/sp_techniques.html#from_this
	 */
//	base_ptr get_base_ptr() {
//		return base_class::get_ptr();
//	};

	/** Return a shared_ptr to const *this.
	 */
//	const_base_ptr get_const_base_ptr() const {
//		return base_class::get_const_ptr();
//	};

	/** Return a shared_ptr to *this.
	 */
	static typename ptr_interface<T>::ptr get_ptr(const base_ptr& b) {
		return boost::static_pointer_cast<T, base_class>(b);
	};

	/** Return a shared_ptr to *this.
	 */
	static typename ptr_interface<T>::const_ptr get_const_ptr(const base_ptr& b) {
		return boost::static_pointer_cast<const T, base_class>(b);
	};

	/** Return true if *this can be cast to an object of class T (test is done using const).
	 */
	static bool is_ptr(const base_ptr& b) {
		return dynamic_cast<const T>(b.get());
	};
};

#endif /*SHARED_PTR_USER_H_*/

