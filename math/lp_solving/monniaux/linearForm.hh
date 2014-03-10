#ifndef _LINEAR_FORM_HH
#define _LINEAR_FORM_HH

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <assert.h>

/**
 LINEAR FORMS
 **/
#define NOTHROW throw()

#ifdef BOOST_VECTORS
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

//#include <ext/bitmap_allocator.h>
//#include <ext/mt_allocator.h>

using namespace boost::numeric::ublas;
#define LINEAR_FORM_INLINE inline
#else
#define LINEAR_FORM_INLINE
#endif

template <class coefficient> class LinearForm {
public:
#ifdef BOOST_VECTORS
	typedef compressed_vector<coefficient> V;
#else
	typedef std::vector<coefficient> V;
#endif

private:
	V coefficients;
	const std::vector<std::string>* variableNames;

public:
#ifdef BOOST_VECTORS
	LinearForm(const std::vector<std::string>& variableNames0) :
	coefficients(variableNames0.size()), variableNames(&variableNames0) {
	}

	LinearForm(const std::vector<std::string>& variableNames0, const V& coeffs) :
	coefficients(coeffs), variableNames(&variableNames0) {
	}

	void erase(unsigned i) {
		coefficients.erase_element(i);
	}
#else
	void erase(unsigned i) {
		coefficients[i] = 0;
	}
#endif

	V& getCoefficients() {
		return coefficients;
	}

	const V& getCoefficients() const {
		return coefficients;
	}

	typedef typename V::const_iterator const_iterator;
	typedef typename V::iterator iterator;

	const_iterator begin() const {
		return coefficients.begin();
	}
	const_iterator end() const {
		return coefficients.end();
	}
	const_iterator find(unsigned i) const {
		return coefficients.find(i);
	}
	iterator begin() {
		return coefficients.begin();
	}
	iterator end() {
		return coefficients.end();
	}
	iterator find(unsigned i) {
		return coefficients.find(i);
	}
	bool hasNonZero(unsigned i) const {
		const_iterator it=find(i);
		// gff 
		// return (it!=end() && it.index()==i && !(*it).isZero());
		return (it!=end() && it.index()==i && (*it)!=coefficient(0));
	}

	bool isZero() const {
		for (const_iterator it=begin(); it!=end(); it++) {
			if (*it != 0)
				return false;
		}
		return true;
	}

	void clear() {
		//coefficients.clear();
	}

	const std::vector<std::string>& getVariableNames() const {
		return *variableNames;
	}

	unsigned dimension() const {
		return coefficients.size();
	}

	unsigned filled() const {
		return coefficients.filled();
	}

	LinearForm<coefficient>
			& operator=(const LinearForm<coefficient>& from) NOTHROW;

#ifdef BOOST_VECTORS
	template <class from> LinearForm(const LinearForm<from>& f) :
	coefficients(f.getVariableNames().size()),
	variableNames(& f.getVariableNames()) {
		for(typename LinearForm<from>::const_iterator it=f.begin();
				it != f.end(); it++) {
			coefficients[it.index()] = numericConversion(*it);
		}
	}
#else
	template <class from> LinearForm(const LinearForm<from>& f) :
		variableNames(&f.getVariableNames()) {
		for (unsigned k=0; k<variableNames->size(); k++) {
			coefficients.push_back(numericConversion(f[k]));
		}
	}
#endif

	unsigned size() const {
		return coefficients.size();
	}

#ifdef BOOST_VECTORS
	void resize(unsigned x) {
		coefficients.resize(x);
	}

	typename V::const_reference
#else
	const coefficient&
#endif
	operator[](unsigned i) const {
		assert(i < size());
		return coefficients[i];
	}

#ifdef BOOST_VECTORS
	typename V::reference
#else
	coefficient&
#endif
	operator[](unsigned i) {
		assert(i < size());
		return coefficients[i];
	}

	const std::string& variableName(unsigned i) const {
		assert(i < size());
		return (*variableNames)[i];
	}

	void addMonomial(unsigned i, coefficient c) NOTHROW {
		assert(i < size());
		coefficients[i] += c;
	}

	void addTimes(const coefficient& k, const LinearForm<coefficient>& u) NOTHROW;
	void
			assignTimes(const coefficient& k, const LinearForm<coefficient>& u) NOTHROW;

	bool operator==(const LinearForm<coefficient>& x) const NOTHROW;

	LinearForm<coefficient>& getForm() {
		return *this;
	}
	const LinearForm<coefficient>& getForm() const {
		return *this;
	}
};

template <class coefficient> bool LinearForm<coefficient>::operator==(const LinearForm<coefficient>& x) const NOTHROW {
	const_iterator aIt = begin(), bIt = x.begin();
	while (true) {
		while (aIt != end() && *aIt == coefficient(0)) aIt++;
		while (bIt != x.end() && *bIt == coefficient(0)) bIt++;

		if (aIt == end() || bIt == x.end()) break;

		if (aIt.index() != bIt.index()) return false;
		if (*aIt != *bIt) return false;
		aIt++; bIt++;
	}
	return (aIt == end() && bIt == x.end());
}

template <class coefficient> LINEAR_FORM_INLINE void LinearForm<coefficient>::addTimes(
		const coefficient& k, const LinearForm<coefficient>& u) NOTHROW {
#ifdef BOOST_VECTORS
#ifdef BOOST_VECTOR_OPS
	noalias(coefficients) += k*u.coefficients;
#else
	if (size() < u.size()) {
		resize(u.size());
	}
	for(typename LinearForm<coefficient>::const_iterator it=u.begin();
			it != u.end(); it++) {
		coefficients[it.index()] += k*(*it);
	}
#endif
#else
	for (unsigned i=0; i<size(); i++) {
		coefficients[i] += k*u.coefficients[i];
	}
#endif
}

template <class coefficient> LINEAR_FORM_INLINE void LinearForm<coefficient>::assignTimes(
		const coefficient& k, const LinearForm<coefficient>& u) NOTHROW {
#ifdef BOOST_VECTORS
#ifdef BOOST_VECTOR_OPS
	noalias(coefficients) = k*u.coefficients;
#else
	coefficients.clear();
	for(typename LinearForm<coefficient>::const_iterator it=u.begin();
			it != u.end(); it++) {
		coefficients[it.index()] = k*(*it);
	}
#endif
#else
	for (unsigned i=0; i<size(); i++) {
		coefficients[i] = k*u.coefficients[i];
	}
#endif
}

template <class coefficient> std::ostream& operator<<(std::ostream& out,
		const LinearForm<coefficient>& form) {
	out << "(+ ";
#ifdef BOOST_VECTORS
	for(typename LinearForm<coefficient>::const_iterator xi=form.begin();
			xi!=form.end(); xi++) {
		const coefficient& coef = *xi;
		unsigned i = xi.index();
#else
	for (unsigned i=0; i<form.size(); i++) {
		const coefficient& coef = form[i];
#endif
		if (coef != 0) {
			out << "(* " << coef << " " << form.variableName(i) << ")";
		}
	}
	out << ")";
	return out;
}

namespace std {
template<class coefficient> struct less<LinearForm<coefficient> > {
	bool operator()(const LinearForm<coefficient> & x, const LinearForm<coefficient> & y) NOTHROW {
		assert(x.size() == y.size());

#ifdef BOOST_VECTORS
		typename LinearForm<coefficient>::const_iterator xIt = x.begin(),
		yIt = y.begin();
		while(xIt != x.end() && yIt != y.end()) {
			if (xIt.index() < yIt.index()) return true;
			if (xIt.index()> yIt.index()) return false;
			if (*xIt < *yIt) return true;
			if (*xIt> *yIt) return false;
			xIt++;
			yIt++;
		}
		if (yIt == y.end()) return false;
		return true;

#else
		for (unsigned i=0; i<x.size(); i++) {
			if (coefficient(x[i]) < coefficient(y[i]))
				return true;
			if (coefficient(x[i]) > coefficient(y[i]))
				return false;
		}
		return false;
#endif
	}
};
}

template <class coefficient> LinearForm<coefficient>& LinearForm<coefficient>::operator=(const LinearForm<coefficient>& from) NOTHROW {
	if (size() < from.size()) {
		std::cerr << "form assignment error: " << size() << "!=" << from.size()
		<< std::endl;
		abort();
	}
	//coefficients.clear();

#ifdef BOOST_VECTORS
	for(const_iterator it=from.begin(); it!=from.end(); it++) {
		coefficients[it.index()] = *it;
	}
#else
	for(unsigned i=0; i<from.size(); i++) {
		coefficients[i] = from.coefficients[i];
	}
#endif
	return *this;
}

template <class coefficient> LINEAR_FORM_INLINE
LinearForm<coefficient> operator+(
		const LinearForm<coefficient> &a, const LinearForm<coefficient> &b) NOTHROW {
	assert(&a.getVariableNames() == &b.getVariableNames());

#ifdef BOOST_VECTORS
	return LinearForm<coefficient>(a.getVariableNames(),
			a.getCoefficients() + b.getCoefficients());
#else
	LinearForm<coefficient> r(a.getVariableNames());
	for (unsigned i=0; i<a.size(); i++) {
		r[i] = a[i]+b[i];
	}
	return r;
#endif
}

template <class coefficient> LINEAR_FORM_INLINE
LinearForm<coefficient> operator-(
		const LinearForm<coefficient> &a, const LinearForm<coefficient> &b) NOTHROW {
	assert(&a.getVariableNames() == &b.getVariableNames());

#ifdef BOOST_VECTORS
	return LinearForm<coefficient>(a.getVariableNames(),
			a.getCoefficients() - b.getCoefficients());
#else
	LinearForm<coefficient> r(a.getVariableNames());
	for (unsigned i=0; i<a.size(); i++) {
		r[i] = a[i]-b[i];
	}
	return r;
#endif
}

template <class coefficient, class scalar> LINEAR_FORM_INLINE
LinearForm<coefficient> operator*(
		const LinearForm<coefficient> &a, const scalar& k) NOTHROW {

#ifdef BOOST_VECTORS
#ifdef BOOST_VECTOR_OPS
	return LinearForm<coefficient>(a.getVariableNames(),
			a.getCoefficients() * k);
#else
	LinearForm<coefficient> r(a.getVariableNames());
	for(typename LinearForm<coefficient>::const_iterator it=a.begin();
			it!=a.end(); it++) {
		r[it.index()] = (*it)*k;
	}
	return r;
#endif
#else
	LinearForm<coefficient> r(a.getVariableNames());
	for (unsigned i=0; i<a.size(); i++) {
		r[i] = a[i]*k;
	}
	return r;
#endif
}

template <class coefficient, class scalar> LINEAR_FORM_INLINE
LinearForm<coefficient> operator/(
		const LinearForm<coefficient> &a, const scalar& k) NOTHROW {
#ifdef BOOST_VECTORS
#ifdef BOOST_VECTOR_OPS
	return LinearForm<coefficient>(a.getVariableNames(),
			a.getCoefficients() / k);
#else
	LinearForm<coefficient> r(a.getVariableNames());
	for(typename LinearForm<coefficient>::const_iterator it=a.begin();
			it!=a.end(); it++) {
		r[it.index()] = (*it)/k;
	}
	return r;
#endif
#else
	LinearForm<coefficient> r(a.getVariableNames());
	for (unsigned i=0; i<a.size(); i++) {
		r[i] = a[i]/k;
	}
	return r;
#endif
}

template <class coefficient> LINEAR_FORM_INLINE
LinearForm<coefficient> operator*(
		const coefficient& k, const LinearForm<coefficient> &a) NOTHROW {
	return a*k;
}

template <class coefficient> LINEAR_FORM_INLINE
LinearForm<coefficient>& operator+=(
		LinearForm<coefficient>& a, const LinearForm<coefficient>& b) NOTHROW {
	assert(&a.getVariableNames() == &b.getVariableNames());
#ifdef BOOST_VECTORS
#ifdef BOOST_VECTOR_OPS
	noalias(a.getCoefficients()) += b.getCoefficients();
	return a;
#else
	if (a.size() < b.size()) {
		a.resize(b.size());
	}
	for(typename LinearForm<coefficient>::const_iterator it=b.begin();
			it != b.end(); it++) {
		a[it.index()] += *it;
	}
	return a;
#endif
#else
	for (unsigned i=0; i<a.size(); i++) {
		a[i] += b[i];
	}
	return a;
#endif
}

template <class coefficient> LINEAR_FORM_INLINE
LinearForm<coefficient>& operator-=(
		LinearForm<coefficient>& a, const LinearForm<coefficient>& b) NOTHROW {
	assert(&a.getVariableNames() == &b.getVariableNames());
#ifdef BOOST_VECTOR_OPS
	noalias(a.getCoefficients()) -= b.getCoefficients();
	return a;
#else
	for (unsigned i=0; i<a.size(); i++) {
		a[i] -= b[i];
	}
	return a;
#endif
}

template <class coefficient> LinearForm<coefficient>& operator*=(
		LinearForm<coefficient>& a, const coefficient& k) {
#ifdef BOOST_VECTOR_OPS
	return (a.getCoefficients() *= k);
#else
	for (unsigned i=0; i<a.size(); i++) {
		a[i] *= k;
	}
	return a;
#endif
}

template <class coefficient> LINEAR_FORM_INLINE
LinearForm<coefficient> operator-(
		LinearForm<coefficient>& a) NOTHROW {
#ifdef BOOST_VECTORS
	return LinearForm<coefficient>(a.getVariableNames(), -a.getCoefficients());
#else
	LinearForm<coefficient> r(a.getVariableNames());
	for (unsigned i=0; i<a.size(); i++) {
		r[i] = -coefficient(a[i]);
	}
	return r;
#endif
}

#ifdef HAS_HASH_TEMPLATE
namespace std {
	template<class coefficient> class hash<LinearForm<coefficient> > {
		static inline size_t djb2(size_t a, size_t b) {
			return (a << 5) + a + b;
		}

	public:
		size_t operator()(const LinearForm<coefficient>& form) const NOTHROW {
			size_t ret=1;
			for(typename LinearForm<coefficient>::const_iterator it=form.begin();
					it != form.end(); it++) {
				ret = djb2(djb2(ret, it.index()), std::hash<coefficient>()(*it));
			}
			return ret;
		}
	};
}

template<class coefficient> class HashedLinearForm : LinearForm<coefficient> {
	size_t hashV;

public:
	const LinearForm<coefficient>& getForm() const {
		return *(static_cast<const LinearForm<coefficient>*>(this));
	}

	size_t hash() const {
		return hashV;
	}

	HashedLinearForm(const LinearForm<coefficient>& form) :
	LinearForm<coefficient>(form),
	hashV(std::hash<LinearForm<coefficient> >() (form)) {
	}

	HashedLinearForm(const HashedLinearForm<coefficient>& form) :
	LinearForm<coefficient>(form),
	hashV(form.hashV) {
	}
	typedef typename LinearForm<coefficient>::const_iterator const_iterator;

	const_iterator begin() const {
		return LinearForm<coefficient>::begin();
	}

	const_iterator end() const {
		return LinearForm<coefficient>::end();
	}

	const coefficient& operator[](unsigned i) const {
		return LinearForm<coefficient>::operator[](i);
	}

	bool operator==(const HashedLinearForm& x) const {
		return getForm() == x.getForm();
	}

	operator LinearForm<coefficient>() {
		return *this;
	}
};

namespace std {
	template<class coefficient> struct less<HashedLinearForm<coefficient> > {
		bool operator()(const HashedLinearForm<coefficient>& a,
				const HashedLinearForm<coefficient>& b) const {
			return std::less<LinearForm<coefficient> >() (a.getForm(), b.getForm());
		}
	};

	template<class coefficient> struct hash<HashedLinearForm<coefficient> > {
		size_t operator()(const HashedLinearForm<coefficient>& a) const {
			return a.hash();
		}
	};
}
#endif // HAS_HASH_TEMPLATE

template <class coefficient, class affine=coefficient> class AffineForm :
		LinearForm<coefficient> {
	affine constant;

public:
	typedef LinearForm<coefficient> linearPartType;

	linearPartType& linearPart() {
		return *this;
	}

	const linearPartType& linearPart() const {
		return *this;
	}

	affine& constantPart() {
		return constant;
	}

	const affine& constantPart() const {
		return constant;
	}

	bool isConstant() const {
		return linearPartType::isZero();
	}

	bool isZero() const {
		return linearPartType::isZero() && constantPart==0;
	}

	AffineForm(const std::vector<std::string>& variableNames) :
		linearPartType(variableNames), constant(0) {
	}

	AffineForm(const LinearForm<coefficient>& linear) :
		linearPartType(linear), constant(0) {
	}

	AffineForm(const LinearForm<coefficient>& linear, const affine& ct) :
		linearPartType(linear), constant(ct) {
	}

	AffineForm(const AffineForm& org) :
		linearPartType(org.linearPart()), constant(org.constantPart()) {
	}

	const std::vector<std::string>& getVariableNames() const {
		return linearPartType::getVariableNames();
	}

	template <class fromC, class fromA> AffineForm(
			const AffineForm<fromC, fromA>& f) :
		linearPartType(f), constant(f.constant) {
	}

	const coefficient& operator[](unsigned i) const {
		return (*((linearPartType*)this))[i];
	}

	coefficient& operator[](unsigned i) {
		return (*((linearPartType*)this))[i];
	}

	void addMonomial(unsigned i, const coefficient& x) {
		LinearForm<coefficient>::addMonomial(i, x);
	}
};

template <class coefficient, class affine> std::ostream& operator<<(
		std::ostream& out, const AffineForm<coefficient, affine>& form) {
	out << "(+ " << (form.linearPart()) << " " << form.constantPart() << ")";
	return out;
}

template <class coefficient, class affine> inline AffineForm<coefficient, affine> operator+(
		const AffineForm<coefficient, affine> &a,
		const AffineForm<coefficient, affine> &b) {
	return AffineForm<coefficient, affine>(a.linearPart()+b.linearPart(),
			a.constantPart()+b.constantPart());
}

template <class coefficient, class affine> inline AffineForm<coefficient, affine> operator-(
		const AffineForm<coefficient, affine> &a,
		const AffineForm<coefficient, affine> &b) {
	return AffineForm<coefficient, affine>(a.linearPart()-b.linearPart(),
			a.constantPart()-b.constantPart());
}

template <class coefficient, class affine, class scalar> inline AffineForm<coefficient, affine> operator*(
		const AffineForm<coefficient, affine> &a, const scalar &k) {
	return AffineForm<coefficient, affine>(a.linearPart()*k, a.constantPart()*k);
}

template <class coefficient, class affine, class scalar> inline AffineForm<coefficient, affine> operator/(
		const AffineForm<coefficient, affine> &a, const scalar &k) {
	return AffineForm<coefficient, affine>(a.linearPart()/k, a.constantPart()/k);
}

template <class coefficient, class affine> inline AffineForm<coefficient, affine>& operator+=(
		AffineForm<coefficient, affine> &a,
		const AffineForm<coefficient, affine> &b) {
	a.linearPart() += b.linearPart();
	a.constantPart() += b.constantPart();
	return a;
}

template <class coefficient, class affine, class scalar> inline AffineForm<coefficient, affine>& operator*=(
		AffineForm<coefficient, affine> &a, const scalar &k) {
	a.linearPart() *= k;
	a.constantPart() *= k;
	return a;
}

struct NonlinearProduct : std::exception {
	NonlinearProduct() {
	}

	const char *what() {
		return "nonlinear product in linear expression";
	}
};

template <class coefficient> inline AffineForm<coefficient, coefficient> operator*(
		const AffineForm<coefficient, coefficient> &a,
		const AffineForm<coefficient, coefficient> &b) {
	if (a.isConstant()) {
		return b * a.constantPart();
	} else if (b.isConstant()) {
		return a * b.constantPart();
	} else {
		throw NonlinearProduct();
	}
}

template <class coefficient> inline AffineForm<coefficient, coefficient> operator/(
		const AffineForm<coefficient, coefficient> &a,
		const AffineForm<coefficient, coefficient> &b) {
	if (b.isConstant()) {
		return a / b.constantPart();
	} else {
		throw NonlinearProduct();
	}
}

template <class coefficient> inline AffineForm<coefficient, coefficient> operator*=(
		AffineForm<coefficient, coefficient> &a,
		const AffineForm<coefficient, coefficient> &b) {
	return (a = a*b);
}

template <class coefficient, class affine> inline AffineForm<coefficient, affine> operator-(
		AffineForm<coefficient, affine> a) {
	return AffineForm<coefficient, affine>(-a.linearPart(), -a.constantPart());
}
#endif
