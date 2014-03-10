
#ifndef _INTERVAL_H
#define _INTERVAL_H

#include <list>

#include "math/numeric/comp.h"
#include "math/scalar_types/scalar_with_infinity.h"

namespace math {
namespace numeric {

template<class scalar_type>
class interval{
public:
	typedef scalar_with_infinity<scalar_type> gen_type;
	interval() {
		my_lower = scalar_with_infinity<scalar_type>::neg_infty();
		my_upper = scalar_with_infinity<scalar_type>::pos_infty();
	}
	interval(gen_type lu) :
		my_lower(lu), my_upper(lu) {
	}
	interval(gen_type l, gen_type u) :
		my_lower(l), my_upper(u) {
	}
	interval(scalar_type l, scalar_type u) :
		my_lower(scalar_with_infinity<scalar_type> (l)),
				my_upper(scalar_with_infinity<scalar_type> (u)) {
	}
	~interval(){}
	gen_type lower() const{
		return my_lower;
	}
	gen_type upper() const{
		return my_upper;
	}
	void set_lower(scalar_type l){
		my_lower = scalar_with_infinity<scalar_type>(l);
	}
	void set_upper(scalar_type u){
		my_upper = scalar_with_infinity<scalar_type>(u);
	}
	void set_empty(){
		my_lower = scalar_with_infinity<scalar_type>::pos_infty();
		my_upper = scalar_with_infinity<scalar_type>::neg_infty();
	}
	scalar_type size() const{
		if(my_upper.is_infinity() || my_lower.is_infinity())
			throw std::runtime_error("interval: interval size infinite");
		if (is_empty())
			throw std::runtime_error("interval: requested size of empty interval");
		return my_upper.get_val() - my_lower.get_val();
	}
	bool is_empty() const{
		if(my_lower.is_infinity() || my_upper.is_infinity())
			return (my_lower > my_upper);
		if(is_GT(my_lower.get_val(),my_upper.get_val()))
			return true;
		else
			return false;
	}
	bool is_finite() const {
		if(my_lower.is_infinity() || my_upper.is_infinity())
			return false;
		else
			return true;
	}
	void intersect(const interval<scalar_type>& intv) {
		// use strict comparisons if either one is infinite and
		// numerical comp otherwise
		if (my_lower.is_infinity() || intv.my_lower.is_infinity()) {
			if (my_lower < intv.my_lower)
				my_lower = intv.my_lower;
		} else {
			if (is_LT(my_lower.get_val(), intv.my_lower.get_val()))
				my_lower = intv.my_lower;
		}
		if (my_upper.is_infinity() || intv.my_upper.is_infinity()) {
			if (my_upper > intv.my_upper)
				my_upper = intv.my_upper;
		} else {
			if (is_GT(my_upper.get_val(), intv.my_upper.get_val()))
				my_upper = intv.my_upper;
		}
	}
private:
	scalar_with_infinity<scalar_type> my_lower;
	scalar_with_infinity<scalar_type> my_upper;
};

/** Intersect two lists of intervals.
 *
 * Does not return any empty intervals. */
template<typename scalar_type>
std::list<interval<scalar_type> > intersect(const std::list<interval<scalar_type> >& L1, const std::list<interval<scalar_type> >& L2) {
	typedef std::list<interval<scalar_type> > intv_list;
	intv_list L;

	for (typename intv_list::const_iterator i1=L1.begin();i1!=L1.end();++i1) {
		for (typename intv_list::const_iterator i2=L2.begin();i2!=L2.end();++i2) {
			interval<scalar_type> I=*i1;
			I.intersect(*i2);
			if (!I.is_empty())
				L.push_back(I);
		}
	}

	return L;
}

}
} // end of namespace math::numeric

template<typename scalar_type>
std::ostream& operator<<(std::ostream& os, const math::numeric::interval<
		scalar_type>& intv) {
   os << "[" << intv.lower() << "," << intv.upper() << "]";
   return os;
}

#endif /* _INTERVAL_H */
