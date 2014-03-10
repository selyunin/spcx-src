#include "core/discrete/discrete_set_stl_set.h"

#include <algorithm>
#include "core/discrete/discrete_set_operators.h"

namespace discrete {

using namespace std;

discrete_set_stl_set::discrete_set_stl_set() {
}

discrete_set_stl_set::~discrete_set_stl_set() {
	// Note : Don't destruct _myset here by hand -- the destructor is called automatically.
}

discrete_set_stl_set* discrete_set_stl_set::create_empty() const {
	return new discrete_set_stl_set();
}

discrete_set_stl_set* discrete_set_stl_set::clone() const {
	// Create a new object of the same type
	discrete_set_stl_set* ip = create_empty();

	for (container_type::const_iterator it = _myset.begin(); it != _myset.end(); ++it) {
		object_type element_clone(*it); // no need to clone, it's just a location_constraint_set
		ip->add(element_clone);
	}
	return ip;
}

const discrete_set_stl_set::container_type& discrete_set_stl_set::get_stl_set() const {
	return _myset;
}

void discrete_set_stl_set::add(const object_type& e) {
	_myset.push_back(e);
}

std::size_t discrete_set_stl_set::get_size() const {
	return _myset.size();
}

/** Returns true if and only if *this contains loc.
 */
bool discrete_set_stl_set::contains(const object_type& loc) const {
	discrete_set_stl_set::ptr dp = discrete_set_stl_set::ptr(new discrete_set_stl_set());
	dp->add(loc);
	discrete_set::ptr d = compute_or_assign_difference(dp, get_const_ptr());
	return d->is_empty();
}

bool discrete_set_stl_set::is_empty() const {
	if (_myset.empty())
		return true;
	else {
		for (container_type::const_iterator it = _myset.begin(); it != _myset.end(); ++it) {
			if (!it->is_empty()) {
				return false; // found a nonempty subset
			}
		}
		return true;
	}
}

//bool discrete_set_stl_set::is_disjoint_from(const discrete_set::const_ptr& dsp) const {
//	throw std::runtime_error("missing implementation discrete_set_stl_set::is_disjoint_from");
//
//	/*	assert(dsp->get_discrete_set_type()==stl_set);
//	 // only compare sets of type stl_set
//
//	 const discrete_set_stl_set* dsssp = (discrete_set_stl_set*)dsp.get();
//	 set<object_type>::iterator i1 = _myset.begin();
//	 set<object_type>::iterator e1 = _myset.end();
//	 set<object_type>::iterator i2 = dsssp->get_stl_set().begin();
//	 set<object_type>::iterator e2 = dsssp->get_stl_set().end();
//
//	 set<object_type> out_set;
//	 set_intersection(i1, e1, i2, e2, inserter(out_set, out_set.begin()));
//
//	 return out_set.empty();
//	 */
//
////	set<object_type> out_set;
////	set_intersection(begin(), end(), dsp->begin(), dsp->end(), inserter(out_set, out_set.begin()));
////
////	return out_set.empty();
//}

/** Add loc to *this */
void discrete_set_stl_set::union_assign(const object_type& loc) {
	add(loc);
}

void discrete_set_stl_set::union_assign(const discrete_set::const_ptr& dsp) {
	//throw std::runtime_error("missing implementation discrete_set_stl_set::union_assign");

		container_type out_set;
//		set_union(begin(), end(), dsp->begin(), dsp->end(), inserter(out_set, out_set.begin()));
		copy(dsp->begin(), dsp->end(), inserter(out_set, out_set.begin()));
		_myset.swap(out_set);
}

void discrete_set_stl_set::remove_empty() {
	container_type::iterator it = _myset.begin();
	while (it != _myset.end()) {
		if (it->is_empty()) {
			it = _myset.erase(it);
		} else {
			++it;
		}
	}
}

void discrete_set_stl_set::intersection_assign(const object_type& loc) {
	for (container_type::iterator it = _myset.begin(); it != _myset.end(); ++it) {
		it->intersection_assign(loc);
	}
	remove_empty();
}

void discrete_set_stl_set::intersection_assign(const discrete_set_stl_set& dsp) {
	//throw std::runtime_error("missing implementation discrete_set_stl_set::intersection_assign");

	container_type new_set;
		// Intersect *this with every object in dsp
		for (container_type::const_iterator it = _myset.begin(); it != _myset.end(); ++it) {
			for (container_type::const_iterator jt = dsp._myset.begin(); jt != dsp._myset.end(); ++jt) {
				object_type x = *it;
				x.intersection_assign(*jt);
				if (!x.is_empty())
					new_set.push_back(x);
			}
		}
		_myset.swap(new_set);
}

void discrete_set_stl_set::existentially_quantify(automaton_id aut_id) {
	for (container_type::iterator it = _myset.begin(); it != _myset.end(); ++it) {
		it->existentially_quantify(aut_id);
	}
}

/** Remove loc from *this */
void discrete_set_stl_set::difference_assign(const object_type& loc) {
	throw std::runtime_error("missing implementation discrete_set_stl_set::difference_assign");

}

void discrete_set_stl_set::difference_assign(const discrete_set::const_ptr& dsp) {
	throw std::runtime_error("missing implementation discrete_set_stl_set::difference_assign");

	/*
	 assert(dsp->get_discrete_set_type()==stl_set);
	 // only compare sets of type stl_set

	 const discrete_set_stl_set* dsssp = (discrete_set_stl_set*)dsp.get();
	 set<object_type>::iterator i1 = _myset.begin();
	 set<object_type>::iterator e1 = _myset.end();
	 set<object_type>::iterator i2 = dsssp->get_stl_set().begin();
	 set<object_type>::iterator e2 = dsssp->get_stl_set().end();

	 std::set<object_type> tmp_set;
	 set_difference(i1, e1, i2, e2, inserter(tmp_set, tmp_set.begin()));
	 _myset.swap(tmp_set); // replace the current set with the result computed above
	 */
	//	container_type out_set;
	//	set_difference(begin(), end(), dsp->begin(), dsp->end(), inserter(out_set, out_set.begin()));
	//
	//	_myset.swap(out_set);

}

void discrete_set_stl_set::create_empty() {
	_myset.clear(); // makes the set empty.
}

void discrete_set_stl_set::print(ostream& os) const {
	if (_myset.size() > 1)
		os << "(";
	for (container_type::const_iterator j = _myset.begin(); j != _myset.end(); ++j) {
		if (j != _myset.begin())
			os << ") | (" << *j;
		else
			os << *j;
	};
	if (_myset.size() > 1)
		os << ")";
}

void discrete_set_stl_set::accept(dispatching::dispatcher<discrete_set_typelist>& d) const {
	d.dispatch(this);
}

}

// END
