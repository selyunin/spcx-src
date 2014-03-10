/***************************************************************************
 *   Copyright (C) 2004 by Goran Frehse                                    *
 *   gfrehse@localhost                                                     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "math/vdom/index_to_index_bimap.h"

#include <iostream>
#include <vector>
#include <stdexcept>

//using namespace Parma_Polyhedra_Library;
//using namespace IO_Operators;
using namespace std;

// --------------------------------------------------------------

// Mapping function for shifting dimensions of ccvs and cvs
// taken from index_to_index_bimap.hh


index_to_index_bimap::index_to_index_bimap() :
	mymax(0), mydommax(0), _removes_space_dimensions(false), _removes_space_dimensions_is_up_to_date(false) {
}

void index_to_index_bimap::swap_assign(size_t x1, size_t x2, size_t y, size_t n) {
	// swaps the variables in positions x1 ... x2 with y ... y2=x2-x1+y)
	// in a space of size n
	mymap.clear();
	// initialize map
	for (size_t i=0; i<n; ++i) {
		mymap.insert(Map::value_type(i, i));
	}
	for (size_t x=x1; x<=x2; ++x) {
		mymap[x]=y;
		mymap[y]=x;
		++y;
	}
	mymax=n-1;
	mydommax=mymax;
	_removes_space_dimensions_is_up_to_date=false;
}

void index_to_index_bimap::move_assign(size_t x1, size_t x2, size_t y, size_t n) {
	// moves the variables x1 ... x2 to position y in a map of size n of the original state set
	mymap.clear();
	// initialize map
	for (size_t i=0; i<n; ++i) {
		mymap.insert(Map::value_type(i, i));
	}
	// move the x values to the position of y
	size_t y1=y;
	for (size_t x=x1; x<=x2; ++x) {
		mymap[x]=y1;
		++y1;
	}
	if (y>=x1) {
		// shift the remaining ones to the left
		y1=x1;
		for (size_t i=x2+1; i+x1<=y+x2; ++i) {
			mymap[i]=y1;
			++y1;
		};
	} else {
		// shift the remaining ones to the right
		y1=y+x2-x1+1;
		for (size_t i=y; i<x1; ++i) {
			mymap[i]=y1;
			++y1;
		};
	}
	mymax=n-1;
	mydommax=mymax;
	_removes_space_dimensions_is_up_to_date=false;
}

bool index_to_index_bimap::has_empty_codomain() const {
	return mymap.empty();
}

size_t index_to_index_bimap::max_in_codomain() const {
	if (has_empty_codomain())
		throw std::runtime_error("PFunction::max_in_codomain() called when has_empty_codomain()");
	return mymax;
}

bool index_to_index_bimap::maps(size_t x, size_t& y) const {
	Map::const_iterator i = mymap.find(x);
	if (i != mymap.end()) {
		y = (*i).second;
		return true;
	} else
		return false;
}

//  size_t
//	Parma_Polyhedra_Library::dimension_type // NOTE:: size_t yiels compilation error!
size_t index_to_index_bimap::get_map(size_t x) const {
	if (has_empty_codomain())
		throw std::runtime_error("PFunction::get_map() called"
				" when has_empty_codomain()");
	Map::const_iterator i = mymap.find(x);
	if (i==mymap.end())
		throw std::runtime_error("PFunction::get_map() called"
				" with non-existent x");
	return i->second;
}

size_t // NOTE:: size_t yiels compilation error!
index_to_index_bimap::get_premap(size_t y) const {
	if (mymap.empty())
		throw std::runtime_error("PFunction::get_premap() called"
				" when map is empty");
	for (Map::const_iterator i = mymap.begin(); i!=mymap.end(); ++i) {
		if (i->second==y)
			return i->first;
	}
	throw std::runtime_error("PFunction::get_premap() called "
			" with non-existent y");
	return mymap.begin()->first; // just return a dummy value
}

bool index_to_index_bimap::in_domain(const size_t& x) const {
	return mymap.find(x) != mymap.end();
}

bool index_to_index_bimap::in_domain(const size_t& x, size_t& y) const {
	Map::const_iterator it = mymap.find(x);
	if (it != mymap.end()) {
		y = it->second;
		return true;
	} else
		return false;
}

bool index_to_index_bimap::in_codomain(const size_t& y) const {
	for (Map::const_iterator i = mymap.begin(); i!=mymap.end(); ++i) {
		if (i->second==y)
			return true;
	}
	return false;
}

void index_to_index_bimap::fill_up_to(size_t newdim) {
	// add the remaining variables that are not already mapped by *this to map
	// in no particular order (first free slots)
	// (they have to be in it, because otherwise they will be removed)
	// just fill in the slots that are not yet taken

	for (size_t i=0; i<newdim; ++i) {
		//      if (!in_domain(i)) // map those that are not yet in the map
		if (!in_codomain(i)) // map those that are not yet in the map
		{
			for (size_t jj=0; jj<newdim; ++jj) {
				//          if (!in_codomain(jj))
				if (!in_domain(jj)) {
					//            insert(i,jj);
					insert(jj, i);
					break;
				}
			}
		}
	}
	_removes_space_dimensions_is_up_to_date=false;
}

void index_to_index_bimap::print(std::ostream& s) const {

	if (has_empty_codomain())
		s << "empty" << std::endl;
	else
		for (Map::const_iterator i = mymap.begin(); i != mymap.end(); ++i) {
			if (i != mymap.begin() ) {
				s << ",";
			}
			s << (*i).first << "->" << (*i).second;
		}
}

void index_to_index_bimap::insert(size_t x, size_t y) {
	//	mymap[x]=y;
	std::pair<Map::iterator, bool> stat = mymap.insert(Map::value_type(x, y));
	if (!stat.second)
		throw std::runtime_error("PFunction::insert(x, y) called with `x' already in domain");
	if (y > mymax)
		mymax = y;
	if (x > mydommax)
		mydommax = x;
	_removes_space_dimensions_is_up_to_date=false;
}

bool index_to_index_bimap::removes_space_dimensions() const {
	// return true of variables are disappearing
	if (_removes_space_dimensions_is_up_to_date) {
		return _removes_space_dimensions;
	} else {
		index_to_index_bimap& pf = const_cast<index_to_index_bimap&>(*this);

		if (mydommax > mymax)
			pf._removes_space_dimensions=true;
		else {
			pf._removes_space_dimensions=false;
			// check if all dimensions in between are present
			// and if index_to_index_bimap is many-to-one
			vector <int> source_id(mymax+1, -1);
			size_t y;
			for (size_t x = 0; x<=mydommax && !pf._removes_space_dimensions; ++x) {
				if (maps(x, y)) {
					if (source_id[y]==-1)
						source_id[y]=x;
					else if (source_id[y]!=(int)x)
						pf._removes_space_dimensions=true;
				} else
					pf._removes_space_dimensions=true;
			}
		}
		pf._removes_space_dimensions_is_up_to_date=true;
		return pf._removes_space_dimensions;
	}
}

size_t index_to_index_bimap::max_in_domain() const {
	return mydommax;
}

index_to_index_bimap index_to_index_bimap::identity(size_t d) {
	index_to_index_bimap M;
	/*
	Map::iterator it=M.mymap.begin();
	for (size_t i=0;i<d;++i) {
		it=M.mymap.insert(it,make_pair(i,i));
		//++it;
	}
	if (d>0)
		M.mymax=d-1;
	else
		M.mymax=0;
	M.mydommax=M.mymax;
	M._removes_space_dimensions=false;
	M._removes_space_dimensions_is_up_to_date=true;
	*/
	for (size_t i=0;i<d;++i) {
		M.insert(i,i);
	}
	return M;
}

std::ostream& operator<<(std::ostream& os,
		const index_to_index_bimap& map) {
	map.print(os);
	return os;
}

