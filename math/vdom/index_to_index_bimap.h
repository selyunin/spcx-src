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

#ifndef GUARD_index_to_index_bimap_h
#define GUARD_index_to_index_bimap_h

#include <map>
#include <iostream>

typedef unsigned int dimension_t;

/** Bidirectional map from index to index, i.e., \f$f:\{0,...,n-1\} \rightarrow \{0,...,m-1\}\f$ for some n and m.
 * If \f$y_i=f(x_i)\f$, we call \f$x_i\f$ the source dimension and \f$y_i\f$ the target dimension.
 * The main application of the map is in reordering the dimensions of sets in \f$R^n\f$.
 * Example: Let \f$x=(x_0,x_1,x_2...,x_{n-1})\f$ be a point in \f$R^n\f$, then swapping
 * dimensions 0 and 1 would result in the point \f$x=(x_1,x_0,x_2,...,x_{n-1})\f$.
 * \note Implements a "Partial Function" as used in the Parma Polyhedra Library (PPL).
 * \see Parma Polyhedra Library documentation.
 */
class index_to_index_bimap {

public:
	index_to_index_bimap(); // : max(0)
	virtual ~index_to_index_bimap() {
	}
	;

	/** Swaps the dimensions in positions x1,...,x2 with y,...,y2, where y2=x2-x1+y,
	 * in a space of dimension n. */
	void swap_assign(size_t x1, size_t x2, size_t y, size_t n);

	/** Moves the dimensions in positions x1,...,x2 to y,...,y2, where y2=x2-x1+y,
	 * shifting the remaining variables (left or right) such that their order is preserved,
	 * in a space of dimension n. */
	void move_assign(size_t x1, size_t x2, size_t y, size_t n);

	/** Returns true if the map is empty. */
	bool has_empty_codomain() const;

	/** Returns the highest dimension resulting from applying the map. */
	size_t max_in_codomain() const;

	/** Returns the highest dimension for which the map is defined. */
	size_t max_in_domain() const;

	/** Returns true if the map is defined for argument \p x.*/
	bool in_domain(const size_t& x) const;

	/** Returns true if the map is defined for argument \p x, and assigns
	 * its mapped value to y.*/
	bool in_domain(const size_t& x, size_t& y) const;

	/** Returns true if some dimension in the map is mapped to \p y.*/
	bool in_codomain(const size_t& y) const;

	/** Returns the dimension to which \p x is mapped. */
	size_t get_map(size_t x) const;

	/** Returns the dimension that is mapped to \p y, i.e., the inverse of the map. */
	size_t get_premap(size_t y) const;

	/** Returns true if \p x is mapped to some target dimension and assigns it to \p y.
	 * If \p x has no target, then \p y remains unchanged. */
	bool maps(size_t x, size_t& y) const;

	/** Adds mappings for the dimensions out of the domain 0,...,n-1 that are not already in the map.
	 * The target for each dimension \f$x_i\f$ is iteratively chosen to be the lowest dimension
	 * that is not yet a target of the map.
	 * \note This function is useful when working with the PPL because any dimension that is not in the codomain of the
	 * map is removed by existential quantification, reducing the dimension of the polyhedron.
	 */
	void fill_up_to(size_t newdim);

	/** Output the map as a stream of characters in the form \f$x_0->f(x_0),x_1->f(x_1),...\f$. */
	void print(std::ostream& s) const;

	/** Insert \f$x \rightarrow y\f$ into the map.
	 * An exception is thrown if x is already in the codomain of the map, i.e., already has a target. */
	void insert(size_t x, size_t y);

	/** Returns true if there is a dimension in 0,...,d that has no target, or whose target is also the target of another dimension.
	 * d is the highest dimension that has a target in the map, i.e., the highest dimension of
	 * the domain.
	 * \note When applying a reordering based on *this to a PPL object, the number of dimensions does
	 * not decrease if \p removes_space_dimensions() returns true.
	 */
	bool removes_space_dimensions() const;

	/** Returns the identity map for d indices. */
	static index_to_index_bimap identity(size_t d);

private:
	typedef std::map<size_t, size_t, std::less<size_t> > Map;
	Map mymap;
	size_t mymax;
	size_t mydommax;
	bool _removes_space_dimensions;
	bool _removes_space_dimensions_is_up_to_date;
};

/**
 * Output as a stream of characters. Calls the print method.
 */
std::ostream& operator<<(std::ostream& os, const index_to_index_bimap& map);

#endif

