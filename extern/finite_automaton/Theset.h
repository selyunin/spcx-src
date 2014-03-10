#ifndef THESET_H_
#define THESET_H_

/*********************************************************************
 * Theset.h													 		 *
 * 																	 *
 * Copyright (C) 2012 by Verimag Research Lab, Grenoble, France.  	 *
 *																	 *
 * This program is free software; you can redistribute it and/or  	 *
 * modify it under the terms of the GNU Lesser General Public     	 *
 * License as published by the Free Software Foundation; either   	 *
 * version 2.1 of the License, or (at your option) any later version.*
 * 												  					 *
 * This program is distributed in the hope that it will be useful,	 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 	 *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU *
 * Lesser General Public License for more details.			  	 	 *
 *																	 *
 * You should have received a copy of the GNU Lesser General Public	 *
 * License along with this library; if not, write to the 			 *
 * Free	Software Foundation, Inc., 									 *
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.	 *
 *
 * @cite http://www.daniweb.com/software-development/cpp/threads/    *
 * 159859 															 *
 *																	 *
 * @authors: Manish Goyal (manish.goyal@imag.fr)					 *
 * 			 Goran Frehse (goran.frehse@imag.fr)					 *
 *********************************************************************/

#include <iostream>
#include <vector>

namespace finite_automaton {
/**
 * Theset class has a set, represented by an array of T type
 * elements of dimension Dim. It is used to compute the
 * power-set (all subsets of a set) of this set.
 */

template<typename T>
class Theset {
private:

	/** Dimension of Arr */
	int Dim;

	/** Arr represents the set */
	T *Arr;

	/** Assists subset computation */
	bool *BitMask;
public:
	/** Create a set from an array of type T */
	Theset(T *arr, int dim);

	/** Destructor */
	~Theset(void);

	/** Compute the next subset of the set Arr */
	std::set<T> NextSubset(void);
};

template<typename T>
std::set<T> Theset<T>::NextSubset(void) {
	std::set<T> nextSubSet;
	int i = Dim - 1;
	while (!BitMask[i] && i >= 0) {
		i--;
		if (i < 0) {
			break;
		}
	}

	if (i >= 0) {
		BitMask[i] = !BitMask[i];
	}

	for (int j = i + 1; j < Dim; j++) {
		BitMask[j] = !BitMask[j];
	}

	for (int j = 0; j < Dim; j++)
		if (BitMask[j])
			nextSubSet.insert(Arr[j]);

	return nextSubSet;
}

template<typename T>
Theset<T>::Theset(T *arr, int dim) {
	Arr = new T[dim];
	BitMask = new bool[dim];
	Dim = dim;
	for (int i = 0; i < dim; i++) {
		Arr[i] = arr[i];
		BitMask[i] = true;
	}
}

template<typename T>
Theset<T>::~Theset(void) {
	delete[] Arr;
	delete[] BitMask;
}
}
#endif /* H_THESET */
