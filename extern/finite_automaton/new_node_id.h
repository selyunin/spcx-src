#ifndef NEW_NODE_ID_H_
#define NEW_NODE_ID_H_

/*********************************************************************
 * new_node_id.h													 *
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
 *																	 *
 * @authors: Manish Goyal (manish.goyal@imag.fr)					 *
 * 			 Goran Frehse (goran.frehse@imag.fr)					 *
 *********************************************************************/

namespace finite_automaton {

/** new_node_id class is used to assign id's to the nodes. */

class new_node_id {

public:

	typedef unsigned int node_id_type;

	/** Create a new 'within the range' id for the node */
	node_id_type static create_id() {
		node_id_type max = highestid + 1;
		if (max > highestid) {
			++highestid;
			return highestid;
		} else
			throw std::runtime_error(
					"@create_id: The maximum limit for node id has reached.");
	}

private:

	/** Sets the highest id to the new id */
	static node_id_type set_id(const node_id_type id) {
		if (highestid < id) {
			highestid = id;
			return highestid;
		} else
			return create_id();
	}

	/** highest-id that has been assigned to the last node. */
	static node_id_type highestid;
};

}

#endif /* NEW_NODE_ID_H_ */
