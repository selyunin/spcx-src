#ifndef LABEL_CACHE_H_
#define LABEL_CACHE_H_

/*********************************************************************
 * label_cache.h													 *
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
 * 			Goran Frehse (goran.frehse@imag.fr)					  	 *
 *********************************************************************/

#include <map>
#include <set>

namespace exploration_graph {
class EG_NFA_interface;
}

namespace finite_automaton {

/**
 * label_cache assigns label id's to the labels of an automaton alphabet,
 * retrieves label id's from the labels and vice-versa. The functionality
 * is carried out using a map, label_id_to_str_map.
 */
class label_cache {

public:

	/** type of label_id */
	typedef int label_id_type;

	/** Get label for label-id. */
	static std::string getLabel(const label_id_type label_id);

	/** Get label-id for the label. */
	static label_id_type getLabelId(const std::string label);

	/** Return all label-id's in the map. */
	static std::vector<label_id_type> getLabelIds();

private:

	static label_id_type label_id; /*!< label-id declaration */

	static std::map<label_id_type, std::string> label_id_to_str_map; /*!< map declaration */

	/** Initialize label_cache. */
	static std::map<label_id_type, std::string> init_label_cache();

	/** Add label to label_id_to_str_map and return the label-id. */
	static label_id_type addLabel(std::string label);

	/** It uses addLabel() method. */
	friend class Automaton;

	friend class exploration_graph::EG_NFA_interface;

};

std::map<label_cache::label_id_type, std::string> label_cache::init_label_cache() {
	std::map<label_id_type, std::string> l_cache;
	l_cache.insert(std::pair<label_id_type, std::string>(0, "EPS"));
	return l_cache;
}

std::string label_cache::getLabel(label_id_type label_id) {
	std::map<label_id_type, std::string>::const_iterator it =
			label_id_to_str_map.find(label_id);
	if (it != label_id_to_str_map.end())
		return it->second;
	else
		return "";
}

label_cache::label_id_type label_cache::getLabelId(std::string label) {
	for (std::map<label_id_type, std::string>::const_iterator it =
			label_id_to_str_map.begin(); it != label_id_to_str_map.end();
			++it) {
		if (it->second == label)
			return it->first;
	}
	return -1;
}

label_cache::label_id_type label_cache::addLabel(std::string label) {
	std::map<label_id_type, std::string>::iterator it;
	label_id++;
	label_id_to_str_map.insert(
			std::pair<label_id_type, std::string>(label_id, label));
	return label_id;
}

std::vector<label_cache::label_id_type> label_cache::getLabelIds() {
	std::vector<label_id_type> labelIds_vector;
	for (std::map<label_id_type, std::string>::const_iterator it =
			label_id_to_str_map.begin(); it != label_id_to_str_map.end();
			++it) {
		labelIds_vector.push_back(it->first);
	}
	return labelIds_vector;
}
}

#endif /* LABEL_CACHE_H_ */
