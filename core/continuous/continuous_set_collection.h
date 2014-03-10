#ifndef GUARD_continuous_set_collection_h
#define GUARD_continuous_set_collection_h
/***************************************************************************
 *   Copyright (C) 2008 by Goran Frehse   *
 *   goran.frehse@imag.fr   *
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

#include <list>
#include "core/continuous/continuous_set.h"

namespace continuous {

/**
 * Container class for list of continuous sets. Helps to store the results
 * of post continuous and post discrete operators.
 * @author Rajarshi
 *
 * @todo some member functions needs to be documented
 */
class continuous_set_collection : public continuous_set {
public:
	typedef boost::shared_ptr<continuous_set_collection> ptr;
	typedef boost::shared_ptr<const continuous_set_collection> const_ptr;

	typedef continuous_set::ptr element_type;
	typedef std::list<element_type> container_type;
	typedef container_type::iterator iterator;
	typedef container_type::const_iterator const_iterator;

	/** Intermediate constructor -- not in a consistent state. */
	continuous_set_collection();

	/** Initialize with an element. */
	explicit continuous_set_collection(continuous_set::ptr element);

	virtual ~continuous_set_collection();

	/** Create a deep copy. */
	virtual continuous_set_collection* clone() const;

	virtual continuous_set_collection* create_universe() const;
	virtual continuous_set_collection* create_empty() const;

	virtual continuous_set::const_ptr default_element() const;
	virtual iterator begin();
	virtual iterator end();
	virtual const_iterator begin() const;
	virtual const_iterator end() const;

	virtual unsigned int size() const;

	virtual int get_memory() const;
	virtual const variable_id_set& get_variable_ids() const;

	virtual unsigned int get_dim() const;

	/** The collection is by construction empty iff the first element is empty. */
	virtual math::tribool is_empty() const;

	virtual void embed_variables(const variable_id_set& id_set);
	virtual void existentially_quantify_variables(const variable_id_set& id_set);

	virtual math::tribool element_wise_contains(continuous_set::const_ptr p) const;
	virtual void delete_if_contained_in(continuous_set::const_ptr p);

	/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
	 */
	virtual void reassign_primedness(unsigned int d, unsigned int p = 0);

	/** Increase the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, increase all. */
	virtual void increase_primedness(unsigned int d = 0);

	/** Decrease the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, decrease all. */
	virtual void decrease_primedness(unsigned int d = 0);



	/** Insert p into the container. By default or with redundancy_check=true,
	 * the element is only inserted if it is not contained already; also,
	 * all elements contained in p are removed. */
	virtual void insert(continuous_set::ptr p, bool redundancy_check=true);

	/** Returns the set in predicate form.*/
	virtual continuous_set_predicate::ptr get_predicate() const;

	/** Accept a visitor. */
	virtual void accept(dispatching::dispatcher<continuous_set_typelist>& d) const;

	class output_format {
	public:
		output_format() { // default values
			preamble="{";
			epilogue="}";
			element_separator=",";
		}
		;
		static output_format matlab() {
			output_format f;
			f.preamble="[";
			f.epilogue="]";
			f.element_separator=",";
			return f;
		}
		;
		static output_format space_separated() {
			output_format f;
			f.preamble="";
			f.epilogue="";
			f.element_separator=" ";
			return f;
		}
		;
		std::string preamble; // written before matrix
		std::string epilogue; // written after matrix
		std::string element_separator; // written between elements
	};
	static const output_format& get_output_format() {
		return my_output_format;
	}
	;
	static void set_output_format(output_format new_format) {
		my_output_format=new_format;
	}
	;


	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const;

protected:
	/** Insert element at end of list without emptiness checking. */
	virtual void push_back(continuous_set::ptr p);

private:
	container_type my_container;
	static output_format my_output_format;
};

}

#endif

