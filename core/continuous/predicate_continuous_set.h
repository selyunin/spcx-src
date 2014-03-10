#ifndef GUARD_predicate_continuous_set_h
#define GUARD_predicate_continuous_set_h
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

#include <boost/enable_shared_from_this.hpp>
#include "core/continuous/continuous_set.h"

/** Forward declarations */
namespace tree {
	class node;
}

namespace continuous {

class predicate_continuous_set : public continuous_set,
	public boost::enable_shared_from_this<predicate_continuous_set> {
public:
	typedef boost::shared_ptr<predicate_continuous_set> ptr;
	typedef boost::shared_ptr<const predicate_continuous_set> const_ptr;

	typedef tree::node predicate_type;
	typedef boost::shared_ptr<predicate_type> predicate_type_ptr;

	/** Return a shared_ptr to *this. */
	virtual continuous_set::ptr get_ptr();

	/** Return a shared_ptr to const *this. */
	virtual continuous_set::const_ptr get_const_ptr() const;

	/** Initialize with an element. */
	explicit predicate_continuous_set(predicate_type_ptr root_node);

	virtual ~predicate_continuous_set();

	virtual void set_predicate(predicate_type_ptr new_pred);

	/** Creates a shallow copy of *this. */
	virtual predicate_continuous_set* clone() const;

	virtual predicate_continuous_set* create_universe() const;
	virtual predicate_continuous_set* create_empty() const;

	virtual int get_memory() const;
	virtual const variable_id_set& get_variable_ids() const;

	virtual unsigned int get_dim() const;

	virtual math::tribool is_empty() const;

	virtual void embed_variables(const variable_id_set& id_set);
	virtual void existentially_quantify_variables(const variable_id_set& id_set);

	/** Set the primedness of the variables with primedness of degree \p d to degree \p p.
	 */
	virtual void reassign_primedness(unsigned int d, unsigned int p = 0);

	/** Increase the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, increase all. */
	virtual void increase_primedness(unsigned int d = 0);

	/** Decrease the primedness of the variables with primedness of degree \p d by 1.
	 * If d is 0, decrease all. */
	virtual void decrease_primedness(unsigned int d = 0);

	/** Intersect with cset. */
	virtual void intersection_assign(const predicate_continuous_set& cset);

	/** Returns the set in predicate form.*/
	virtual continuous_set_predicate::ptr get_predicate() const;

	/** Accept a visitor. */
	virtual void accept(dispatching::dispatcher<continuous_set_typelist>& d) const;

	class output_format {
public:
		output_format() { // default values
			preamble="";
			epilogue="";
			nullstring="true";
		}
		;
		static output_format matlab() {
			output_format f;
			f.preamble="[";
			f.epilogue="]";
			f.nullstring="";
			return f;
		}
		;
		static output_format space_separated() {
			output_format f;
			f.preamble="";
			f.epilogue="";
			f.nullstring="true";
			return f;
		}
		;
		std::string preamble; // written before matrix
		std::string epilogue; // written after matrix
		std::string nullstring; // written if null pointer
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
	predicate_continuous_set();

private:
	predicate_type_ptr my_predicate;
	static output_format my_output_format;
};

}

#endif

