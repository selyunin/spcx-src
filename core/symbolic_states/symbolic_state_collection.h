#ifndef SYMBOLIC_STATE_COLLECTION_H_
#define SYMBOLIC_STATE_COLLECTION_H_

#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include "math/vdom/variable.h"
#include "utility/shared_ptr_output.h"
#include "utility/printable.h"
#include "utility/simple_iterators/collection_base.h"

//#include "symbolic_state.h"

/** Forward declarations. */
namespace hybrid_automata {
class symbolic_state;
typedef boost::shared_ptr<symbolic_state> symbolic_state_ptr;
typedef boost::shared_ptr<const symbolic_state> symbolic_state_const_ptr;
class adapt_discrete_set_visitor;
class adapt_continuous_set_visitor;
}

namespace hybrid_automata {

/**
 * An enumerable collection of symbolic states.
 * The collection can be traversed using a forward iterator called symbolic_state_collection::const_iterator
 * (in STL-speak an "input iterator"). It supports dereferencing (*i), incrementation (++i), assignment (j=i),
 * and testing for equality (i==j, i!=j).
 *
 * \note { For symbolic representations like bdds etc. this class would be nonsense).
 * \todo { Think about what needs to be done for purely symbolic representations. }
 */
class symbolic_state_collection: public simple_iterators::collection_const_base<symbolic_state_ptr>,
		public printable,
		public boost::enable_shared_from_this<symbolic_state_collection> {
public:
	typedef boost::shared_ptr<symbolic_state_collection> ptr;
	typedef boost::shared_ptr<const symbolic_state_collection> const_ptr;

	/** Return a shared_ptr to *this. */
	ptr get_ptr();

	/** Return a shared_ptr to const *this. */
	const_ptr get_const_ptr() const;

	/** Virtual constructor, initialized to be empty. */
	virtual symbolic_state_collection* create() const = 0;

	/** Creates an identical copy of *this
	 *
	 * \return pointer to be new clone object.
	 */
	virtual symbolic_state_collection* clone() const = 0;

	/** Deep copy constructor. Adds copies of the symbolic states in \p sstate_set to *this.
	 */
	virtual void copy(const symbolic_state_collection::const_ptr& sstate_set);

	/** Virtual destructor. */
	virtual ~symbolic_state_collection() {
	}
	;

	/** Returns the ids of all variables over which the set is defined. */
	virtual variable_id_set get_variable_ids() const;

	/** Adds \p sstate to \p *this.
	 */
	virtual void add(const symbolic_state_ptr& sstate) = 0;

	/**
	 * \brief Utility function
	 *
	 * \return The memory used for representing this symbolic state collection.
	 */
	virtual std::size_t get_memory() const = 0;

	/**
	 * \brief Basic boolean set manipulation function
	 *
	 * checks for the emptyness of the set.
	 * \return \p true if and only if the symbolic_state_collection is empty
	 */
	virtual bool is_empty() const = 0;

	/**
	 * Returns the size of *this.
	 */
	virtual unsigned int size() const = 0;

	/**
	 * \brief Auxiliary set manipulation function.
	 *
	 * \param sstate_set pointer to a symbolic_state_collection object.
	 * \return \p true if and only if this set is disjoint
	 * with the set pointed by the passed parameter.
	 */
	virtual bool is_disjoint_from(const symbolic_state_collection::ptr& sstate_set) const;

	/**
	 * \brief Auxiliary set manipulation function.
	 *
	 * The default implementation corresponds to difference assign followed by checking emtpiness,
	 * but it may be overridden by a more efficient implementation.
	 *
	 * \return \p true if and only if \p *this contains \p *ps.
	 */
	virtual bool contains(const symbolic_state_collection::ptr& sstate_set) const= 0;

	/** \brief Makes the collection empty.
	 */
	virtual void clear() = 0;

	/**
	 * \brief Assigns to \p *this the union of \p *this and \p *sstate_set.
	 *
	 * The default implementation enumerates all symbolic states in sstate_set
	 * and adds them to *this. */
	virtual void union_assign(const symbolic_state_collection::const_ptr& sstate_set);

	/**
	 * \brief Assigns to \p *this the intersection of \p *this and \p *sstate_set.
	 */
	virtual void intersection_assign(const symbolic_state_collection::ptr& sstate_set) = 0;

	/**
	 * \brief Assigns to \p *this the symbolic states of \p *this that are
	 * not in \p *sstate_set.
	 */
	virtual void difference_assign(const symbolic_state_collection::const_ptr& sstate_set) = 0;

	/** Adapt the discrete sets in *this according to v.
	 * Return true if adaptation successful. */
	virtual bool accept(adapt_discrete_set_visitor& v) = 0;

	/** Adapt the continuous sets in *this according to v.
	 * Return true if adaptation successful. */
	virtual bool accept(adapt_continuous_set_visitor& v) = 0;

	/** Prints the sequence of symbolic states, separated by a string
	 * according to the output_format.
	 * The default format is
	 * s1
	 * s2
	 * ...
	 */
	virtual void print(std::ostream& os) const;

	class output_format {
	public:
		output_format() { // default values
			preamble = "{";
			epilogue = "}";
			element_separator = ",";
			empty_signal = "false";
		}
		;
		static output_format lf_separated() {
			output_format f;
			f.preamble = "";
			f.epilogue = "";
			f.element_separator = "\n";
			f.empty_signal = "";
			return f;
		}
		;
		static output_format matlab() {
			output_format f = lf_separated();
			// @todo
			return f;
		}
		;
		static output_format JVX() {
			output_format f = lf_separated();
			return f;
		}
		;
		std::string preamble; // written before matrix
		std::string epilogue; // written after matrix
		std::string element_separator; // written between elements
		std::string empty_signal; // written if empty
	};

	static const output_format& get_output_format() {
		return my_output_format;
	}
	;
	static void set_output_format(output_format new_format) {
		my_output_format = new_format;
	}
	;

private:
	static output_format my_output_format;
};

}

#endif /*SYMBOLIC_STATE_COLLECTION_H_*/
