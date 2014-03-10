#ifndef SYMBOLIC_STATE_H_
#define SYMBOLIC_STATE_H_

#include "boost/shared_ptr.hpp"
#include "math/vdom/variable.h"
//#include "../../utility/shared_ptr_output.h"
//#include "../continuous/continuous_set.h"
//#include "../discrete/discrete_set.h"
#include "utility/printable.h"

/** Forward declaration of classes used in header file. */
namespace discrete {
class discrete_set;
typedef boost::shared_ptr<discrete_set> discrete_set_ptr;
typedef boost::shared_ptr<const discrete_set> discrete_set_const_ptr;
}
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;
}

namespace hybrid_automata {

/**
 * A symbolic state is a pair of a discrete_set \f$D\f$ and continuous_set \f$C\f$.
 * The semantics is the set \f$D \times C\f$.
 * For general sets of states (not the cross product of discrete and continuous parts)
 * use symbolic_state_collection.
 */
class symbolic_state : public printable
{
public:
	typedef boost::shared_ptr<symbolic_state> ptr;
	typedef boost::shared_ptr<const symbolic_state> const_ptr;

	/**
	 * \brief Constructor.
	 */
	symbolic_state();
	/**
	 * \brief Constructor (adopting, doesn't clone ds and cs).
	 */
	symbolic_state(const discrete::discrete_set_ptr& ds,const continuous::continuous_set_ptr& cs);

	/** \brief Constructor that clones ds and cs.
	 */
	void clone(const discrete::discrete_set_ptr& ds,const continuous::continuous_set_ptr& cs);

	/** \brief Constructor that clones ds and cs.
	 */
	symbolic_state::ptr clone();

	/**
	 * \brief Returns true iff the symbolic state is empty.
	 *
	 * A symbolic state is empty
	 * if both the continuous and discrete set associated with it are empty.
	 *
	 * \return \p true if and only if both the continuous and discrete sets are empty.
	 */
	virtual bool is_empty() const;

	/**
	 * \brief Member field mutator function.
	 *
	 * Sets the continuous set member of the \p *this to the object pointed by passed parameter \p cs
	 * \param cs continuous_set_ptr pointing to a continous set implementation object which is to
	 * be set as \p *this continuous set member field.
	 */
	virtual void set_continuous_set(const continuous::continuous_set_ptr& cs);

	/**
	 * \brief Member field mutator function.
	 *
	 * Sets the discrete set member of \p *this to the discrete set pointed by the parameter \p ds.
	 * \param ds discrete_set::ptr pointing to a discrete set implementaion object which is to be set
	 * as \p *this discrete set member field.
	 */
	virtual void set_discrete_set(const discrete::discrete_set_ptr& ds);

	/**
	 * \brief Member field accessor function.
	 *
	 * \return The continuous set member of \p *this
	 */
	virtual const continuous::continuous_set_ptr& get_continuous_set() const;

	/**
	 * \brief Member field accessor function.
	 *
	 * \return The discrete set field of \p *this.
	 */
	virtual const discrete::discrete_set_ptr& get_discrete_set() const;

	/**
	 * Computes the intersection set between this and sstate efficiently ans
	 * assigns the result to \code this \endcode.
	 *
	 * \param sstate symbolic_state::ptr with which intersection has to be computed.
	 */
	virtual void intersection_assign(const symbolic_state::ptr& sstate);

	/**
	 * Output as a stream of characters.
	 */
	void print(std::ostream& os) const;

	/**
	 * \brief Virtual Destructor.
	 * C++ specific function.
	 */
	virtual ~symbolic_state();

	class output_format {
	public:
		output_format() { // default values
			preamble = "";
			element_separator = " & (";
			epilogue = ")";
			skip_continuous = false;
			skip_discrete = false;
		}
		;
		static output_format lf_separated() {
			output_format f;
			f.preamble = "";
			f.epilogue = "";
			f.element_separator = "\n";
			return f;
		}
		;
		static output_format lf_separated_continuous_only() {
			output_format f = lf_separated();
			f.skip_discrete = true;
			return f;
		}
		;
		static output_format matlab() {
			output_format f = lf_separated_continuous_only();
			// @todo
			return f;
		}
		;
		static output_format JVX() {
			output_format f = lf_separated_continuous_only();
			return f;
		}
		;
		std::string preamble; // written before matrix
		std::string epilogue; // written after matrix
		std::string element_separator; // written between elements
		bool skip_continuous; // don't output continuous
		bool skip_discrete; // don't output discrete
		variable_id_set output_variables; // if empty output all
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
	discrete::discrete_set_ptr d_set;
	continuous::continuous_set_ptr c_set;
	static output_format my_output_format;
};

}

#endif /*SYMBOLIC_STATE_H_*/
