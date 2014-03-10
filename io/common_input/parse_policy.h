/*
 * parse_policy.h
 *
 *  Created on: Mar 21, 2011
 *      Author: frehse
 */

#ifndef PARSE_POLICY_H_
#define PARSE_POLICY_H_

namespace parser {

/** Properties of parsing */
class parse_policy {
public:
	/** Default parse policy
	 *
	 * Scalars have dimension 0 and
	 * unknown variables are added.
	 * */
	parse_policy() :
		scalar_dim(0), add_unknown_vars(true), use_symbol_table_cache(true) {
	}

	/** Dimension of scalars. If 1, scalars
	 * are interpreted as vectors of dimension 1.
	 */
	unsigned int scalar_dim;

	/** If true, add variables even they are
	 * not in the symbol table. If false, throw.
	 */
	bool add_unknown_vars;

	/** If true, add symbol tables to symbol_table_cache,
	 * so that they can be reused by other components.
	 */
	bool use_symbol_table_cache;

	/** SX format parse policy
	 *
	 * Scalars have dimension 1 and
	 * throw if variables are not in the
	 * symbol table.
	 */
	static parse_policy SX_policy() {
		parse_policy SXp;
		SXp.scalar_dim = 1;
		SXp.add_unknown_vars = false;
		return SXp;
	}
};

}

#endif /* PARSE_POLICY_H_ */
