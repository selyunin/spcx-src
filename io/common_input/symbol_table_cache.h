/*
 * symbol_table_cache.h
 *
 *  Created on: Mar 22, 2011
 *      Author: frehse
 */

#ifndef SYMBOL_TABLE_CACHE_H_
#define SYMBOL_TABLE_CACHE_H_

#include "symbol_table.h"

namespace parser {

/** A class storing symbol tables that are associated to components. */
class symbol_table_cache {
public:
	static symbol_table& get_symbol_table(const std::string& component) {
		return my_tables[component];
	}

	static bool has_symbol_table(const std::string& component) {
		container_type::iterator it = my_tables.find(component);
		return (it != my_tables.end());
	}

private:
	typedef std::map<std::string, symbol_table> container_type;
	static container_type my_tables;
};

}

#endif /* SYMBOL_TABLE_CACHE_H_ */
