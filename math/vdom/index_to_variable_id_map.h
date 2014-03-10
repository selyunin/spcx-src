#ifndef INDEX_TO_VARIABLE_ID_MAP_H_
#define INDEX_TO_VARIABLE_ID_MAP_H_

#include <vector>
#include <map>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "math/vdom/index_to_index_bimap.h"
#include "math/vdom/variable.h"

/** Forward declarations */
class index_to_variable_id_map;
typedef boost::shared_ptr<const index_to_variable_id_map> index_to_variable_id_map_ptr;
class positional_vdomain;

typedef unsigned int index_type;

/** A bidirectional map from \p variable_id to \p index_type.
 * It assigns indices to \p variable_id based on the assumptions
 * that the indices always start with 0, and that indices are
 * assigned without "holes" (if i1,i2 are indices in the map
 * then all i with i1 < i < i2 are indices in the map).
 *
 * @todo What does this mean?
 */
class index_to_variable_id_map : public boost::enable_shared_from_this<index_to_variable_id_map> {
	friend void get_common_map(const index_to_variable_id_map_ptr& M1,
			const index_to_variable_id_map_ptr& M2, index_to_variable_id_map_ptr& M,
			index_type& newdim, index_to_index_bimap& iimap2);
public:
	typedef unsigned int size_type;

	index_to_variable_id_map() {
	}
	;

	/** Returns the number of ids in the map.
	 */
	size_type dimensions() const;

	/** Returns the id belonging to the index \p i.
	 * Throws if the index is not in the map.
	 */
	const variable_id& get_id(const index_type& i) const;

	/** Returns the index belonging to the id \p id.
	 * Throws if the id is not in the map.
	 */
	const index_type& get_index(const variable_id& id) const;

	/** Returns true if id is in the map and false otherwise. */
	bool has_id(const variable_id& id) const;

	/** Returns the index belonging to the id \p id if there is one.
	 * has_id is true if the id was found, and false otherwise.
	 * If has_id is false, then the returned index_type defaults to 0.
	 */
	index_type check_for_index(const variable_id& id, bool& has_id) const;

	/** Returns all variable ids in the map. */
	const variable_id_set& get_ids() const;

	/** Returns the vector of variables in the map, ordered by
	 * their index.
	 */
	const std::vector<variable_id>& get_id_vector() const;

	/** Returns a pointer to a new map constructed as follows:
	 * Add every \p variable_id in \p id_set to the map,
	 * assigning the next higher index in sequence.
	 * If a \p variable_id already exists, it keeps its old index.
	 */
	static index_to_variable_id_map_ptr get_map_with_ids(const variable_id_set& id_set);

	/** Returns a pointer to a new map constructed as follows:
	 * Add \p variable_id in \p id_set to the map,
	 * assigning the next higher index in sequence.
	 * If \p variable_id already exists, it keeps its old index.
	 */
	index_to_variable_id_map_ptr get_map_with_id_added(const variable_id& id) const;

	/** Returns a pointer to a new map constructed as follows:
	 * Add every \p variable_id in \p id_set to the map,
	 * assigning the next higher index in sequence.
	 * If a \p variable_id already exists, it keeps its old index.
	 */
	index_to_variable_id_map_ptr get_map_with_ids_added(const variable_id_set& id_set) const;

	/** Returns a pointer to a new map constructed as follows:
	 * Remove every \p variable_id in \p id_set from the map,
	 * filling the created "holes" between indices by moving the higher indices down.
	 * I.e., for any index i the id at i+1 is the id that was at the smallest index j>i such that
	 * j is not in \p id_set.
	 *
	 * \todo{The current implementation is somewhat inefficient, this should be rewritten if it turns out
	 * to be slowing down things in practice.}
	 */
	index_to_variable_id_map_ptr get_map_with_ids_removed(const variable_id_set& id_set) const;

	/** Returns a pointer to a new map in which the primedness of each variable with primedness d has been set to p. */
	index_to_variable_id_map_ptr get_map_with_primedness_reassigned(unsigned int d, unsigned int p = 0) const;

	/** Returns a pointer to a new map in which the primedness of each variable with primedness d has been increased by 1.
	 * If d is 0, increase all. */
	index_to_variable_id_map_ptr get_map_with_primedness_increased(unsigned int d = 0) const;

	/** Returns a pointer to a new map in which the primedness of each variable with primedness d  has been decreased by 1.
	 * If d is 0, decrease all. */
	index_to_variable_id_map_ptr get_map_with_primedness_decreased(unsigned int d = 0) const;

	/** Returns true iff all variables in iimap are in *this. */
	bool contains_variables(const index_to_variable_id_map_ptr& iimap) const;

	/** Output as a chain of characters. */
	void print(std::ostream& os) const;

	/** Get the \p index_to_variable_id_map_ptr of the empty map (which has no variables). */
	static const index_to_variable_id_map_ptr& empty_map();

private:
	/** Create a new index_to_variable_id_map from id_vector and id_to_index_map. */
	index_to_variable_id_map(std::vector<variable_id> id_vector,
			std::map<variable_id, index_type> id_to_index_map, variable_id_set ids) :
		my_id_vector(id_vector), my_id_to_index_map(id_to_index_map), my_ids(ids) {
	}
	;

	/** Obtain an index_to_variable_id_map_ptr corresponding to the id_vector \p new_v.
	 * First, the cache is searched for the corresponding pointer. If none is found,
	 * a new map is instantiated (involving the computation of the corresponding \p id_to_index_map),
	 * added to the cache, and its pointer returned. */
	static index_to_variable_id_map_ptr get_index_to_variable_id_map_ptr(const std::vector<
			variable_id>& new_v);

	static variable_id_set compute_id_set(const std::vector<variable_id>& id_vec);

	// Properties
	std::vector<variable_id> my_id_vector;
	std::map<variable_id, index_type> my_id_to_index_map;
	variable_id_set my_ids;
	static std::map<std::vector<variable_id>, index_to_variable_id_map_ptr>
			index_to_variable_id_map_cache;

	friend class positional_vdomain;
};

/**
 * Output as a stream of characters. Calls the print method.
 */
std::ostream& operator<<(std::ostream& os, const index_to_variable_id_map_ptr dsp);

/** Constructs a common index_to_variable_id_map from two maps M1 and M2.
 * The resulting map consists of M1 plus a map for the variables of M2
 * that are not in M1, in the same order as in M2 but in consecutive
 * indices.
 * Let M1 = 0->i_0,1->i_1,...,m-1->i_{m-1},
 * M2 = 0->j_0,1->j_1,...n-1->j_{n-1},
 * and let J_1,...,J_z be the variable_ids of M2 that are not in M1.
 * Then the result is
 * M = 0->i_0,1->i_1,...,m-1->i_{m-1},m->J_1,m+z-1->J_z.
 * The size of the new index domain (m+z) is returned in newdim.
 * It also produces an index_to_index_bimap iimap2, attributing the
 * indices of M2 to their new place in M.
 * Note: The domain of iimap2 is over indices 0,...,newdim-1, because the
 * function applying the mapping (e.g. in the PPL)
 * may not be capable of augmenting the dimension. Thus we assume that
 * the space of the set corresponding to M2 is augmented (by embedding) to
 * dimension newdim, and then iimap2 is applied to remap the indices.
 * In iimap2 the indices n,...,newdim-1 are mapped to arbitrary places,
 * since they are assumed to be a result of the embedding, and therefore
 * irrelevant.
 */
void get_common_map(const index_to_variable_id_map_ptr& M1, const index_to_variable_id_map_ptr& M2,
		index_to_variable_id_map_ptr& M, index_type& newdim, index_to_index_bimap& iimap2);

/** Returns the map necessary to remap dimensions from a space from_M to a space
 * to_M.
 *
 * Space dimensions are considered to be removed if they do not exist in to_M. */
index_to_index_bimap get_index_to_index_mapping(const index_to_variable_id_map_ptr& to_M, const index_to_variable_id_map_ptr& from_M);

#include "index_to_variable_id_map.hpp"

#endif /*INDEX_TO_VARIABLE_ID_MAP_H_*/
