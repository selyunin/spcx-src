/*
 * glpk_cleaner.h
 *
 *  Created on: Mar 2, 2011
 *      Author: frehse
 */

#ifndef GLPK_CLEANER_H_
#define GLPK_CLEANER_H_

#include <glpk.h>

/** Factory for glp_prob objects.
 *
 * This factory should be used to create and destroy all glp_prob objects.
 * It cleans up glpk memory upon program exit
 *
 * See http://old.nabble.com/Memory-leak--td28048181.html:
 > As it can be seen, there seems to remain 1.108bytes still reachable,
 >  maybe due to a missing free call. I guess that the environment block
 > env (glplib01.c:67) has not been deallocated, that is,
 > lib_free_env(void) may have not been called properly.

 Yes, those 1108 bytes are occupied by the glpk environment block. It
 is allocated implicitly on the very first call to any glpk api routine.
 To deallocate it you need to call glp_free_env. Note, however, that
 glp_free_env frees all the memory allocated by glpk and, thus,
 invalidates all problem objects which are still exist. For more details
 please see the glpk reference manual, Section "GLPK environment
 routines".
 *
 */
class glpk_cleaner {
public:
	static glp_prob* create_glp_prob() {
		++prb_count;
		return glp_create_prob();
	}
	static void delete_glp_prob(glp_prob*& p) {
		--prb_count;
		glp_delete_prob(p);
		if (prb_count == 0) {
//			glp_free_env();
		}
	}
private:
	static unsigned int prb_count;
};


#endif /* GLPK_CLEANER_H_ */
