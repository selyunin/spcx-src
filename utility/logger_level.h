/*
 * logger_level.h
 *
 *  Created on: Oct 16, 2010
 *      Author: frehse
 */

#ifndef LOGGER_LEVEL_H_
#define LOGGER_LEVEL_H_

struct logger_level {
	/** These are the logging levels.
	 *
	 * @attention OFF is reserved for setting the active
	 * level of the logger so that no messages are logged.
	 * ALWAYS is reserved for messages that should be output
	 * unless the active level is OFF. The active level
	 * can not be ALWAYS.
	 * NEVER is reserved for messages that should never
	 * be output. The active level can not be NEVER.
	 *
	 * DEBUG is intended for profiling purposes, so avoid
	 * costly output messages.
	 * Debug levels DEBUG1,...,DEBUG4 are intended for
	 * debug output of increasing detail (and irrelevance).
	 */
	typedef enum {
		OFF, ALWAYS, LOW, MEDIUM, HIGH, DEBUG, DEBUG1, DEBUG2, DEBUG3, DEBUG4, DEBUG5, DEBUG6, DEBUG7, NEVER
	} level;
};

#endif /* LOGGER_LEVEL_H_ */
