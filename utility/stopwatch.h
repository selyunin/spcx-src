#ifndef GUARD_stopwatch_h
#define GUARD_stopwatch_h

/***************************************************************************
 stopwatch.h  -  description
 -------------------
 begin                : Thu Feb 5 2004
 copyright            : (C) 2004 by Goran Frehse
 email                : goran.frehse@gmx.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <string>
#include <boost/date_time/posix_time/posix_time.hpp> //include all types plus i/o

class stopwatch {
public:
	typedef boost::posix_time::ptime time_type;
	typedef boost::posix_time::time_duration duration_type;

	stopwatch();

	~stopwatch();

	/** Returns the time since construction in microsecond accuracy. */
	duration_type value() const;

	/** Returns the value as a string in HH:MM:SS.fff. */
	std::string string_value() const;

	/** Returns the time since last report or call to delta. */
	duration_type delta();

	/** Returns the current universal time. */
	static time_type universal_time();

	/** Convert to seconds at millisecond accuracy. */
	static double seconds(const duration_type& t);

private:
	time_type start;
	time_type report_time;
};

#endif
