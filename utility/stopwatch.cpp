/***************************************************************************
 stopwatch.cpp  -  description
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

#include "stopwatch.h"

using namespace std;

stopwatch::stopwatch() {
	start = universal_time();
	report_time = universal_time();
}

stopwatch::~stopwatch() {
}

stopwatch::duration_type stopwatch::value() const {
	// returns the time since the start
	stopwatch::duration_type dur = universal_time() - start;
	return dur;
}

std::string stopwatch::string_value() const {
	return to_simple_string(value()).substr(0, 12);
}

stopwatch::duration_type stopwatch::delta() {
	// returns the time since the last report or delta call
	stopwatch::duration_type dur = universal_time() - report_time;
	report_time = universal_time();
	return dur;
}

stopwatch::time_type stopwatch::universal_time() {
	return boost::posix_time::microsec_clock::universal_time();
}

double stopwatch::seconds(const duration_type& t) {
	return double(t.total_milliseconds()) / 1000.0;
}
