#ifndef __TIMINGS_H
#define __TIMINGS_H

#ifdef TIMINGS
#include <sys/time.h>
#endif

#ifdef TIMINGS
#ifdef TIMINGS_REALTIME
#define TIMING_DECLARE() struct timeval t1, t2;
#define TIMING_BEGINS() gettimeofday(&t1, 0)
#define TIMING_ENDS(x) do {\
    gettimeofday(&t2, 0); \
    x += (t2.tv_usec-t1.tv_usec)*1E-6 + (t2.tv_sec-t1.tv_sec); \
    } while(0)
#else
#include <sys/resource.h>

#define TIMING_DECLARE() struct rusage t1, t2;
#define TIMING_BEGINS() getrusage(RUSAGE_SELF, &t1)
#define TIMING_ENDS(x) do {\
    getrusage(RUSAGE_SELF, &t2); \
    x += (t2.ru_utime.tv_usec-t1.ru_utime.tv_usec)*1E-6 + (t2.ru_utime.tv_sec-t1.ru_utime.tv_sec); \
    } while(0)
#endif
#else

#define TIMING_DECLARE()
#define TIMING_BEGINS()
#define TIMING_ENDS(x)
#endif

#endif
