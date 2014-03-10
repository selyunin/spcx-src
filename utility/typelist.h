#ifndef TYPELIST_H_
#define TYPELIST_H_

/** Typelist Template
 * Taken from Andrei Alexandrescu. Generic Programming: Typelists and Applications 
 * http://www.ddj.com/cpp/184403813 */

template <class H, class T>
struct typelist
{
    typedef H head;
    typedef T tail;
};

class null_typelist {};

#endif /*TYPELIST_H_*/
