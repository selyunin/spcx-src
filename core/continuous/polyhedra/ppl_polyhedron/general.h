#ifndef GUARD_general_h
#define GUARD_general_h
/***************************************************************************
                          general.h  -  description
                             -------------------
    begin                : Wed Feb 11 2004
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

//#include <stdio.h>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

namespace ppl_polyhedron {

//static unsigned int VERBOSE_LEVEL;
//static unsigned int line_number;

//using namespace std;

void pause_key();

void throw_error(std::string s);

void throw_warning(std::string s);

void progress_dot();
void progress_dot_end();

void message(unsigned int level, std::string s);
std::string message_prefix(unsigned int l);

// greatest common divisor for cost_fun_type calculations
// use gcd_assign instead
int gcd( int num1, int num2 );


// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------

// ATTENTION: the following version seems to return a different iterator, so that the comparison
//            find(x,y)!=x.end() can always be true!
//template <typename key_class, typename value_class>
//typename map <key_class,value_class>::const_iterator find(map < key_class , value_class> obj, value_class val)
////void find(map < key_class , value_class> obj, value_class val)
//{
//  typename map <key_class,value_class>::const_iterator i=obj.end();
//  for (i=obj.begin(); i!=obj.end(); ++i)
//    if (i->second==val)
//      return i;
//  return i;
//};

//template <typename key_class, typename value_class>
////typename map <key_class,value_class>::iterator map_val_find(typename map <key_class,value_class>::iterator i, typename map <key_class,value_class>::iterator iend, value_class val)
//typename map <key_class,value_class>::const_iterator map_val_find(typename map <key_class,value_class>::const_iterator i, typename map <key_class,value_class>::const_iterator iend, value_class val)
////void find(map < key_class , value_class> obj, value_class val)
//{
//  while (i!=iend)
//  {
//    if (i->second==val)
//      return i;
//    ++i;
//  };
//  return iend;
//};

template <typename key_class, typename value_class>
typename std::map <key_class,value_class>::const_iterator find(std::map < key_class , value_class> obj, typename std::map <key_class,value_class>::const_iterator i, typename std::map <key_class,value_class>::const_iterator iend, value_class val)
//void find(map < key_class , value_class> obj, value_class val)
{
  while (i!=iend)
  {
    if (i->second==val)
      return i;
    ++i;
  };
  return iend;
}

}

#endif
