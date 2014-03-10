/***************************************************************************
                          general.cpp  -  description
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

#include "stdlib.h" /* for abort() -wsc */
#include "core/continuous/polyhedra/ppl_polyhedron/general.h"

namespace ppl_polyhedron {

unsigned int VERBOSE_LEVEL=600011;
unsigned int line_number=0;

using namespace std;

// --------------------------------------------------------------

void pause_key()
{
  cout << "waiting for key:" << endl;
  char c;
  cin >> c;
}

void throw_error(string s)
{
  cout << "ERROR: " << s << endl << flush;
  cout << "       (while parsing line " << line_number << ")." << endl << flush;
  //cout << "ERROR: " << s << flush << endl;
  //cerr << "ERROR: " << s << flush << endl;
  abort();
//  pause_key();
}

void throw_warning(string s)
{
  cout << "WARNING: " << s << flush << endl;
  cout << "Press a key and then press Enter to continue:"<< flush << endl;
  pause_key();
}

unsigned int PROGRESS_DOT_COUNT = 0;

void progress_dot()
{
  if ((VERBOSE_LEVEL & 10) !=0)
  {
    if (PROGRESS_DOT_COUNT==0)
      cout << message_prefix(VERBOSE_LEVEL);
    ++PROGRESS_DOT_COUNT;
    cout << "." << flush;
    if (PROGRESS_DOT_COUNT>=50)
    {
      PROGRESS_DOT_COUNT=0;
      cout << endl;
    };
  };
}

void progress_dot_end()
{
  if ((VERBOSE_LEVEL & 10) !=0)
    cout << endl;
  PROGRESS_DOT_COUNT=0;
}

void message(unsigned int level, string s)
{
  if (VERBOSE_LEVEL >= level)
  {
    if (level<4000)
      cout << endl;

    string prae=message_prefix(level);
    cout << prae << s << flush << endl;
    if (level<4000)
    {
    cout << prae;
      for (unsigned int i=0;i<s.size();++i)
        cout << "-";
    cout << endl << flush;
    };
  };
}

string message_prefix(unsigned int l)
{
  string prae="";
  l=l/1000;
  while (l>1)
  {
    l=l/2;
    prae+="  ";
  };
  return prae;
}


// greatest common divisor for cost_fun_type calculations
// use gcd_assign instead
int gcd( int num1, int num2 )
{
  int remainder = num2 % num1;
  if ( remainder != 0 )
    return gcd( remainder,num1 );
  return num1;
}

}

