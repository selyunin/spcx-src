#ifndef GUARD_bidirectional_map_h
#define GUARD_bidirectional_map_h
/***************************************************************************
 *   Copyright (C) 2005 by Goran Frehse   *
 *   goran.frehse@imag.fr   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/* Boost.MultiIndex example of a bidirectional map.
 *
 * Copyright 2003-2005 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#if !defined(NDEBUG)
#define BOOST_MULTI_INDEX_ENABLE_INVARIANT_CHECKING
#define BOOST_MULTI_INDEX_ENABLE_SAFE_MODE
#endif

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <iostream>
#include <string>
#include <utility>

using boost::multi_index_container;
using namespace boost::multi_index;

/* tags for accessing both sides of a bidirectional map */

struct from{};
struct to{};

/* The class template bidirectional_map wraps the specification
 * of a bidirectional map based on multi_index_container.
 */

template<typename FromType,typename ToType>
struct bidirectional_map
{
  typedef std::pair<FromType,ToType> value_type;

#if defined(BOOST_NO_POINTER_TO_MEMBER_TEMPLATE_PARAMETERS) ||\
    defined(BOOST_MSVC)&&(BOOST_MSVC<1300) ||\
    defined(BOOST_INTEL_CXX_VERSION)&&defined(_MSC_VER)&&\
           (BOOST_INTEL_CXX_VERSION<=700)

/* see Compiler specifics: Use of member_offset for info on member<> and
 * member_offset<>
 */

  BOOST_STATIC_CONSTANT(unsigned,from_offset=offsetof(value_type,first));
  BOOST_STATIC_CONSTANT(unsigned,to_offset  =offsetof(value_type,second));

  typedef multi_index_container<
    value_type,
    indexed_by<
      ordered_unique<
        tag<from>,member_offset<value_type,FromType,from_offset> >,
      ordered_unique<
        tag<to>,  member_offset<value_type,ToType,to_offset> >
    >
  > type;

#else

  /* A bidirectional map can be simulated as a multi_index_container
   * of pairs of (FromType,ToType) with two unique indices, one
   * for each member of the pair.
   */

  typedef multi_index_container<
    value_type,
    indexed_by<
      ordered_unique<
        tag<from>,member<value_type,FromType,&value_type::first> >,
      ordered_unique<
        tag<to>,  member<value_type,ToType,&value_type::second> >
    >
  > type;

#endif
};

#endif
