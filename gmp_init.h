/* Copyright (c) 1997-2007
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
   http://www.math.tu-berlin.de/polymake,  mailto:polymake@math.tu-berlin.de

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*/

#ifndef _POLYMAKE_GMP_INIT_H
#define _POLYMAKE_GMP_INIT_H "$Project: polymake $$Id: gmp_init.h 7556 2007-01-12 17:36:36Z gawrilow $"

//#include <ext/defines.h>
#include <gmp.h>

struct gmp_alloc_init {
   gmp_alloc_init();
};

namespace {
   gmp_alloc_init do_init;
}

namespace std_ext {
   template <typename _Key> struct hash;

   template <> struct hash<MP_INT> {
   protected:
      size_t _do(mpz_srcptr a) const {
	 size_t result=0;
	 for (int i=0, n=mpz_size(a); i<n; ++i)
	    (result <<= 1) ^= mpz_getlimbn(a, i);
	 return result;
      }
   public:
      size_t operator() (const MP_INT& a) const {
	 return _do(&a);
      }
   };
}

#endif // _POLYMAKE_GMP_INIT_H

// Local Variables:
// mode:C++
// End:
