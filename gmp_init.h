/* Copyright (c) 1997-2002
   Technische Universitaet Berlin, Germany
   Department of Mathematics,
   Research Group Algorithmic and Discrete Mathematics
   mailto:polymake@math.tu-berlin.de

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, 59 Temple Place - Suite 330 Boston, MA 02111-1307, USA.
*/

#ifndef _POLYMAKE_GMP_INIT_H
#define _POLYMAKE_GMP_INIT_H "$Project: polymake $$Id: gmp_init.h,v 1.8 2003/01/06 17:07:06 gawrilow Exp $"

#ifdef __GNUC__
#if __GNUC_MINOR__<3 && !defined(__APPLE__)
#pragma interface
#endif
#endif

#include <gmp.h>
#include <tweaks.h>

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
