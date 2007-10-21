/* Copyright (c) 1997-2006
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

#ident "$Project: polymake $$Id: gmp_init.cc 7315 2006-04-02 21:37:53Z gawrilow $"

#include <memory>
#include "gmp_init.h"

#if defined(__GNUC__)
#if __GNUC__==3 && __GNUC_MINOR__==3
namespace {
   void* pm_gmp_allocate(size_t n)
   {
      return __builtin_expect(n==0,0) ? 0 : std::__alloc::allocate(n);
   }
}
# define pm_gmp_deallocate std::__alloc::deallocate
# define pm_gmp_reallocate std::__alloc::reallocate
#endif // gcc 3.3

#if __GNUC__==3 && __GNUC_MINOR__==4 || __GNUC__==4
# include <ext/pool_allocator.h>

namespace {
# if __GNUC__==3 && __GNUC_PATCHLEVEL__ == 0
   typedef __gnu_cxx::__pool_alloc<true,0> raw_alloc;

   void* pm_gmp_allocate(size_t n)
   {
      return __builtin_expect(n==0,0) ? 0 : raw_alloc::allocate(n);
   }

   void pm_gmp_deallocate(void* p, size_t n) { raw_alloc::deallocate(p,n); }
# else  // >3.4.0
   __gnu_cxx::__pool_alloc<char> gmp_allocator;

   void* pm_gmp_allocate(size_t n) { return gmp_allocator.allocate(n); }
   void pm_gmp_deallocate(void* p, size_t n) { gmp_allocator.deallocate(reinterpret_cast<char*>(p), n); }
# endif

   void* pm_gmp_reallocate(void* p, size_t old_sz, size_t new_sz)
   {
      static const bool use_new=getenv("GLIBCPP_FORCE_NEW") || getenv("GLIBCXX_FORCE_NEW");
      const size_t align=8, limit=128;
      if (!use_new && ((old_sz+align-1)&~(align-1)) == ((new_sz+align-1)&~(align-1)) && new_sz < limit)
	 return p;
      void* new_p=pm_gmp_allocate(new_sz);
      if (new_p) {
	 memcpy(new_p, p, old_sz < new_sz ? old_sz : new_sz);
	 pm_gmp_deallocate(p,old_sz);
      }
      return new_p;
   }
}
#endif // gcc >=3.4
#elif defined(__INTEL_COMPILER)
namespace {
   std::allocator<char> gmp_allocator;
   void* pm_gmp_allocate(size_t n) { return gmp_allocator.allocate(n); }
   void pm_gmp_deallocate(void* p, size_t n) { gmp_allocator.deallocate(reinterpret_cast<char*>(p), n); }
   void* pm_gmp_reallocate(void* p, size_t old_sz, size_t new_sz)
   {
      void* new_p=pm_gmp_allocate(new_sz);
      if (new_p) {
	 memcpy(new_p, p, old_sz < new_sz ? old_sz : new_sz);
	 pm_gmp_deallocate(p,old_sz);
      }
      return new_p;
   }
}
#endif // __INTEL_COMPILER

gmp_alloc_init::gmp_alloc_init()
{
   static bool done=false;
   if (!done) {
      mp_set_memory_functions(pm_gmp_allocate, pm_gmp_reallocate, pm_gmp_deallocate);
      done=true;
   }
}
