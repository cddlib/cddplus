#ifndef _POLYMAKE_TWEAKS_H
#define _POLYMAKE_TWEAKS_H "$Project: polymake $$Id: tweaks.h,v 1.4 2003/01/06 17:07:06 gawrilow Exp $"

#ifdef std_ext
#  error std_ext already defined!
#endif

#if defined(__GNUC__) && __GNUC__ == 3 && __GNUC_MINOR__ > 0
#  define std_ext __gnu_cxx
#  define STD_EXT_HEADER(name) <ext/name>
#else
# if !defined(__INTEL_COMPILER)
#  define std_ext std
# endif
#  define STD_EXT_HEADER(name) <name>
#endif

#if defined(__GNUC__) || defined(__COMO__)
#  define _POLYMAKE_ALIGN(what,n) what __attribute__ ((aligned (n)))
#elif defined(__INTEL_COMPILER)
# if __INTEL_COMPILER<=700
#  define _POLYMAKE_ALIGN(what,n) __declspec(align(16)) what
# else
#  define _POLYMAKE_ALIGN(what,n) __declspec(align(n)) what
# endif
#endif

#endif // _POLYMAKE_TWEAKS_H

// Local Variables:
// mode:C++
// End:
