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

#ident "$Project: polymake $$Id: Rational.cc,v 1.15 2002/12/11 14:09:07 gawrilow Exp $"

#if defined(__GNUC__) && __GNUC_MINOR__<3 && !defined(__APPLE__)
#pragma implementation
#endif

#include <cctype>
#include "Rational.h"
#ifdef __APPLE__
# include <cstdlib>
#else
# include <alloca.h>
#endif

Rational& Rational::set(const char* s) throw (gmp_error) {
   const char* digit=s;
   while (*digit && *digit!='/') ++digit;
   int numerator_digits=digit-s;
   if (!numerator_digits)
      throw gmp_error("Rational: syntax error in string");

   char *numerator=(char*)alloca(numerator_digits+1);
   memcpy(numerator,s,numerator_digits);
   numerator[numerator_digits]=0;
   if (mpz_set_str(mpq_numref(rep), numerator, 0) < 0)
      throw gmp_error("Rational: syntax error in numerator");

   if (*digit) {
      ++digit;
      if (!isdigit(*digit) || mpz_set_str(mpq_denref(rep), digit, 0) < 0)
	 throw gmp_error("Rational: syntax error in denominator");
      canonicalize();
   } else {
      mpz_set_ui(mpq_denref(rep),1);
   }
   return *this;
}

std::istream& operator>> (std::istream& is, Rational& a) {
   Integer::read(is,mpq_numref(a.rep));
   if (is.peek() == '/') {
      is.ignore();
      Integer::read(is,mpq_denref(a.rep),false);
      a.canonicalize();
   } else {
      mpz_set_ui(mpq_denref(a.rep), 1);
   }
   return is;
}

std::ostream& operator<< (std::ostream &os, const Rational& a) {
   if (!mpz_cmp_ui(mpq_denref(a.rep),1))
      return os << numerator(a);
   size_t sn=Integer::strsize(os,mpq_numref(a.rep)),
          sd=Integer::strsize(os,mpq_denref(a.rep));
   char* buf=(char*)alloca(sn+sd);
   Integer::putstr(os,buf,mpq_numref(a.rep));
   char *end=(char*)memchr(buf,0,sn);
   *end++='/';
   Integer::putstr(os,end,mpq_denref(a.rep));
   return os << buf;
}

std::string Rational::to_string(int base) const {
   if (!mpz_cmp_ui(mpq_denref(rep),1))
      return numerator(*this).to_string();
   std::string s(mpz_sizeinbase(mpq_numref(rep), base)+2 +	// numerator with possible sign and slash
		 mpz_sizeinbase(mpq_denref(rep), base)+1,	// denominator and terminating '\0'
		 '\0');
   char *buf=const_cast<char*>(s.data());
   mpz_get_str(buf,base,mpq_numref(rep));
   buf+=strlen(buf);
   *buf++='/';
   mpz_get_str(buf,base,mpq_denref(rep));
   s.resize(s.find('\0'));
   return s;
}
