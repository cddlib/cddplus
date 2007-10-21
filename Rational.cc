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

#ident "$Project: polymake $$Id: Rational.cc 7565 2007-01-16 16:29:43Z gawrilow $"

#include <cctype>
#include "Rational.h"

Rational& Rational::set(const char* s) throw (gmp_error)
{
   const char* digit=s;
   while (*digit && *digit!='/') ++digit;
   int numerator_digits=digit-s;
   if (!numerator_digits)
      throw gmp_error("Rational: syntax error in string");

   Integer::little_buffer num(numerator_digits+1);
   memcpy(num,s,numerator_digits);
   num[numerator_digits]=0;
   if (mpz_set_str(mpq_numref(rep), num, 0) < 0)
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

std::string Rational::to_string(int base) const
{
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

// defined in GMP, but not public
extern "C" {
   void __gmp_divide_by_zero(void) __attribute__((noreturn));
}

void Rational::zero_division()
{
   __gmp_divide_by_zero();
}
