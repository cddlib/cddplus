/* Copyright (c) 1998
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

#ident "$Project: polymake $$Id: gmp_rational.cc,v 1.6 1999/02/08 14:20:23 gawrilow Exp $"

#pragma implementation

#include "gmp_rational.h"
#include <cmath>

Rational& Rational::operator= (double d)
{
   double x=d;
   MP_INT *num=mpq_numref(rep),
          *den=mpq_denref(rep);
   mpz_set_ui(num, 0);
   mpz_set_ui(den, 1);
   if (x != 0.0) {
      int neg = x < 0;
      if (neg)
	 x = -x;

      const long shift = 15;         // a safe shift per step
      const double width = 32768.0;  // = 2^shift
      const int maxiter = 20;        // ought not be necessary, but just in case,
                                     // max 300 bits of precision
      int expt;
      double mantissa = frexp(x, &expt);
      long exponent = expt;
      double intpart;
      int k = 0;
      while (mantissa != 0.0 && k++ < maxiter) {
	  mantissa *= width;
	  mantissa = modf(mantissa, &intpart);
	  mpz_mul_2exp(num, num, shift);         // num <<= shift;
	  mpz_add_ui(num, num, static_cast<long>(intpart));   // num += (long)intpart;
	  exponent -= shift;
      }
      if (exponent > 0)
	 mpz_mul_2exp(num, num, exponent); //num <<= exponent;
      else if (exponent < 0)
	 mpz_mul_2exp(den, den, -exponent); //den <<= -exponent;
      if (neg)
	 mpz_neg(num, num);
   }
   mpq_canonicalize(rep);
   return *this;
}

istream& operator>> (istream& is, Rational& a)
{
   read(is,mpq_numref(a.rep),'/');
   if (is.peek() == '/') {
      is.ignore();
      read(is,mpq_denref(a.rep),0,false);
      mpq_canonicalize(a.rep);
   } else {
      mpz_set_ui(mpq_denref(a.rep), 1);
   }
   return is;
}

ostream& operator<< (ostream &os, const Rational& a)
{
   static const Integer One(1);
   os << a.numerator();
   if (a.denominator()!=One) os << "/" << a.denominator();
   return os;
}
