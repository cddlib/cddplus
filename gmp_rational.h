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

#ifndef _POLYMAKE_GMP_RATIONAL_H
#define _POLYMAKE_GMP_RATIONAL_H

//! project = "$Project: polymake $"
//!   rcsid = "$Id: gmp_rational.h,v 1.1 1999/02/08 14:18:27 gawrilow Exp $"

#pragma interface

#include <gmp_integer.h>

/** Multiple Precision Rational class wrapping GMP rational type mpq_t
    Based on the GNU Multiple Precision Library (GMP)
 
    Introduction to GNU MP
    
    GNU MP is a portable library written in C for arbitrary precision
    arithmetic on integers, rational numbers, and floating-point numbers.
    It aims to provide the fastest possible arithmetic for all applications
    that need higher precision than is directly supported by the basic C
    types.
    
    Many applications use just a few hundred bits of precision; but some
    applications may need thousands or even millions of bits.  MP is
    designed to give good performance for both, by choosing algorithms
    based on the sizes of the operands, and by carefully keeping the
    overhead at a minimum.
    
    The speed of MP is achieved by using fullwords as the basic
    arithmetic type, by using sophisticated algorithms, by including
    carefully optimized assembly code for the most common inner loops for
    many different CPUs, and by a general emphasis on speed (as opposed to
    simplicity or elegance).
    
    There is carefully optimized assembly code for these CPUs: DEC
    Alpha, Amd 29000, HPPA 1.0 and 1.1, Intel Pentium and generic x86,
    Intel i960, Motorola MC68000, MC68020, MC88100, and MC88110,
    Motorola/IBM PowerPC, National NS32000, IBM POWER, MIPS R3000, R4000,
    SPARCv7, SuperSPARC, generic SPARCv8, and DEC VAX.  Some optimizations
    also for ARM, Clipper, IBM ROMP (RT), and Pyramid AP/XP.

    @index main
*/
class Rational {
private:
   /// GMP's representation
   mpq_t rep;
public:
   /// Initializes to 0.
   Rational()
   {
      mpq_init(rep);
   }

   /// Uses \\mpq_set\\.
   Rational(const Rational& r)
   {
      mpq_init(rep);
      mpq_set(rep, r.rep);
   }

   /// @group Given numerator, denominator set to 1.
   Rational(const Integer& num)
   {
      mpz_init_set(mpq_numref(rep), num.rep);
      mpz_init_set_ui(mpq_denref(rep), 1);
   }

   Rational(long num)
   {
      mpz_init_set_si(mpq_numref(rep), num);
      mpz_init_set_ui(mpq_denref(rep), 1);
   }

   Rational(int num)
   {
      mpz_init_set_si(mpq_numref(rep), num);
      mpz_init_set_ui(mpq_denref(rep), 1);
   }

   /// @group Given both numerator and denominator, of the same type
   Rational(const Integer& num, const Integer& den)
   {   
      mpz_init_set(mpq_numref(rep), num.rep);
      mpz_init_set(mpq_denref(rep), den.rep);
      mpq_canonicalize(rep);
   }

   Rational(long num, long den)
   {
      mpz_init_set_si(mpq_numref(rep), num);
      mpz_init_set_si(mpq_denref(rep), den);
      mpq_canonicalize(rep);
   }

   Rational(int num, int den)
   {
      mpz_init_set_si(mpq_numref(rep), num);
      mpz_init_set_si(mpq_denref(rep), den);
      mpq_canonicalize(rep);
   }

   /// @group Mixed argument types.
   Rational(const Integer& num, long den)
   {
      mpz_init_set(mpq_numref(rep), num.rep);
      mpz_init_set_si(mpq_denref(rep), den);
      mpq_canonicalize(rep);
   }

   Rational(long num, const Integer& den)
   {   
      mpz_init_set_si(mpq_numref(rep), num);
      mpz_init_set(mpq_denref(rep), den.rep);
      mpq_canonicalize(rep);
   }

   Rational(const Integer& num, int den)
   {
      mpz_init_set(mpq_numref(rep), num.rep);
      mpz_init_set_si(mpq_denref(rep), den);
      mpq_canonicalize(rep);
   }

   Rational(int num, const Integer& den)
   {   
      mpz_init_set_si(mpq_numref(rep), num);
      mpz_init_set(mpq_denref(rep), den.rep);
      mpq_canonicalize(rep);
   }

   Rational(long num, int den)
   {
      mpz_init_set_si(mpq_numref(rep), num);
      mpz_init_set_si(mpq_denref(rep), den);
      mpq_canonicalize(rep);
   }

   /// @endgroup
   Rational(int num, long den)
   {
      mpz_init_set_si(mpq_numref(rep), num);
      mpz_init_set_si(mpq_denref(rep), den);
      mpq_canonicalize(rep);
   }

   Rational(const double d)
   {
      mpq_init(rep);
      *this=d;
   }

   friend Integer& Integer::operator= (const Rational& r);

   ~Rational()
   {
      mpq_clear (rep);
   }

   /** @group Get the numerator rsp. denominator.
       A new \\Integer\\ object is returned, !!not!! the reference to the existing component.
   */
   Integer numerator() const
   {
      return mpq_numref(rep);
   }

   Integer denominator() const
   {
      return mpq_denref(rep);
   }

   /// @group Set the numerator, denominator, or both.
   void set(const Integer& num, const Integer& den)
   {
      mpq_set_num(rep, num.rep);
      mpq_set_den(rep, den.rep);
      mpq_canonicalize(rep);
   }

   void set_numerator(const Integer& num)
   {
      mpq_set_num(rep, num.rep);
      mpq_canonicalize(rep);
   }

   /// @endgroup
   void set_denominator(const Integer& den)
   {
      mpq_set_den(rep, den.rep);
      mpq_canonicalize(rep);
   }

   /// Uses \\mpq_set\\.
   Rational& operator= (const Rational& a)
   {
      mpq_set(rep, a.rep);
      return *this;
   }

   /// @group Assigns an integer number.
   Rational& operator= (const Integer& a)
   {
      mpz_set(mpq_numref(rep), a.rep);
      mpz_set_ui(mpq_denref(rep), 1);
      return *this;
   }

   /// Uses \\mpq_set_si\\.
   Rational& operator= (long i)
   {
      mpq_set_si(rep, i, 1);
      return *this;
   }

   /** Uses \\mpq_set_si\\.
       @endgroup */
   Rational& operator= (int i)
   {
      mpq_set_si(rep, i, 1);
      return *this;
   }

   /// Assigns the next approximation of a floating-point number.
   Rational& operator= (double d);

   /// Swaps the values.
   void swap (Rational& b)
   {
      __STD::swap(*rep, *b.rep);
   }

   /// Uses \\mpz_neg\\.
   Rational& negate()
   {
      mpz_neg(mpq_numref(rep), mpq_numref(rep));
      return *this;
   }

   /// Uses \\mpq_add\\.
   Rational& operator+= (const Rational& a) 
   {
      mpq_add(rep, rep, a.rep);
      return *this;
   }

   /// uses \\mpq_sub\\.
   Rational& operator-= (const Rational& a) 
   {
      mpq_sub(rep, rep, a.rep);
      return *this;
   }

   /// Uses \\mpq_mul\\.
   Rational& operator*= (const Rational& a) 
   {
      mpq_mul(rep, rep, a.rep);
      return *this;
   }

   /// Uses \\mpq_div\\.
   Rational& operator/=(const Rational& a) 
   {
      mpq_div(rep, rep, a.rep);
      return *this;
   }
   
   /// Multiplies with a power of 2, uses \\mpz_mul_2exp\\.
   Rational& operator<<= (unsigned long a)
   { 
      mpz_mul_2exp(mpq_numref(rep), mpq_numref(rep), a);
      mpq_canonicalize(rep);
      return *this;
   }
   
   /// Divides thru a power of 2, uses \\mpz_mul_2exp\\.
   Rational& operator>>= (unsigned long a)
   {
      mpz_mul_2exp(mpq_denref(rep), mpq_denref(rep), a);
      mpq_canonicalize(rep);
      return *this;
   }

   friend Rational operator+ (const Rational& a, const Rational& b);
   friend Rational operator- (const Rational& a, const Rational& b); 
   friend Rational operator- (const Rational& a); 
   friend Rational operator* (const Rational& a, const Rational& b);
   friend Rational operator/ (const Rational& a, const Rational& b); 
   friend Rational operator<<(const Rational& a, unsigned long b);
   friend Rational operator>>(const Rational& a, unsigned long b);

   /// Uses \\mpq_sgn\\.
   friend int sign(const Rational& a)
   {
      return mpq_sgn(a.rep);
   }

   /// Comparison with 0, uses \\mpq_sgn\\.
   bool operator!() const
   {
      return !mpq_sgn(rep);
   }

   /// Comparison, uses \\mpq_cmp\\.
   int compare(const Rational& b) const { return mpq_cmp(rep, b.rep); }

   friend Rational abs(const Rational& a);
   friend Rational inv(const Rational& a);

   /// Converts Rational to double.
   operator double() const
   {
      return mpq_get_d(rep);
   }

   friend istream& operator>> (istream &input, Rational& a);
};

/// Uses \\mpq_add\\.
inline Rational operator+ (const Rational& a, const Rational& b) return c;
{
   mpq_add(c.rep, a.rep, b.rep);
}

/// Uses \\mpq_sub\\.
inline Rational operator- (const Rational& a, const Rational &b) return c;
{
   mpq_sub(c.rep, a.rep, b.rep);
}

/// Uses \\mpq_neg\\.
inline Rational operator- (const Rational& a) return b;
{
   mpq_neg(b.rep, a.rep);
}

/// Uses \\mpq_mul\\.
inline Rational operator* (const Rational& a, const Rational& b) return c;
{
   mpq_mul(c.rep, a.rep, b.rep);
}

/// Uses \\mpq_div\\.
inline Rational operator/ (const Rational& a, const Rational& b) return c;
{
   mpq_div(c.rep, a.rep, b.rep); 
}

/// Multiplies with a power of 2, uses \\mpz_mul_2exp\\.
inline Rational operator<< (const Rational& a, unsigned long b) return c(a);
{
   mpz_mul_2exp(mpq_numref(c.rep), mpq_numref(c.rep), b);
   mpq_canonicalize(c.rep);
}

/// Divides thru a power of 2, uses \\mpz_mul_2exp\\.
inline Rational operator>> (const Rational& a, unsigned long b) return c(a);
{
   mpz_mul_2exp(mpq_denref(c.rep), mpq_denref(c.rep), b);
   mpq_canonicalize(c.rep);
}

/** @group Relational operators, using \\mpq_cmp\\.
    @belongs_to Rational */
inline bool operator== (const Rational& a, const Rational& b)
{
   return a.compare(b)==0;
}

/** @endgroup
    @belongs_to Rational */
inline bool operator< (const Rational& a, const Rational& b)
{
   return a.compare(b)<0; 
}

/// Absolute value, uses \\mpz_abs\\.
inline Rational abs(const Rational& a) return b(a);
{
   mpz_abs(mpq_numref(b.rep), mpq_numref(b.rep));
}

/// Inversion, uses \\mpq_inv\\.
inline Rational inv(const Rational& a) return b;
{
   mpq_inv(b.rep, a.rep);
}

/** Input from a stream.
    Numerator and denominator are expected to be separated by `/'.
    Omitted denominator is assumed being equal to 1.
*/
istream& operator>> (istream& is, Rational& a);

/** Output in a stream.
    Numerator and denominator are separated by `/'.
    Denominators of integral numbers are suppressed.
    The numeric base and its display are controlled by the \\ios::basefield\\ and
    \\ios::showbase\\ flags of the stream.
    @belongs_to Rational
*/
ostream& operator<< (ostream &os, const Rational& a);

inline Integer& Integer::operator= (const Rational& r)
{
   mpz_tdiv_q(rep, mpq_numref(r.rep), mpq_denref(r.rep));
   return *this;
}

inline Integer::Integer(const Rational& r)
{
   mpz_init(rep);
   *this=r;
}

__STL_BEGIN_NAMESPACE
inline void swap(Rational& r1, Rational& r2) { r1.swap(r2); }
__STL_END_NAMESPACE

#endif // _POLYMAKE_GMP_RATIONAL_H

// Local Variables:
// mode:C++
// End:
