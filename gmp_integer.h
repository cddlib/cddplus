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

#ifndef _POLYMAKE_GMP_INTEGER_H
#define _POLYMAKE_GMP_INTEGER_H

//! project = "$Project: polymake $"
//!   rcsid = "$Id: gmp_integer.h,v 1.1 1999/02/08 14:15:20 gawrilow Exp $"

#pragma interface

#include <iostream>
#include <string>
#include <algobase.h>
extern "C" {
#include <gmp.h>
}

class Rational;

/** @group Reads an \\Integer\\ value from an input stream.
    @param is input stream
    @param dst storage for the value
    @param delim delimiter character (additionally to whitespaces)
    @param allow_sign whether leading whitespaces and sign are expected
    @belongs_to Integer
*/
void read(istream& is, mpz_t dst, char delim=0, bool allow_sign=true);

/** Multiple Precision Integer class wrapping GMP integer type mpz_t
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
class Integer {
private:
   /// GMP's representation
   mpz_t rep;

   /// Gets a ready \\mpz_t\\ instance.
   Integer(const MP_INT *arg)
   {
      mpz_init_set(rep, arg);
   }
public:
  
   /// Initializes to 0.
   Integer()
   {
      mpz_init(rep);
   }

   /// Copy constructor uses \\mpz_init_set\\.
   Integer(const Integer& i)
   {
      mpz_init_set(rep, i.rep);
   }

   /// Uses \\mpz_init_set_si\\.
   explicit Integer(long i)
   {   
      mpz_init_set_si(rep, i);
   }

   /// Uses \\mpz_init_set_si\\.
   explicit Integer(int i)
   {   
      mpz_init_set_si(rep, i);
   }

   /// Uses \\mpz_init_set_d\\.
   explicit Integer(double d)
   {   
      mpz_init_set_d(rep, d);
   }

   /// Integer division - uses \\mpz_tdiv_q\\.
   inline explicit Integer(const Rational& r);

   /// Integer division - uses \\mpz_tdiv_q\\.
   inline Integer& operator= (const Rational& r);

   ~Integer()
   {
      mpz_clear(rep);
   }

   /// Uses \\mpz_set\\.
   Integer& operator= (const Integer& a)
   {
      mpz_set (rep, a.rep);
      return *this;
   }

   /// Uses \\mpz_set_si\\.
   Integer& operator= (long i)
   {
      mpz_set_si(rep, i);
      return *this;
   }

   /// Uses \\mpz_set_si\\.
   Integer& operator= (int i)
   {
      mpz_set_si(rep, i);
      return *this;
   }

   /// Uses \\mpz_set_d\\.
   Integer& operator= (double d)
   {
      mpz_set_d(rep, d);
      return *this;
   }

   /// Swaps the values.
   void swap(Integer& b)
   {
      __STD::swap(*rep, *b.rep);
   }

   /// Uses \\mpz_add_ui\\.
   Integer& operator++()
   {
      mpz_add_ui(rep, rep, 1);
      return *this;
   }

   friend Integer operator++ (Integer& a, int);

   /// Uses \\mpz_sub_ui\\.
   Integer& operator--()
   {
      mpz_sub_ui(rep, rep, 1);
      return *this;
   }

   friend Integer operator-- (Integer& a, int);

   /// Uses \\mpz_neg\\.
   Integer& negate()
   {
      mpz_neg(rep, rep);
      return *this;
   }

   /// Uses \\mpz_add\\.
   Integer& operator+=(const Integer& a) 
   {
      mpz_add(rep, rep, a.rep);
      return *this;
   }

   /// Uses \\mpz_sub\\.
   Integer& operator-=(const Integer& a) 
   {
      mpz_sub(rep, rep, a.rep);
      return *this;
   }

   /// Uses \\mpz_mul\\.
   Integer& operator*=(const Integer& a) 
   {
      mpz_mul(rep, rep, a.rep);
      return *this;
   }

   /// Uses \\mpz_tdiv_q\\.
   Integer& operator/=(const Integer& a) 
   {
      mpz_tdiv_q(rep, rep, a.rep);
      return *this;
   }

   /// Uses \\mpz_tdiv_r\\.
   Integer& operator%=(const Integer& a) 
   {
      mpz_tdiv_r(rep, rep, a.rep);
      return *this;
   }

   /// Multiplies with a power of 2, uses \\mpz_mul_2exp\\.
   Integer& operator<<= (unsigned long a)
   {
      mpz_mul_2exp(rep, rep, a);
      return *this;
   }

   /// Divides thru a power of 2, uses \\mpz_tdiv_q_2exp\\.
   Integer& operator>>= (unsigned long a)
   {
      mpz_tdiv_q_2exp(rep, rep, a);
      return *this;
   }

   friend Integer operator+ (const Integer& a, const Integer& b);
   friend Integer operator- (const Integer& a, const Integer& b); 
   friend Integer operator- (const Integer& a); 
   friend Integer operator* (const Integer& a, const Integer& b); 
   friend Integer operator/ (const Integer& a, const Integer& b); 
   friend Integer operator% (const Integer& a, const Integer& b); 
   friend Integer operator<< (const Integer& a, unsigned long b);
   friend Integer operator>> (const Integer& a, unsigned long b);

   /// Uses \\mpz_sgn\\.
   friend int sign(const Integer& a)
   {
      return mpz_sgn(a.rep);
   }

   /// Compares with 0, uses \\mpz_sgn\\.
   bool operator!() const
   {
      return !mpz_sgn(rep);
   }

   /// Comparison, uses \\mpz_cmp\\.
   int compare(const Integer& b) const { return mpz_cmp(rep, b.rep); }

   friend Integer fac(unsigned long k);
   friend Integer pow(const Integer& a, unsigned long k); 
   friend Integer sqrt(const Integer& a);
   friend Integer abs(const Integer& a);
   friend Integer gcd(const Integer&a, const Integer& b);
   friend void gcd_ext(const Integer&a, const Integer& b, Integer& c, Integer& p, Integer& q);

   /// Output in a stream, uses \\mpz_get_str\\.
   friend ostream& operator<<(ostream& os, const Integer& a);
   friend istream& operator>>(istream& is, Integer& a);
   friend class Rational;
};

/// Uses \\mpz_add_ui\\.
inline Integer operator++ (Integer& a, int) return temp(a);
{
   mpz_add_ui(a.rep, a.rep, 1);
}

/// Uses \\mpz_sub_ui\\.
inline Integer operator-- (Integer& a, int) return temp(a);
{
   mpz_sub_ui(a.rep, a.rep, 1);
}

/// Uses \\mpz_add\\.
inline Integer operator+ (const Integer& a, const Integer& b) return c;
{
   mpz_add(c.rep, a.rep, b.rep);
}

/// Uses \\mpz_sub\\.
inline Integer operator- (const Integer& a, const Integer &b) return c;
{
   mpz_sub(c.rep, a.rep, b.rep);
}

/// Uses \\mpz_neg\\.
inline Integer operator- (const Integer& a) return b;
{
   mpz_neg(b.rep, a.rep);
}

/// Uses \\mpz_mul\\.
inline Integer operator* (const Integer& a, const Integer& b) return c;
{
   mpz_mul(c.rep, a.rep, b.rep);
}

/// Uses \\mpz_tdiv_q\\. Quotient is rounded towards 0.
inline Integer operator/ (const Integer& a, const Integer& b) return c;
{
   mpz_tdiv_q(c.rep, a.rep, b.rep);
}

/// Uses \\mpz_tdiv_r\\. Quotient is rounded towards 0.
inline Integer operator% (const Integer& a, const Integer& b) return c;
{
   mpz_tdiv_r(c.rep, a.rep, b.rep);
}

/// Uses \\mpz_mul_2exp\\.
inline Integer operator<< (const Integer& a, unsigned long b) return c;
{
   mpz_mul_2exp(c.rep, a.rep, b);
}

/// Uses \\mpz_tdiv_q_2exp\\.
inline Integer operator>> (const Integer& a, unsigned long b) return c;
{
   mpz_tdiv_q_2exp(c.rep, a.rep, b);
}

/** @group Relational operators, using \\mpz_cmp\\.
    @belongs_to Integer */
inline bool operator== (const Integer& a, const Integer& b)
{
   return a.compare(b)==0;
}

/** @endgroup
    @belongs_to Integer */
inline bool operator< (const Integer& a, const Integer& b)
{
   return a.compare(b)<0;
}

/// Factorial function, uses \\mpz_fac_ui\\.
inline Integer fac(unsigned long k) return b;
{
   mpz_fac_ui(b.rep, k);
}

/// Exponentional function, uses \\mpz_pow_ui\\.
inline Integer pow(const Integer& a, unsigned long k) return b;
{
   mpz_pow_ui(b.rep, a.rep, k);
}

/// Integral part of square root, uses \\mpz_sqrt\\.
inline Integer sqrt(const Integer& a) return b;
{
   mpz_sqrt(b.rep, a.rep); 
}

/// Absolute value, uses \\mpz_abs\\.
inline Integer abs(const Integer& a) return b;
{
   mpz_abs(b.rep, a.rep);
}

/// Greatest common divisor, uses \\mpz_gcd\\.
inline Integer gcd(const Integer&a, const Integer& b) return c;
{
   mpz_gcd(c.rep, a.rep, b.rep);
}

/// Extended gcd algorithm: c=a*p+b*q
inline void gcd_ext(const Integer&a, const Integer& b, Integer& c, Integer& p, Integer& q)
{
   mpz_gcdext(c.rep, p.rep, q.rep, a.rep, b.rep);
}

/** Input from a stream, uses \\mpz_set_str\\.
    The numeric base is chosen according to the \\ios::basefield\\ flags of the stream.
    If neither is set, it is detected automatically on the "0" or "0x" prefix.
*/
inline istream& operator>> (istream& is, Integer& a)
{
   read(is,a.rep);
   return is;
}

__STL_BEGIN_NAMESPACE
inline void swap(Integer& i1, Integer& i2) { i1.swap(i2); }
__STL_END_NAMESPACE

struct gmp_alloc_init {
   gmp_alloc_init();
};

//namespace {
//   gmp_alloc_init do_init;
//};

static gmp_alloc_init do_init; //modified by K. Fukuda

#endif // _POLYMAKE_GMP_INTEGER_H

// Local Variables:
// mode:C++
// End:
