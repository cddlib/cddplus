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

#ifndef _POLYMAKE_GMP_INTEGER_H
#define _POLYMAKE_GMP_INTEGER_H "$Project: polymake $$Id: Integer.h,v 1.32 2003/01/06 17:07:06 gawrilow Exp $"

#if defined(__GNUC__) && __GNUC_MINOR__<3 && !defined(__APPLE__)
#pragma interface
#endif

#include <iostream>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <gmp_init.h>

class Integer; class Rational;

namespace std {
   Integer abs(const Integer&);
   Rational abs(const Rational&);
}

/** Exception type
    A constructor of @@link Integer @@end or @@link Rational @@end from \\const char*\\ throws an exception
    of this type in case of a syntax error.
*/
class gmp_error : public std::domain_error {
public:
   gmp_error(const std::string& what_arg) : std::domain_error(what_arg) { }
};

class temp_Integer : public MP_INT {
protected:
   /// never instantiate this class: it is a pure masquerade
   temp_Integer();
   ~temp_Integer();
};

/** Arbitrary precision integer number.
    It is a wrapper around GMP (GNU Multiple Precision Library) type \\mpz_t\\.
    Developed and tested with GMP versions 3.1.1 and 4.0
    See the GMP Home Pages at `http://www.swox.com/gmp/'.
    @index main
*/
class Integer {
private:
   /// GMP's representation
   mpz_t rep;

public:
   /// Initialize to 0.
   Integer() { mpz_init(rep); }

   Integer(const Integer& i) {
      mpz_init_set(rep, i.rep);
   }

   Integer(long i) {
      mpz_init_set_si(rep, i);
   }

   Integer(int i) {
      mpz_init_set_si(rep, i);
   }

   Integer(double d) {
      mpz_init_set_d(rep, d);
   }

   /// Recognizes automatically number base 10, 8, or 16.
   explicit Integer(const char* s) {
      mpz_init(rep);
      try {
	 set(s);
      }
      catch (const gmp_error&) {
	 mpz_clear(rep);
	 throw;
      }
   }

   explicit Integer(mpz_srcptr src) {
      mpz_init_set(rep, src);
   }

   /// take over the GMP structure without copying
   explicit Integer(temp_Integer& tmp) {
      rep[0]=tmp;
   }

   /// Performs division with rounding via truncation.
   inline Integer(const Rational& r);

   /// Performs division with rounding via truncation.
   inline Integer& operator= (const Rational& r);

   template <class Arg>
   Integer(void (*f)(mpz_ptr,Arg), Arg a) {
      mpz_init(rep);
      f(rep,a);
   }

   template <class Arg1, class Arg2>
   Integer(void (*f)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b) {
      mpz_init(rep);
      f(rep,a,b);
   }

   template <class Arg1, class Arg2>
   Integer(Arg2 (*f)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b) {
      mpz_init(rep);
      f(rep,a,b);
   }

   static
   mpz_srcptr _tmp_negate(mpz_ptr temp, mpz_srcptr src) {
      *temp=*src; mpz_neg(temp,temp); return temp;
   }

   ~Integer() {
      mpz_clear(rep);
   }

   Integer& operator= (const Integer& b) {
      mpz_set(rep, b.rep);
      return *this;
   }

   Integer& operator= (long b) {
      mpz_set_si(rep, b);
      return *this;
   }

   Integer& operator= (int b) {
      return operator=(long(b));
   }

   Integer& operator= (double d) {
      mpz_set_d(rep, d);
      return *this;
   }

   /// Recognizes automatically number base 10, 8, or 16.
   Integer& set(const char *s) throw(gmp_error) {
      if (mpz_set_str(rep, s, 0) < 0) throw gmp_error("Integer: syntax error in string");
      return *this;
   }

   /// for the seldom case of unwrapped GMP objects coexisting with us
   Integer& set(mpz_srcptr src) {
      mpz_set(rep, src);
      return *this;
   }

   operator double() const {
      return mpz_get_d(rep);
   }

   operator long() const throw(gmp_error) {
      if (!mpz_fits_slong_p(rep)) throw gmp_error("Integer: value too big");
      return mpz_get_si(rep);
   }

   operator int() const throw(gmp_error) {
      if (!mpz_fits_sint_p(rep)) throw gmp_error("Integer: value too big");
      return mpz_get_si(rep);
   }

   /// Converts integer to string.
   std::string to_string(int base=10) const;

   /// Swaps the values.
   void swap(Integer& b) {
      mpz_swap(rep, b.rep);
   }

   /** Accelerated combination of copy constructor and destructor.
       Aimed to be used in container classes only! */
   friend void relocate(Integer* from, Integer* to) {
      to->rep[0] = from->rep[0];
   }

   /// Increment by 1.
   Integer& operator++() {
      mpz_add_ui(rep, rep, 1);
      return *this;
   }

   /// Decrement by 1.
   Integer& operator--() {
      mpz_sub_ui(rep, rep, 1);
      return *this;
   }

   /// In-place negation.
   Integer& negate() {
      mpz_neg(rep, rep);
      return *this;
   }

   Integer& operator+= (const Integer& b) {
      mpz_add(rep, rep, b.rep);
      return *this;
   }

   Integer& operator+= (long b) {
      if (b>=0) mpz_add_ui(rep, rep, b);
      else mpz_sub_ui(rep, rep, -b);
      return *this;
   }

   Integer& operator-= (const Integer& b) {
      mpz_sub(rep, rep, b.rep);
      return *this;
   }

   Integer& operator-= (long b) {
      if (b>=0) mpz_sub_ui(rep, rep, b);
      else mpz_add_ui(rep, rep, -b);
      return *this;
   }

   Integer& operator*= (const Integer& b) {
      mpz_mul(rep, rep, b.rep);
      return *this;
   }

   Integer& operator*= (long b) {
      mpz_mul_si(rep, rep, b);
      return *this;
   }

   /// @group Division with rounding via truncation.
   Integer& operator/= (const Integer& b) {
      mpz_tdiv_q(rep, rep, b.rep);
      return *this;
   }

   Integer& operator/= (long b) {
      if (b>=0) {	// the case b==0 is handled within GMP and converted to ZERO DIV exception
	 mpz_tdiv_q_ui(rep, rep, b);
      } else {
	 mpz_tdiv_q_ui(rep, rep, -b);
	 negate();
      }
      return *this;
   }

   Integer& operator%= (const Integer& b) {
      mpz_tdiv_r(rep, rep, b.rep);
      return *this;
   }

   /// @endgroup
   Integer& operator%= (long b) {
      mpz_tdiv_r_ui(rep, rep, b>=0 ? b : -b);
      return *this;
   }

   /// Multiplies with 2<sup>k</sup>.
   Integer& operator<<= (unsigned long k) {
      mpz_mul_2exp(rep, rep, k);
      return *this;
   }

   /// Divides thru 2<sup>k</sup>, rounds via truncation.
   Integer& operator>>= (unsigned long k) {
      mpz_tdiv_q_2exp(rep, rep, k);
      return *this;
   }

   friend Integer operator+ (const Integer& a, const Integer& b) {
      return Integer(mpz_add, a.rep, b.rep);
   }

   friend Integer operator+ (const Integer& a, long b) {
      return Integer(b>=0 ? mpz_add_ui : mpz_sub_ui, a.rep, (unsigned long)(b>=0 ? b : -b));
   }

   friend Integer operator- (const Integer& a, const Integer& b) {
      return Integer(mpz_sub, a.rep, b.rep);
   }

   friend Integer operator- (const Integer& a, long b) {
      return Integer(b>=0 ? mpz_sub_ui : mpz_add_ui, a.rep, (unsigned long)(b>=0 ? b : -b));
   }

   friend Integer operator- (long a, const Integer& b) {
      mpz_t minus_b;
      return Integer(a>=0 ? mpz_add_ui : mpz_sub_ui,
		     _tmp_negate(minus_b,b.rep),  (unsigned long)(a>=0 ? a : -a));
   }

   friend Integer operator- (const Integer& a) {
      return Integer(mpz_neg, a.rep);
   }

   friend Integer operator* (const Integer& a, const Integer& b) {
      return Integer(mpz_mul, a.rep, b.rep);
   }

   friend Integer operator* (const Integer& a, long b) {
      return Integer(mpz_mul_si, a.rep, b);
   }

   /// @group Division with rounding via truncation.
   friend Integer operator/ (const Integer& a, const Integer& b) {
      return Integer(mpz_tdiv_q, a.rep, b.rep);
   }

   friend Integer operator/ (const Integer& a, long b) {
      if (b>=0)
	 return Integer(mpz_tdiv_q_ui, a.rep, (unsigned long)b);
      mpz_t minus_a;
      return Integer(mpz_tdiv_q_ui, _tmp_negate(minus_a,a.rep), (unsigned long)(-b));
   }

   friend int operator/ (int a, const Integer& b) {
      return mpz_fits_sint_p(b.rep) ? a/mpz_get_si(b.rep) : 0;
   }

   friend long operator/ (long a, const Integer& b) {
      return mpz_fits_slong_p(b.rep) ? a/mpz_get_si(b.rep) : 0;
   }

   friend Integer operator% (const Integer& a, const Integer& b) {
      return Integer(mpz_tdiv_r, a.rep, b.rep);
   }

   friend long operator% (const Integer& a, long b) {
      long r=mpz_tdiv_ui(a.rep, b);
      return mpz_sgn(a.rep)>=0 ? r : -r;
   }

   friend Integer operator% (int a, const Integer& b) {
      return mpz_fits_sint_p(b.rep) ? Integer(a%mpz_get_si(b.rep)) : b;
   }

   friend Integer operator% (long a, const Integer& b) {
      return mpz_fits_slong_p(b.rep) ? Integer(a%mpz_get_si(b.rep)) : b;
   }

   /// Multiplies with 2<sup>k</sup>.
   friend Integer operator<< (const Integer& a, unsigned long k) {
      return Integer(mpz_mul_2exp, a.rep, k);
   }

   /// Divides through 2<sup>k</sup>, truncates to 0.
   friend Integer operator>> (const Integer& a, unsigned long k) {
      return Integer(mpz_tdiv_q_2exp, a.rep, k);
   }

   /// Compares with 0.
   bool operator!() const {
      return !mpz_sgn(rep);
   }

   /// Compares with 0.
   operator bool() const {
      return mpz_sgn(rep);
   }

   /// Comparison. The magnitude of the return value is arbitrary, only its sign is relevant.
   int compare(const Integer& b) const {
      return mpz_cmp(rep, b.rep);
   }

   friend bool operator== (const Integer& a, long b) {
      return mpz_fits_slong_p(a.rep) && mpz_get_si(a.rep)==b;
   }

   friend bool operator< (const Integer& a, long b) {
      return mpz_fits_slong_p(a.rep) ? mpz_get_si(a.rep)<b : mpz_sgn(a.rep)<0;
   }

   friend bool operator> (const Integer& a, long b) {
      return mpz_fits_slong_p(a.rep) ? mpz_get_si(a.rep)>b : mpz_sgn(a.rep)>0;
   }

   friend bool abs_equal(const Integer& a, const Integer& b) {
      return !mpz_cmpabs(a.rep, b.rep);
   }

   friend bool abs_equal(const Integer& a, long b) {
      return !mpz_cmpabs_ui(a.rep, std::abs(b));
   }

   /// Factorial function.
   friend Integer fac(unsigned long k) {
      return Integer(mpz_fac_ui, k);
   }

   friend Integer pow(const Integer& a, unsigned long k) {
      return Integer(mpz_pow_ui, a.rep, k);
   }

   /// Integral part of square root.
   friend Integer sqrt(const Integer& a) {
      return Integer(mpz_sqrt, a.rep);
   }

   /// Absolute value.
   friend Integer std::abs(const Integer& a);

   /// Greatest common divisor.
   friend Integer gcd(const Integer& a, const Integer& b) {
      return Integer(mpz_gcd, a.rep, b.rep);
   }

   friend Integer gcd(const Integer& a, long b) {
      return Integer(mpz_gcd_ui, a.rep, (unsigned long)b);
   }

   friend Integer gcd(long a, const Integer& b) {
      return Integer(mpz_gcd_ui, b.rep, (unsigned long)a);
   }

   /// Least common multiple.
   friend Integer lcm(const Integer& a, const Integer& b) {
      return Integer(mpz_lcm, a.rep, b.rep);
   }
#if __GNU_MP_VERSION >=4 
   friend Integer lcm(const Integer& a, long b) {
      return Integer(mpz_lcm_ui, a.rep, (unsigned long)b);
   }

   friend Integer lcm(long a, const Integer& b) {
      return Integer(mpz_lcm_ui, b.rep, (unsigned long)a);
   }
#endif
   /// Extended gcd algorithm: $g=a*p+b*q$.
   friend void gcd_ext(const Integer& a, const Integer& b, Integer& g, Integer& p, Integer& q) {
      mpz_gcdext(g.rep, p.rep, q.rep, a.rep, b.rep);
   }

   /// Divide a/b with presumption that b is a multiple of a.
   friend Integer div_exact(const Integer& a, const Integer& b) {
      return Integer(mpz_divexact, a.rep, b.rep);
   }

   /// @group binomial coefficient \\n\\ choose \\k\\.
   friend Integer binom(const Integer& n, unsigned long k) {
      return Integer(mpz_bin_ui, n.rep, k);
   }

   friend Integer binom(unsigned long n, unsigned long k) {
      return Integer(mpz_bin_uiui, n, k);
   }

protected:
   /** Reads an \\Integer\\ value from an input stream.
       @param is input stream
       @param dst storage for the value
       @param allow_sign whether leading whitespaces and sign are expected
   */
   static void read(std::istream& is, mpz_ptr dst, bool allow_sign=true);

   /** Calculates the size of the buffer needed to store an ASCII representation of an \\Integer\\.
       @param os an output stream whose formatting flags must be taken into account.
   */
   static size_t strsize(const std::ostream& os, mpz_srcptr src);

   /** Produces an ASCII representation of an Integer.
       @param buf buffer of size not less than the return value of strsize().
   */
   static void putstr(const std::ostream& os, char* buf, mpz_srcptr src);
public:
   friend std::istream& operator>> (std::istream& is, Integer& a) {
      read(is,a.rep);
      return is;
   }

   friend std::ostream& operator<< (std::ostream& os, const Integer& a);

   friend std::istream& operator>> (std::istream&, Rational&);
   friend std::ostream& operator<< (std::ostream&, const Rational&);

   friend class Rational;

   mpz_srcptr get_rep() const { return rep; }

   struct div_t;
};

/// Analogous to ::div_t
struct Integer::div_t {
   Integer quot, rem;

   div_t(const Integer& n, const Integer& d) {
      mpz_tdiv_qr(quot.rep, rem.rep, n.rep, d.rep);
   }
};

inline Integer operator+ (const Integer& a) {
   return a;
}

inline const Integer operator++ (Integer& a, int) {
   Integer copy(a); ++a; return copy;
}

inline Integer operator-- (Integer& a, int) {
   Integer copy(a); --a; return copy;
}

inline Integer operator+ (long a, const Integer& b) { return b+a; }
inline Integer operator+ (const Integer& a, int b) { return a+long(b); }
inline Integer operator+ (int a, const Integer& b) { return b+long(a); }

inline Integer operator- (const Integer& a, int b) { return a-long(b); }
inline Integer operator- (int a, const Integer& b) { return long(a)-b; }

inline Integer operator* (long a, const Integer& b) { return b*a; }
inline Integer operator* (const Integer& a, int b) { return a*long(b); }
inline Integer operator* (int a, const Integer& b) { return b*long(a); }

inline Integer operator/ (const Integer& a, int b) { return a/long(b); }
inline Integer operator% (const Integer& a, int b) { return a%long(b); }

inline Integer operator<< (const Integer& a, unsigned int k) {
   return a << static_cast<unsigned long>(k);
}
inline Integer operator>> (const Integer& a, unsigned int k) {
   return a >> static_cast<unsigned long>(k);
}
inline Integer operator<< (const Integer& a, long k) {
   if (k<0) return a >> static_cast<unsigned long>(-k);
   return a << static_cast<unsigned long>(k);
}
inline Integer operator>> (const Integer& a, long k) {
   if (k<0) return a << static_cast<unsigned long>(-k);
   return a >> static_cast<unsigned long>(k);
}
inline Integer operator<< (const Integer& a, int k) { return a << long(k); }
inline Integer operator>> (const Integer& a, int k) { return a >> long(k); }

inline bool operator== (const Integer& a, const Integer& b) { return a.compare(b)==0; }
inline bool operator!= (const Integer& a, const Integer& b) { return a.compare(b)!=0; }
inline bool operator< (const Integer& a, const Integer& b) { return a.compare(b)<0; }
inline bool operator> (const Integer& a, const Integer& b) { return a.compare(b)>0; }
inline bool operator<= (const Integer& a, const Integer& b) { return a.compare(b)<=0; }
inline bool operator>= (const Integer& a, const Integer& b) { return a.compare(b)>=0; }

inline bool operator== (const temp_Integer& a, const temp_Integer& b) { return mpz_cmp(&a,&b)==0; }
inline bool operator!= (const temp_Integer& a, const temp_Integer& b) { return mpz_cmp(&a,&b)!=0; }
inline bool operator< (const temp_Integer& a, const temp_Integer& b) { return mpz_cmp(&a,&b)<0; }
inline bool operator> (const temp_Integer& a, const temp_Integer& b) { return mpz_cmp(&a,&b)>0; }
inline bool operator<= (const temp_Integer& a, const temp_Integer& b) { return mpz_cmp(&a,&b)<=0; }
inline bool operator>= (const temp_Integer& a, const temp_Integer& b) { return mpz_cmp(&a,&b)>=0; }

inline bool operator== (long a, const Integer& b) { return b==a; }
inline bool operator== (const Integer& a, int b) { return a==long(b); }
inline bool operator== (int a, const Integer& b) { return b==long(a); }
inline bool abs_equal  (long a, const Integer& b) { return abs_equal(b,a); }
inline bool abs_equal  (const Integer& a, int b) { return abs_equal(a,long(b)); }
inline bool abs_equal  (int a, const Integer& b) { return abs_equal(b,long(a)); }

inline bool operator!= (const Integer& a, long b) { return !(a==b); }
inline bool operator!= (const Integer& a, int b) { return !(a==long(b)); }
inline bool operator!= (long a, const Integer& b) { return !(b==a); }
inline bool operator!= (int a, const Integer& b) { return !(b==long(a)); }

inline bool operator< (const Integer& a, int b) { return a<long(b); }
inline bool operator< (long a, const Integer& b) { return b>a; }
inline bool operator< (int a, const Integer& b) { return b>long(a); }

inline bool operator> (const Integer& a, int b) { return a>long(b); }
inline bool operator> (long a, const Integer& b) { return b<a; }
inline bool operator> (int a, const Integer& b) { return b<long(a); }

inline bool operator<= (const Integer& a, long b) { return !(a>b); }
inline bool operator<= (const Integer& a, int b) { return !(a>long(b)); }
inline bool operator<= (long a, const Integer& b) { return !(b<a); }
inline bool operator<= (int a, const Integer& b) { return !(b<long(a)); }

inline bool operator>= (const Integer& a, long b) { return !(a>b); }
inline bool operator>= (const Integer& a, int b) { return !(a>long(b)); }
inline bool operator>= (long a, const Integer& b) { return !(b>a); }
inline bool operator>= (int a, const Integer& b) { return !(b>long(a)); }

namespace std {
   inline void swap(Integer& i1, Integer& i2) { i1.swap(i2); }

   inline Integer abs(const Integer& a) {
      return Integer(mpz_abs, a.rep);
   }

   inline Integer::div_t div(const Integer& n, const Integer& d) {
      return Integer::div_t(n,d);
   }
}
namespace std_ext {
   template <> struct hash<Integer> : hash<MP_INT> {
      size_t operator() (const Integer& a) const {
	 return _do(a.get_rep());
      }
   };
}

#endif // _POLYMAKE_GMP_INTEGER_H

// Local Variables:
// mode:C++
// End:
