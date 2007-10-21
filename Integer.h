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

#ifndef _POLYMAKE_GMP_INTEGER_H
#define _POLYMAKE_GMP_INTEGER_H "$Project: polymake $$Id: Integer.h 7556 2007-01-12 17:36:36Z gawrilow $"

#include <iostream>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <gmp_init.h>
#include <cctype>
#include <limits>

class Integer; class Rational;

/** Exception type
    A constructor of Integer or Rational from const char* throws an exception
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

// forward declarations of some friends
namespace std {
   Integer abs(const Integer&);
   Rational abs(const Rational&);
}

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

#if defined(__GMP_PLUSPLUS__)
   // convert to the GMP's own C++ wrapper
   mpz_class gmp() const
   {
      return mpz_class(rep);
   }
   
   // construct an Integer from the GMP's own C++ wrapper
   Integer(const mpz_class& i)
   {
      mpz_init_set(rep,i.get_mpz_t());
   }
#endif

   /// Initialize to 0.
   Integer() { mpz_init(rep); }

   Integer(const Integer& i)
   {
      mpz_init_set(rep, i.rep);
   }

   Integer(long i)
   {
      mpz_init_set_si(rep, i);
   }

   Integer(int i)
   {
      mpz_init_set_si(rep, i);
   }

   Integer(double d)
   {
      mpz_init_set_d(rep, d);
   }

   /// Recognizes automatically number base 10, 8, or 16.
   explicit Integer(const char* s)
   {
      mpz_init(rep);
      try {
	 set(s);
      }
      catch (const gmp_error&) {
	 mpz_clear(rep);
	 throw;
      }
   }

   explicit Integer(mpz_srcptr src)
   {
      mpz_init_set(rep, src);
   }

   /// take over the GMP structure without copying
   explicit Integer(temp_Integer& tmp)
   {
      rep[0]=tmp;
   }

   /// Performs division with rounding via truncation.
   inline Integer(const Rational& r);

   /// Performs division with rounding via truncation.
   inline Integer& operator= (const Rational& r);

   template <class Arg>
   Integer(void (*f)(mpz_ptr,Arg), Arg a)
   {
      mpz_init(rep);
      f(rep,a);
   }

   template <class Arg1, class Arg2>
   Integer(void (*f)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b)
   {
      mpz_init(rep);
      f(rep,a,b);
   }

   template <class Arg1, class Arg2>
   Integer(Arg2 (*f)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b)
   {
      mpz_init(rep);
      f(rep,a,b);
   }

   static
   mpz_srcptr _tmp_negate(mpz_ptr temp, mpz_srcptr src)
   {
      *temp=*src; mpz_neg(temp,temp); return temp;
   }

   ~Integer() { mpz_clear(rep); }

   Integer& operator= (const Integer& b)
   {
      mpz_set(rep, b.rep);
      return *this;
   }

   Integer& operator= (long b)
   {
      mpz_set_si(rep, b);
      return *this;
   }

   Integer& operator= (int b) { return operator=(long(b)); }

   Integer& operator= (double d)
   {
      mpz_set_d(rep, d);
      return *this;
   }

   /// Recognizes automatically number base 10, 8, or 16.
   Integer& set(const char *s) throw(gmp_error)
   {
      if (mpz_set_str(rep, s, 0) < 0)
	 throw gmp_error("Integer: syntax error in string");
      return *this;
   }

   /// for the seldom case of unwrapped GMP objects coexisting with us
   Integer& set(mpz_srcptr src)
   {
      mpz_set(rep, src);
      return *this;
   }

   operator double() const { return mpz_get_d(rep); }

   operator long() const throw(gmp_error)
   {
      if (!mpz_fits_slong_p(rep))
	 throw gmp_error("Integer: value too big");
      return mpz_get_si(rep);
   }

   operator int() const throw(gmp_error)
   {
      if (!mpz_fits_sint_p(rep))
	 throw gmp_error("Integer: value too big");
      return mpz_get_si(rep);
   }

   /// Converts integer to string.
   std::string to_string(int base=10) const;

   void swap(Integer& b) { mpz_swap(rep, b.rep); }

   /** Accelerated combination of copy constructor and destructor.
       Aimed to be used in container classes only! */
   friend void relocate(Integer* from, Integer* to)
   {
      to->rep[0] = from->rep[0];
   }

   Integer& operator++()
   {
      mpz_add_ui(rep, rep, 1);
      return *this;
   }

   Integer& operator--()
   {
      mpz_sub_ui(rep, rep, 1);
      return *this;
   }

   /// In-place negation.
   Integer& negate()
   {
      mpz_neg(rep, rep);
      return *this;
   }

   Integer& operator+= (const Integer& b)
   {
      mpz_add(rep, rep, b.rep);
      return *this;
   }

   Integer& operator+= (long b)
   {
      if (b>=0) mpz_add_ui(rep, rep, b);
      else mpz_sub_ui(rep, rep, -b);
      return *this;
   }

   Integer& operator-= (const Integer& b)
   {
      mpz_sub(rep, rep, b.rep);
      return *this;
   }

   Integer& operator-= (long b)
   {
      if (b>=0) mpz_sub_ui(rep, rep, b);
      else mpz_add_ui(rep, rep, -b);
      return *this;
   }

   Integer& operator*= (const Integer& b)
   {
      mpz_mul(rep, rep, b.rep);
      return *this;
   }

   Integer& operator*= (long b)
   {
      mpz_mul_si(rep, rep, b);
      return *this;
   }

   /// @group Division with rounding via truncation.
   Integer& operator/= (const Integer& b)
   {
      mpz_tdiv_q(rep, rep, b.rep);
      return *this;
   }

   Integer& operator/= (long b)
   {
      if (b>=0) {	// the case b==0 is handled within GMP and converted to ZERO DIV exception
	 mpz_tdiv_q_ui(rep, rep, b);
      } else {
	 mpz_tdiv_q_ui(rep, rep, -b);
	 negate();
      }
      return *this;
   }

   Integer& operator%= (const Integer& b)
   {
      mpz_tdiv_r(rep, rep, b.rep);
      return *this;
   }

   /// @endgroup
   Integer& operator%= (long b)
   {
      mpz_tdiv_r_ui(rep, rep, b>=0 ? b : -b);
      return *this;
   }

   /// Multiplies with 2<sup>k</sup>.
   Integer& operator<<= (unsigned long k)
   {
      mpz_mul_2exp(rep, rep, k);
      return *this;
   }

   /// Divides thru 2<sup>k</sup>, rounds via truncation.
   Integer& operator>>= (unsigned long k)
   {
      mpz_tdiv_q_2exp(rep, rep, k);
      return *this;
   }

   friend Integer operator+ (const Integer& a, const Integer& b)
   {
      return Integer(mpz_add, a.rep, b.rep);
   }

   friend Integer operator+ (const Integer& a, long b)
   {
      return Integer(b>=0 ? mpz_add_ui : mpz_sub_ui, a.rep, (unsigned long)(b>=0 ? b : -b));
   }

   friend Integer operator- (const Integer& a, const Integer& b)
   {
      return Integer(mpz_sub, a.rep, b.rep);
   }

   friend Integer operator- (const Integer& a, long b)
   {
      return Integer(b>=0 ? mpz_sub_ui : mpz_add_ui, a.rep, (unsigned long)(b>=0 ? b : -b));
   }

   friend Integer operator- (long a, const Integer& b)
   {
      mpz_t minus_b;
      return Integer(a>=0 ? mpz_add_ui : mpz_sub_ui,
		     _tmp_negate(minus_b,b.rep),  (unsigned long)(a>=0 ? a : -a));
   }

   friend Integer operator- (const Integer& a)
   {
      return Integer(mpz_neg, a.rep);
   }

   friend Integer operator* (const Integer& a, const Integer& b)
   {
      return Integer(mpz_mul, a.rep, b.rep);
   }

   friend Integer operator* (const Integer& a, long b)
   {
      return Integer(mpz_mul_si, a.rep, b);
   }

   /// @group Division with rounding via truncation.
   friend Integer operator/ (const Integer& a, const Integer& b)
   {
      return Integer(mpz_tdiv_q, a.rep, b.rep);
   }

   friend Integer operator/ (const Integer& a, long b)
   {
      if (b>=0)
	 return Integer(mpz_tdiv_q_ui, a.rep, (unsigned long)b);
      mpz_t minus_a;
      return Integer(mpz_tdiv_q_ui, _tmp_negate(minus_a,a.rep), (unsigned long)(-b));
   }

   friend int operator/ (int a, const Integer& b)
   {
      return mpz_fits_sint_p(b.rep) ? a/mpz_get_si(b.rep) : 0;
   }

   friend long operator/ (long a, const Integer& b)
   {
      return mpz_fits_slong_p(b.rep) ? a/mpz_get_si(b.rep) : 0;
   }

   friend Integer operator% (const Integer& a, const Integer& b)
   {
      return Integer(mpz_tdiv_r, a.rep, b.rep);
   }

   friend long operator% (const Integer& a, long b)
   {
      long r=mpz_tdiv_ui(a.rep, b);
      return mpz_sgn(a.rep)>=0 ? r : -r;
   }

   friend Integer operator% (int a, const Integer& b)
   {
      return mpz_fits_sint_p(b.rep) ? Integer(a%mpz_get_si(b.rep)) : b;
   }

   friend Integer operator% (long a, const Integer& b)
   {
      return mpz_fits_slong_p(b.rep) ? Integer(a%mpz_get_si(b.rep)) : b;
   }

   /// Multiplies with 2<sup>k</sup>.
   friend Integer operator<< (const Integer& a, unsigned long k)
   {
      return Integer(mpz_mul_2exp, a.rep, k);
   }

   /// Divides through 2<sup>k</sup>, truncates to 0.
   friend Integer operator>> (const Integer& a, unsigned long k)
   {
      return Integer(mpz_tdiv_q_2exp, a.rep, k);
   }

   /// Compares with 0.
   bool operator!() const { return !mpz_sgn(rep); }

   /// Compares with 0.
   operator bool() const { return mpz_sgn(rep); }

   /// Comparison. The magnitude of the return value is arbitrary, only its sign is relevant.
   int compare(const Integer& b) const
   {
      return mpz_cmp(rep, b.rep);
   }

   int compare(long b) const
   {
      return mpz_cmp_si(rep, b);
   }

   int compare(double b) const
   {
      return mpz_cmp_d(rep, b);
   }

   friend bool operator== (const Integer& a, long b)
   {
      return mpz_fits_slong_p(a.rep) && mpz_get_si(a.rep)==b;
   }

   friend bool operator< (const Integer& a, long b)
   {
      return mpz_fits_slong_p(a.rep) ? mpz_get_si(a.rep)<b : mpz_sgn(a.rep)<0;
   }

   friend bool operator> (const Integer& a, long b)
   {
      return mpz_fits_slong_p(a.rep) ? mpz_get_si(a.rep)>b : mpz_sgn(a.rep)>0;
   }

   friend bool abs_equal(const Integer& a, const Integer& b)
   {
      return !mpz_cmpabs(a.rep, b.rep);
   }

   friend bool abs_equal(const Integer& a, long b)
   {
      return !mpz_cmpabs_ui(a.rep, std::abs(b));
   }

   static Integer fac(unsigned long k)
   {
      return Integer(mpz_fac_ui, k);
   }

   static Integer pow(const Integer& a, unsigned long k)
   {
      return Integer(mpz_pow_ui, a.rep, k);
   }

   static Integer pow(unsigned long a, unsigned long k)
   {
      return Integer(mpz_ui_pow_ui, a, k);
   }

   friend Integer sqrt(const Integer& a)
   {
      return Integer(mpz_sqrt, a.rep);
   }

   friend Integer std::abs(const Integer& a);

   friend Integer gcd(const Integer& a, const Integer& b)
   {
      return Integer(mpz_gcd, a.rep, b.rep);
   }
   friend Integer gcd(const Integer& a, long b)
   {
      return Integer(mpz_gcd_ui, a.rep, (unsigned long)b);
   }
   friend Integer gcd(long a, const Integer& b)
   {
      return Integer(mpz_gcd_ui, b.rep, (unsigned long)a);
   }

   friend Integer lcm(const Integer& a, const Integer& b)
   {
      return Integer(mpz_lcm, a.rep, b.rep);
   }
#if __GNU_MP_VERSION >=4 
   friend Integer lcm(const Integer& a, long b)
   {
      return Integer(mpz_lcm_ui, a.rep, (unsigned long)b);
   }
   friend Integer lcm(long a, const Integer& b)
   {
      return Integer(mpz_lcm_ui, b.rep, (unsigned long)a);
   }
#endif

   /// Extended gcd algorithm: $g=a*p+b*q$.
   friend void gcd_ext(const Integer& a, const Integer& b, Integer& g, Integer& p, Integer& q)
   {
      mpz_gcdext(g.rep, p.rep, q.rep, a.rep, b.rep);
   }
   friend Integer div_exact(const Integer& a, const Integer& b)
   {
      return Integer(mpz_divexact, a.rep, b.rep);
   }

   static Integer binom(const Integer& n, unsigned long k)
   {
      return Integer(mpz_bin_ui, n.rep, k);
   }

   static Integer binom(unsigned long n, unsigned long k)
   {
      return Integer(mpz_bin_uiui, n, k);
   }

   /// @param allow_sign whether leading whitespaces and sign are expected
   template <typename Traits>
   void read(std::basic_istream<char, Traits>& is, bool allow_sign=true);

   /// Calculates the size of the buffer needed to store an ASCII representation of an \\Integer\\.
   size_t strsize(std::ios::fmtflags flags) const;

   /** Produces a printable representation of an Integer.
       @param buf buffer of size not less than the return value of strsize().
   */
   void putstr(std::ios::fmtflags flags, char* buf) const;

   friend class Rational;

   mpz_srcptr get_rep() const { return rep; }

   struct div_t;

   class little_buffer {
   private:
      static char little[256];
      char *buf;
   public:
      explicit little_buffer(int n)
	 : buf(n<256 ? little : new char[n]) { }

      ~little_buffer() { if (buf!=little) delete buf; }

      operator char* () const { return buf; }
   };
};

/// Analogous to ::div_t
struct Integer::div_t {
   Integer quot, rem;

   div_t(const Integer& n, const Integer& d)
   {
      mpz_tdiv_qr(quot.rep, rem.rep, n.rep, d.rep);
   }
};

inline Integer operator+ (const Integer& a) { return a; }

inline const Integer operator++ (Integer& a, int) { Integer copy(a); ++a; return copy; }
inline const Integer operator-- (Integer& a, int) { Integer copy(a); --a; return copy; }

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

inline Integer operator<< (const Integer& a, unsigned int k)
{
   return a << static_cast<unsigned long>(k);
}

inline Integer operator>> (const Integer& a, unsigned int k)
{
   return a >> static_cast<unsigned long>(k);
}

inline Integer operator<< (const Integer& a, long k)
{
   if (k<0) return a >> static_cast<unsigned long>(-k);
   return a << static_cast<unsigned long>(k);
}

inline Integer operator>> (const Integer& a, long k)
{
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
inline bool abs_equal(long a, const Integer& b) { return abs_equal(b,a); }
inline bool abs_equal(const Integer& a, int b) { return abs_equal(a,long(b)); }
inline bool abs_equal(int a, const Integer& b) { return abs_equal(b,long(a)); }

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

template <typename Traits>
void Integer::read(std::basic_istream<char, Traits>& is, bool allow_sign)
{
   std::ios::iostate exc=is.exceptions();
   is.exceptions(std::ios::goodbit);

   std::string text;
   char c=0;

   if (allow_sign) {
      is >> c;		// skip leading whitespaces, consume the sign or the first digit
      switch (c) {
      case '-':
	 text += c;
	 /* NOBREAK */
      case '+':
	 is >> c;	// there may be some spaces after the sign, too
      }
   } else {
      is.get(c);
   }
   if (is.eof()) {
      is.setstate(std::ios::failbit);
   } else {
      bool valid=false;
      int base=0;
      switch (is.flags() & std::ios::basefield) {
      case std::ios::hex:
	 base=16; break;
      case std::ios::oct:
	 base=8;  break;
      case std::ios::dec:
	 base=10; break;
      default:		// the base is to be guessed from the prefix
	 if (c == '0') {
	    is.get(c);
	    if (c == 'x' || c == 'X') {
	       base=16;
	       is.get(c);
	    } else {
	       text += '0';
	       valid=true;
	       base=8;
	    }
	 } else {
	    base=10;
	 }
      }

      // gather all feasible characters
      while (!is.eof()) {
	 if (isdigit(c)
	     ? (base==8 && c > '7')
	     : !isalpha(c) || base != 16 || (isupper(c) ? c > 'F' : c > 'f')) {
	    is.unget();
	    break;
	 }
	 text += c;
	 valid=true;
	 is.get(c);
      }

      if (valid) {
	 mpz_set_str(rep, text.c_str(), base);
	 is.clear(is.rdstate() & std::ios::eofbit);
      } else {
	 is.setstate(std::ios::failbit);
      }
   }
   is.exceptions(exc);
}

template <typename Traits>
std::basic_ostream<char, Traits>&
operator<< (std::basic_ostream<char, Traits>& os, const Integer& a)
{
   int s=a.strsize(os.flags());
   Integer::little_buffer buf(s);
   a.putstr(os.flags(), buf);
   return os << static_cast<char*>(buf);
}

template <typename Traits> inline
std::basic_istream<char, Traits>&
operator>> (std::basic_istream<char, Traits>& is, Integer& a)
{
   a.read(is);
   return is;
}

namespace std {

inline void swap(Integer& i1, Integer& i2) { i1.swap(i2); }

inline Integer abs(const Integer& a)
{
   return Integer(mpz_abs, a.rep);
}

inline Integer::div_t div(const Integer& n, const Integer& d)
{
   return Integer::div_t(n,d);
}

template <>
struct numeric_limits<Integer> : numeric_limits<long> {
   static double min() throw() { return -numeric_limits<double>::infinity(); }
   static double max() throw() { return +numeric_limits<double>::infinity(); }
   static const int digits=INT_MAX;
   static const int digits10=INT_MAX;
   static const bool is_bounded=false;
};

} // end namespace std

namespace std_ext {

template <> struct hash<Integer> : hash<MP_INT> {
   size_t operator() (const Integer& a) const { return _do(a.get_rep()); }
};

} // end namespace std_ext

#endif // _POLYMAKE_GMP_INTEGER_H

// Local Variables:
// mode:C++
// End:
