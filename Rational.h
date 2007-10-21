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

#ifndef _POLYMAKE_GMP_RATIONAL_H
#define _POLYMAKE_GMP_RATIONAL_H "$Project: polymake $$Id: Rational.h 7565 2007-01-16 16:29:43Z gawrilow $"

#include <Integer.h>

#if __GNU_MP_VERSION < 4
#define _tmp_little_Integer(x) \
   mpz_t x; \
   mp_limb_t x##_limb; \
   x[0]._mp_alloc=1; \
   x[0]._mp_d=&x##_limb
#endif

class temp_Rational : public MP_RAT {
protected:
   /// never instantiate this class: it is a pure masquerade
   temp_Rational();
   ~temp_Rational();
};

Rational inv(const Rational& a);

/** Arbitrary precision rational number.
    It is a wrapper around GMP (GNU Multiple Precision Library) type \\mpq_t\\.
    Developed and tested with GMP Version 3.1 and higher.
    See the GMP Home Pages at `http://www.swox.com/gmp/'.
    @index main
*/
class Rational {
private:
   /// GMP's representation
   mpq_t rep;

   static void zero_division() __attribute__((noreturn));

   void canonicalize()
   {
      if (!mpz_sgn(mpq_denref(rep))) zero_division();
      mpq_canonicalize(rep);
   }

public:
   
#if defined(__GMP_PLUSPLUS__)
   //constructs from gmp's mpz_class as numerator and denominator
   Rational(const mpz_class& num, const mpz_class& den)
   {
      mpz_init_set(mpq_numref(rep), num.get_mpz_t());
      mpz_init_set(mpq_denref(rep), den.get_mpz_t());
      canonicalize();
   }
#endif

   /// Initializes to 0.
   Rational() { mpq_init(rep); }

   Rational(const Rational& a)
   {
      mpz_init_set(mpq_numref(rep), mpq_numref(a.rep));
      mpz_init_set(mpq_denref(rep), mpq_denref(a.rep));
   }

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

   Rational(const Integer& num, const Integer& den)
   {
      mpz_init_set(mpq_numref(rep), num.rep);
      mpz_init_set(mpq_denref(rep), den.rep);
      canonicalize();
   }

   Rational(long num, long den)
   {
      mpz_init_set_si(mpq_numref(rep), num);
      mpz_init_set_si(mpq_denref(rep), den);
      canonicalize();
   }

   Rational(const Integer& num, long den)
   {
      mpz_init_set(mpq_numref(rep), num.rep);
      mpz_init_set_si(mpq_denref(rep), den);
      canonicalize();
   }

   Rational(long num, const Integer& den)
   {
      mpz_init_set_si(mpq_numref(rep), num);
      mpz_init_set(mpq_denref(rep), den.rep);
      canonicalize();
   }

   Rational(const double& b)
   {
      mpq_init(rep);
      mpq_set_d(rep,b);
   }

   explicit Rational(const char* s)
   {
      mpq_init(rep);
      try {
	 set(s);
      }
      catch (const gmp_error&) {
	 mpq_clear(rep);
	 throw;
      }
   }

   explicit Rational(mpq_srcptr src)
   {
      mpz_init_set(mpq_numref(rep), mpq_numref(src));
      mpz_init_set(mpq_denref(rep), mpq_denref(src));
      canonicalize();
   }

   explicit Rational(mpz_srcptr num_src)
   {
      mpz_init_set(mpq_numref(rep), num_src);
      mpz_init_set_ui(mpq_denref(rep), 1);
   }

   Rational(mpz_srcptr num_src, mpz_srcptr den_src)
   {
      mpz_init_set(mpq_numref(rep), num_src);
      mpz_init_set(mpq_denref(rep), den_src);
      canonicalize();
   }

   explicit Rational(temp_Rational& tmp)
   {
      rep[0]=tmp;
      canonicalize();
   }

   explicit Rational(temp_Integer& tmp_num)
   {
      *mpq_numref(rep)=tmp_num;
      mpz_init_set_ui(mpq_denref(rep), 1);
   }

   Rational(temp_Integer& tmp_num, temp_Integer& tmp_den)
   {
      *mpq_numref(rep)=tmp_num;
      *mpq_denref(rep)=tmp_den;
      canonicalize();
   }

   Rational(void (*f)(mpq_ptr,mpq_srcptr), mpq_srcptr a)
   {
      mpq_init(rep);
      f(rep,a);
   }

   template <class Arg>
   Rational(void (*f)(mpq_ptr,mpq_srcptr,Arg), mpq_srcptr a, Arg b)
   {
      mpq_init(rep);
      f(rep,a,b);
   }

   template <class Arg>
   Rational(void (*numf)(mpz_ptr,Arg), Arg a, mpz_srcptr den)
   {
      mpz_init(mpq_numref(rep));
      numf(mpq_numref(rep),a);
      mpz_init_set(mpq_denref(rep),den);
   }

   template <class Arg>
   Rational(mpz_srcptr num, void (*denf)(mpz_ptr,Arg), Arg a)
   {
      mpz_init_set(mpq_numref(rep),num);
      mpz_init(mpq_denref(rep));
      denf(mpq_denref(rep),a);
   }

   template <class Arg1, class Arg2>
   Rational(void (*numf)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b, mpz_srcptr den)
   {
      mpz_init(mpq_numref(rep));
      numf(mpq_numref(rep),a,b);
      mpz_init_set(mpq_denref(rep),den);
   }

   template <class Arg1, class Arg2>
   Rational(mpz_srcptr num, void (*denf)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b)
   {
      mpz_init_set(mpq_numref(rep),num);
      mpz_init(mpq_denref(rep));
      denf(mpq_denref(rep),a,b);
   }

   template <class Arg1, class Arg2>
   Rational(void (*numf)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b, unsigned long den)
   {
      mpz_init(mpq_numref(rep));
      numf(mpq_numref(rep),a,b);
      mpz_init_set_ui(mpq_denref(rep),den);
   }

   template <class Arg1, class Arg2>
   Rational(unsigned long num, void (*denf)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b)
   {
      mpz_init_set_ui(mpq_numref(rep),num);
      mpz_init(mpq_denref(rep));
      denf(mpq_denref(rep),a,b);
   }

   template <class Arg1, class Arg2>
   Rational(void (*numf)(mpz_ptr,Arg1,Arg2), mpz_srcptr num, Arg1 a, Arg2 b, mpz_srcptr den)
   {
      mpz_init_set(mpq_numref(rep),num);
      numf(mpq_numref(rep),a,b);
      mpz_init_set(mpq_denref(rep),den);
   }

   template <class Arg1, class Arg2>
   Rational(mpz_srcptr num, void (*denf)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b, mpz_srcptr den)
   {
      mpz_init_set(mpq_numref(rep),num);
      mpz_init_set(mpq_denref(rep),den);
      denf(mpq_denref(rep),a,b);
   }

   template <class Arg1, class Arg2, class Arg3, class Arg4>
   Rational(void (*numf)(mpz_ptr,Arg1,Arg2), Arg1 a, Arg2 b,
	    void (*denf)(mpz_ptr,Arg3,Arg4), Arg3 c, Arg4 d)
   {
      mpq_init(rep);
      numf(mpq_numref(rep),a,b);
      denf(mpq_denref(rep),c,d);
   }

   friend Integer::Integer(const Rational& r);
   friend Integer& Integer::operator= (const Rational& r);

   ~Rational() { mpq_clear(rep); }

protected:
   enum proxy_kind { num, den };

   template <proxy_kind kind, bool _canonicalize=true>
   class proxy : public Integer {
   private:
      void canonicalize() {
	 if (!_canonicalize) return;

	 // The constant 8 is arbitrarily chosen.
	 // Every non-zero double word aligned address would do as well.
	 enum { fict_addr=8 };

	 reinterpret_cast<Rational*>
	    (// address of the num/den field
	      reinterpret_cast<char*>(this)
	      - (// offset from the begin of a mpq_t to appropriate mpz_t
		 reinterpret_cast<char*>(kind==num ? mpq_numref(reinterpret_cast<Rational*>(fict_addr)->rep)
					           : mpq_denref(reinterpret_cast<Rational*>(fict_addr)->rep))
		 -reinterpret_cast<char*>(fict_addr)))
	    ->canonicalize();
      }

      /// both undefined
      proxy(const proxy&);
      proxy();

      friend class Rational;
   public:

      template <class T>
      proxy& operator= (const T& b)
      {
	 Integer::operator=(b);
	 canonicalize();
	 return *this;
      }

      proxy& operator++()
      {
	 Integer::operator++();
	 canonicalize();
	 return *this;
      }

      proxy& operator--()
      {
	 Integer::operator--();
	 canonicalize();
	 return *this;
      }

      template <class T>
      proxy& operator+= (const T& b)
      {
	 Integer::operator+=(b);
	 canonicalize();
	 return *this;
      }

      template <class T>
      proxy& operator-= (const T& b)
      {
	 Integer::operator-=(b);
	 canonicalize();
	 return *this;
      }

      template <class T>
      proxy& operator*= (const T& b)
      {
	 Integer::operator*=(b);
	 canonicalize();
	 return *this;
      }

      template <class T>
      proxy& operator/= (const T& b)
      {
	 Integer::operator/=(b);
	 canonicalize();
	 return *this;
      }

      template <class T>
      proxy& operator%= (const T& b)
      {
	 Integer::operator%=(b);
	 canonicalize();
	 return *this;
      }

      template <class T>
      proxy& operator<<= (const T& b)
      {
	 Integer::operator<<=(b);
	 canonicalize();
	 return *this;
      }

      template <class T>
      proxy& operator>>= (const T& b)
      {
	 Integer::operator>>=(b);
	 canonicalize();
	 return *this;
      }
   };

public:
   friend const proxy<num>& numerator(const Rational& a)
   {
      return *reinterpret_cast<const proxy<num>*>(mpq_numref(a.rep));
   }

   friend proxy<num>& numerator(Rational& a)
   {
      return *reinterpret_cast<proxy<num>*>(mpq_numref(a.rep));
   }

   friend const proxy<den>& denominator(const Rational& a)
   {
      return *reinterpret_cast<const proxy<den>*>(mpq_denref(a.rep));
   }

   friend proxy<den>& denominator(Rational& a)
   {
      return *reinterpret_cast<proxy<den>*>(mpq_denref(a.rep));
   }

   friend proxy<num,false>& numerator_nocanon(Rational& a)
   {
      return *reinterpret_cast<proxy<num,false>*>(mpq_numref(a.rep));
   }

   void set(const Integer& num, const Integer& den)
   {
      set(num.rep, den.rep);
   }

   void set(long num, long den)
   {
      mpz_set_si(mpq_numref(rep),num);
      mpz_set_si(mpq_denref(rep),den);
      canonicalize();
   }

   /** Conversion from a printable representation.
       Numerator and denominator are expected delimited by `/'.
       Omitted denominator assumed equal to 1.
   */
   Rational& set(const char *s) throw(gmp_error);

   Rational& operator= (const Rational& b)
   {
      mpq_set(rep, b.rep);
      return *this;
   }

   /// for the seldom case of unwrapped GMP objects coexisting with us
   Rational& set(mpz_srcptr num_src, mpz_srcptr den_src)
   {
      mpq_set_num(rep, num_src);
      mpq_set_den(rep, den_src);
      canonicalize();
      return *this;
   }

   Rational& set(mpq_srcptr src)
   {
      mpq_set(rep, src);
      return *this;
   }

   Rational& operator= (const Integer& b)
   {
      mpq_set_z(rep, b.rep);
      return *this;
   }

   Rational& operator= (long b)
   {
      mpq_set_si(rep, b, 1);
      return *this;
   }

   Rational& operator= (int b) { return operator=(long(b)); }

   Rational& operator= (const double& b)
   {
      mpq_set_d(rep, b);
      return *this;
   }

   operator double() const
   {
      return mpq_get_d(rep);
   }

   operator long() const
   {
      return static_cast<Integer>(*this);
   }

   operator int() const
   {
      return static_cast<Integer>(*this);
   }

   /// Convert rational to string.
   std::string to_string(int base=10) const;

   /// Swaps the values.
   void swap(Rational& b) { mpq_swap(rep, b.rep); }

   /** Accelerated combination of copy constructor and destructor.
       Aimed to be used in container classes only! */
   friend void relocate(Rational* from, Rational* to)
   {
      to->rep[0] = from->rep[0];
   }

   /// Comparison with 0.
   bool operator!() const { return !mpq_sgn(rep); }

   /// Comparison with 0.
   operator bool() const { return mpq_sgn(rep); }

   Rational& operator++ ()
   {
      mpz_add(mpq_numref(rep), mpq_numref(rep), mpq_denref(rep));
      return *this;
   }

   Rational& operator-- ()
   {
      mpz_sub(mpq_numref(rep), mpq_numref(rep), mpq_denref(rep));
      return *this;
   }

   /// In-place negation.
   Rational& negate()
   {
      mpq_neg(rep, rep);
      return *this;
   }

   Rational& operator+= (const Rational& b)
   {
      mpq_add(rep, rep, b.rep);
      return *this;
   }

   Rational& operator+= (const Integer& b)
   {
#if __GNU_MP_VERSION >= 4
      mpz_addmul(mpq_numref(rep), mpq_denref(rep), b.rep);
#else
      Integer t=denominator(*this)*b;
      mpz_add(mpq_numref(rep), mpq_numref(rep), t.rep);
#endif
      return *this;
   }

   Rational& operator+= (long b)
   {
#if __GNU_MP_VERSION >= 4
      if (b>=0) mpz_addmul_ui(mpq_numref(rep),  mpq_denref(rep), b);
      else      mpz_submul_ui(mpq_numref(rep),  mpq_denref(rep), -b);
#else
      if (b) {
	 mpz_t minus_den;
	 mpz_srcptr den= b>0 ? mpq_denref(rep)
	                     : (b=-b, Integer::_tmp_negate(minus_den,mpq_denref(rep)));
	 mpz_addmul_ui(mpq_numref(rep), den, b);
      }
#endif
      return *this;
   }

   Rational& operator-= (const Rational& b)
   {
      mpq_sub(rep, rep, b.rep);
      return *this;
   }

   Rational& operator-= (const Integer& b)
   {
#if __GNU_MP_VERSION >= 4
      mpz_submul(mpq_numref(rep), mpq_denref(rep), b.rep);
#else
      Integer t=denominator(*this)*b;
      mpz_sub(mpq_numref(rep), mpq_numref(rep), t.rep);
#endif
      return *this;
   }

   Rational& operator-= (long b)
   {
#if __GNU_MP_VERSION >= 4
      if (b>=0) mpz_submul_ui(mpq_numref(rep), mpq_denref(rep), b);
      else      mpz_addmul_ui(mpq_numref(rep), mpq_denref(rep), -b);
#else
      if (b) {
	 mpz_t minus_den;
	 mpz_srcptr den= b>0 ? Integer::_tmp_negate(minus_den,mpq_denref(rep))
	                     : (b=-b, mpq_denref(rep));
	 mpz_addmul_ui(mpq_numref(rep), den, b);
      }
#endif
      return *this;
   }

   Rational& operator*= (const Rational& b)
   {
      mpq_mul(rep, rep, b.rep);
      return *this;
   }

   Rational& operator*= (const Integer& b);

   Rational& operator*= (long b)
   {
      if (!*this) return *this;
      if (b) {
#if __GNU_MP_VERSION >= 4
	 unsigned long g=mpz_gcd_ui(0, mpq_denref(rep), b>=0 ? b : -b);
	 if (g==1) {
	    mpz_mul_si(mpq_numref(rep), mpq_numref(rep), b);
	 } else {
	    mpz_mul_si(mpq_numref(rep), mpq_numref(rep), b/g);
	    mpz_divexact_ui(mpq_denref(rep), mpq_denref(rep), g);
	 }
#else
	 _tmp_little_Integer(g);
	 if (mpz_gcd_ui(g, mpq_denref(rep), b>=0 ? b : -b) == 1) {
	    mpz_mul_si(mpq_numref(rep), mpq_numref(rep), b);
	 } else {
	    mpz_mul_si(mpq_numref(rep), mpq_numref(rep), b/g_limb);
	    mpz_divexact(mpq_denref(rep), mpq_denref(rep), g);
	 }
#endif
      } else {
	 *this=0;
      }
      return *this;
   }

   Rational& operator/= (const Rational& b)
   {
      mpq_div(rep, rep, b.rep);
      return *this;
   }

   Rational& operator/= (const Integer& b);

   Rational& operator/= (long b)
   {
      if (!b) zero_division();
      if (*this) {
	 const long babs=b>0 ? b : -b;
#if __GNU_MP_VERSION >= 4
	 unsigned long g=mpz_gcd_ui(0, mpq_numref(rep), babs);
	 if (g==1) {
	    mpz_mul_ui(mpq_denref(rep), mpq_denref(rep), babs);
	 } else {
	    mpz_mul_ui(mpq_denref(rep), mpq_denref(rep), babs/g);
	    mpz_divexact_ui(mpq_numref(rep), mpq_numref(rep), g);
	 }
#else
	 _tmp_little_Integer(g);
	 if (mpz_gcd_ui(g, mpq_numref(rep), babs) == 1) {
	    mpz_mul_ui(mpq_denref(rep), mpq_denref(rep), babs);
	 } else {
	    mpz_mul_ui(mpq_denref(rep), mpq_denref(rep), babs/g_limb);
	    mpz_divexact(mpq_numref(rep), mpq_numref(rep), g);
	 }
#endif
	 if (b<0) mpz_neg(mpq_numref(rep), mpq_numref(rep));
      }
      return *this;
   }

   /// Multiply with $2^k$.
   Rational& operator<<= (unsigned long k)
   {
#if __GNU_MP_VERSION >= 4
      mpq_mul_2exp(rep, rep, k);
#else
      if (mpq_sgn(rep)) {
	 unsigned long den0s=mpz_scan1(mpq_denref(rep), 0);
	 if (k<=den0s) {
	    mpz_tdiv_q_2exp(mpq_denref(rep), mpq_denref(rep), k);
	 } else {
	    mpz_tdiv_q_2exp(mpq_denref(rep), mpq_denref(rep), den0s);
	    mpz_mul_2exp(mpq_numref(rep), mpq_numref(rep), k-den0s);
	 }
      }
#endif
      return *this;
   }
   
   /// Divide thru $2^k$.
   Rational& operator>>= (unsigned long k)
   {
#if __GNU_MP_VERSION >= 4
      mpq_div_2exp(rep, rep, k);
#else
      if (mpz_sgn(mpq_numref(rep))) {
	 unsigned long num0s=mpz_scan1(mpq_numref(rep), 0);
	 if (k<=num0s) {
	    mpz_tdiv_q_2exp(mpq_numref(rep), mpq_numref(rep), k);
	 } else {
	    mpz_tdiv_q_2exp(mpq_numref(rep), mpq_numref(rep), num0s);
	    mpz_mul_2exp(mpq_denref(rep), mpq_denref(rep), k-num0s);
	 }
      }
#endif
      return *this;
   }

   friend Rational operator+ (const Rational& a, const Rational& b)
   {
      return Rational(mpq_add, a.rep, b.get_rep());
   }

   friend Rational operator+ (const Rational& a, const Integer& b)
   {
#if __GNU_MP_VERSION >= 4
      return Rational(mpz_addmul, mpq_numref(a.rep), mpq_denref(a.rep), b.get_rep(), mpq_denref(a.rep));
#else
      Integer t=denominator(a)*b;
      return Rational(mpz_add, mpq_numref(a.rep), t.get_rep(), mpq_denref(a.rep));
#endif
   }

   friend Rational operator+ (const Rational& a, long b)
   {
#if __GNU_MP_VERSION >= 4
      if (b>=0) return Rational(mpz_addmul_ui, mpq_numref(a.rep), mpq_denref(a.rep), (unsigned long)b,
				mpq_denref(a.rep));
      return Rational(mpz_submul_ui, mpq_numref(a.rep), mpq_denref(a.rep), (unsigned long)(-b),
		      mpq_denref(a.rep));
#else
      mpz_t minus_den;
      mpz_srcptr den= b>=0 ? mpq_denref(a.rep)
	                   : (b=-b, Integer::_tmp_negate(minus_den,mpq_denref(a.rep)));
      return Rational(mpz_addmul_ui, mpq_numref(a.rep), den, (unsigned long)b, mpq_denref(a.rep));
#endif
   }

   friend Rational operator- (const Rational& a, const Rational& b)
   {
      return Rational(mpq_sub, a.rep, b.rep);
   }

   friend Rational operator- (const Rational& a, const Integer& b)
   {
#if __GNU_MP_VERSION >= 4
      return Rational(mpz_submul, mpq_numref(a.rep), mpq_denref(a.rep), b.get_rep(), mpq_denref(a.rep));
#else
      Integer t=denominator(a)*b;
      return Rational(mpz_sub, mpq_numref(a.rep), t.get_rep(), mpq_denref(a.rep));
#endif
   }

   friend Rational operator- (const Integer& a, const Rational& b)
   {
#if __GNU_MP_VERSION >= 4
      mpz_t minus_num;
      return Rational(mpz_addmul, Integer::_tmp_negate(minus_num,mpq_numref(b.rep)), mpq_denref(b.rep), a.get_rep(),
		      mpq_denref(b.rep));
#else
      Integer t=denominator(b)*a;
      return Rational(mpz_sub, t.get_rep(), mpq_numref(b.rep), mpq_denref(b.rep));
#endif
   }

   friend Rational operator- (const Rational& a, long b)
   {
#if __GNU_MP_VERSION >= 4
      if (b>=0) return Rational(mpz_submul_ui, mpq_numref(a.rep), mpq_denref(a.rep), (unsigned long)b,
				mpq_denref(a.rep));
      return Rational(mpz_addmul_ui, mpq_numref(a.rep), mpq_denref(a.rep), (unsigned long)(-b),
		      mpq_denref(a.rep));
#else
      mpz_t minus_den;
      mpz_srcptr den= b>=0 ? Integer::_tmp_negate(minus_den,mpq_denref(a.rep))
	                   : (b=-b, mpq_denref(a.rep));
      return Rational(mpz_addmul_ui, mpq_numref(a.rep), den, (unsigned long)b, mpq_denref(a.rep));
#endif
   }

   friend Rational operator- (long a, const Rational& b)
   {
      mpz_t minus_num;
      Integer::_tmp_negate(minus_num,mpq_numref(b.rep));
#if __GNU_MP_VERSION >= 4
      if (a>=0) return Rational(mpz_addmul_ui, minus_num, mpq_denref(b.rep), (unsigned long)a, mpq_denref(b.rep));
      return Rational(mpz_submul_ui, minus_num, mpq_denref(b.rep), (unsigned long)(-a), mpq_denref(b.rep));
#else
      mpz_t minus_den;
      mpz_srcptr den= a>=0 ? mpq_denref(b.rep)
	                   : (a=-a, Integer::_tmp_negate(minus_den,mpq_denref(b.rep)));
      return Rational(mpz_addmul_ui, minus_num, den, (unsigned long)a, mpq_denref(b.rep));
#endif
   }

   friend Rational operator- (const Rational& a)
   {
      return Rational(mpq_neg,a.rep);
   }

   friend Rational operator* (const Rational& a, const Rational& b)
   {
      return Rational(mpq_mul,a.rep,b.rep);
   }

   friend Rational operator* (const Rational& a, const Integer& b);

   friend Rational operator* (const Rational& a, long b)
   {
      if (!b || !a) return Rational();
#if __GNU_MP_VERSION >= 4
      unsigned long g=mpz_gcd_ui(0, mpq_denref(a.rep), b>=0 ? b : -b);
      if (g==1)
	 return Rational(mpz_mul_si, mpq_numref(a.rep), b, mpq_denref(a.rep));
      return Rational(mpz_mul_si, mpq_numref(a.rep), b/(long)g,
		      mpz_divexact_ui, mpq_denref(a.rep), g);
#else
      _tmp_little_Integer(g);
      if (mpz_gcd_ui(g, mpq_denref(a.rep), b>=0 ? b : -b) == 1)
	 return Rational(mpz_mul_si, mpq_numref(a.rep), b, mpq_denref(a.rep));
      return Rational(mpz_mul_si, mpq_numref(a.rep), b/(long)g_limb,
		      mpz_divexact, mpq_denref(a.rep), const_cast<mpz_srcptr>(g));
#endif
   }

   friend Rational operator/ (const Rational& a, const Rational& b)
   {
      return Rational(mpq_div,a.rep,b.rep);
   }

   friend Rational operator/ (const Rational& a, long b)
   {
      if (!b) zero_division();
      if (!a) return Rational();
      mpz_t minus_num;
      mpz_srcptr num= b>0 ? mpq_numref(a.rep)
	                  : (b=-b, Integer::_tmp_negate(minus_num,mpq_numref(a.rep)));
#if __GNU_MP_VERSION >= 4
      unsigned long g=mpz_gcd_ui(0, mpq_numref(a.rep), b);
      if (g==1)
	 return Rational(num, mpz_mul_ui, mpq_denref(a.rep), (unsigned long)b);
      return Rational(mpz_divexact_ui, num, g,
		      mpz_mul_ui, mpq_denref(a.rep), b/g);
#else
      _tmp_little_Integer(g);
      if (mpz_gcd_ui(g, mpq_numref(a.rep), b) == 1)
	 return Rational(num, mpz_mul_ui, mpq_denref(a.rep), (unsigned long)b);
      return Rational(mpz_divexact, num, const_cast<mpz_srcptr>(g),
		      mpz_mul_ui, mpq_denref(a.rep), b/g_limb);
#endif
   }

   friend Rational operator/ (long a, const Rational& b)
   {
      if (!b) zero_division();
      if (!a) return Rational();
      mpz_t num_abs;
      mpz_srcptr num= mpz_sgn(mpq_numref(b.rep))>0 ? mpq_numref(b.rep)
	                                           : (a=-a, Integer::_tmp_negate(num_abs,mpq_numref(b.rep)));
#if __GNU_MP_VERSION >= 4
      unsigned long g=mpz_gcd_ui(0, mpq_numref(b.rep), a>=0 ? a : -a);
      if (g==1)
	 return Rational(mpz_mul_si, mpq_denref(b.rep), a, num);
      return Rational(mpz_mul_si, mpq_denref(b.rep), a/(long)g,
		      mpz_divexact_ui, num, g);
#else
      _tmp_little_Integer(g);
      if (mpz_gcd_ui(g, mpq_numref(b.rep), a>=0 ? a : -a) == 1)
	 return Rational(mpz_mul_si, mpq_denref(b.rep), a, num);
      return Rational(mpz_mul_si, mpq_denref(b.rep), a/(long)g_limb,
		      mpz_divexact, num, const_cast<mpz_srcptr>(g));
#endif
   }

   friend Rational operator/ (const Rational& a, const Integer& b)
   {
      if (!b) zero_division();
      if (!a) return Rational();
      const Integer g=gcd(numerator(a),b);
      if (g==1)
	 return Rational(mpq_numref(a.rep), mpz_mul, mpq_denref(a.rep), b.get_rep());
      const Integer t=div_exact(b,g);
      return Rational(mpz_divexact, mpq_numref(a.rep), g.get_rep(),
		      mpz_mul, mpq_denref(a.rep), t.get_rep());
   }

   friend Rational operator/ (const Integer& a, const Rational& b)
   {
      if (!b) zero_division();
      if (!a) return Rational();
      const Integer g=gcd(a,numerator(b));
      if (g==1)
	 return Rational(mpz_mul, mpq_denref(b.rep), a.get_rep(), mpq_numref(b.rep));
      const Integer t=div_exact(a,g);
      return Rational(mpz_mul, mpq_denref(b.rep), t.get_rep(),
		      mpz_divexact, mpq_numref(b.rep), g.get_rep());
   }

   friend Rational floor(const Rational& a)
   {
      return Rational(mpz_fdiv_q, mpq_numref(a.rep), mpq_denref(a.rep), 1);
   }

   friend Rational ceil(const Rational& a)
   {
      return Rational(mpz_cdiv_q, mpq_numref(a.rep), mpq_denref(a.rep), 1);
   }

   /// Multiply with $2^k$
   friend Rational operator<< (const Rational& a, unsigned long k)
   {
#if __GNU_MP_VERSION >= 4
      return Rational(mpq_mul_2exp, a.rep, k);
#else
      if (!a) return Rational(0);
      unsigned long den0s=mpz_scan1(mpq_denref(a.rep), 0);
      if (k<=den0s)
	 return Rational(mpq_numref(a.rep), mpz_tdiv_q_2exp, mpq_denref(a.rep), k);
      return Rational(mpz_mul_2exp, mpq_numref(a.rep), k-den0s,
		      mpz_tdiv_q_2exp, mpq_denref(a.rep), den0s);
#endif
   }

   /// Divide thru $2^k$
   friend Rational operator>> (const Rational& a, unsigned long k)
   {
#if __GNU_MP_VERSION >= 4
      return Rational(mpq_div_2exp, a.rep, k);
#else
      if (!a) return Rational(0);
      unsigned long num0s=mpz_scan1(mpq_numref(a.rep), 0);
      if (k<=num0s)
	 return Rational(mpz_tdiv_q_2exp, mpq_numref(a.rep), k, mpq_denref(a.rep));
      return Rational(mpz_tdiv_q_2exp, mpq_numref(a.rep), num0s,
		      mpz_mul_2exp, mpq_denref(a.rep), k-num0s);
#endif
   }

   /// Comparison. The magnitude of the return value is arbitrary, only its sign is relevant.
   int compare(const Rational& b) const
   {
      return mpq_cmp(rep, b.rep);
   }

   int compare(long b) const
   {
      return mpz_cmp_ui(mpq_denref(rep),1) ? numerator(*this).compare(b*denominator(*this))
	                                   : numerator(*this).compare(b);
   }

   int compare(const Integer& b) const
   {
      return mpz_cmp_ui(mpq_denref(rep),1) ? numerator(*this).compare(b*denominator(*this))
	                                   : numerator(*this).compare(b);
   }

   friend bool operator== (const Rational& a, const Rational& b)
   {
#if __GNU_MP_VERSION >= 4
      return mpq_equal(a.rep, b.rep);
#else
      return !mpz_cmp(mpq_denref(a.rep), mpq_denref(b.rep)) &&
	     !mpz_cmp(mpq_numref(a.rep), mpq_numref(b.rep));
#endif
   }

   friend bool operator== (const Rational& a, const Integer& b)
   {
      return !mpz_cmp_ui(mpq_denref(a.rep),1) &&
	     !mpz_cmp(mpq_numref(a.rep), b.get_rep());
   }

   friend bool operator== (const Rational& a, long b)
   {
      return !mpz_cmp_ui(mpq_denref(a.rep),1) &&
	     mpz_fits_slong_p(mpq_numref(a.rep)) && mpz_get_si(mpq_numref(a.rep))==b;
   }

   friend bool abs_equal(const Rational& a, const Rational& b)
   {
      return !mpz_cmp(mpq_denref(a.rep), mpq_denref(b.rep)) &&
	     !mpz_cmpabs(mpq_numref(a.rep), mpq_numref(b.rep));
   }
   friend bool abs_equal(const Rational& a, const Integer& b)
   {
      return !mpz_cmp_ui(mpq_denref(a.rep),1) &&
             !mpz_cmpabs(mpq_numref(a.rep), b.get_rep());
   }
   friend bool abs_equal(const Rational& a, long b)
   {
      return !mpz_cmp_ui(mpq_denref(a.rep),1) &&
             mpz_fits_slong_p(mpq_numref(a.rep)) && mpz_get_si(mpq_numref(a.rep))==std::abs(b);
   }

   friend Rational std::abs(const Rational& a);

   friend Rational inv(const Rational& a)
   {
      return Rational(mpq_inv, a.rep);
   }

   mpq_srcptr get_rep() const { return rep; }

   template <typename Traits>
   void read (std::basic_istream<char, Traits>& is)
   {
      numerator_nocanon(*this).read(is);
      if (!is.eof() && is.peek() == '/') {
	 is.ignore();
	 denominator(*this).read(is,false);
	 canonicalize();
      } else {
	 mpz_set_ui(mpq_denref(rep), 1);
      }
   }
};

inline Rational operator+ (const Rational& a) { return a; }

inline bool operator!= (const Rational& a, const Rational& b) { return !(a==b); }
inline bool operator< (const Rational& a, const Rational& b) { return a.compare(b)<0; }
inline bool operator> (const Rational& a, const Rational& b) { return a.compare(b)>0; }
inline bool operator<= (const Rational& a, const Rational& b) { return a.compare(b)<=0; }
inline bool operator>= (const Rational& a, const Rational& b) { return a.compare(b)>=0; }

inline bool operator== (const temp_Rational& a, const temp_Rational& b) { return mpq_equal(&a,&b); }
inline bool operator!= (const temp_Rational& a, const temp_Rational& b) { return !(a==b); }
inline bool operator< (const temp_Rational& a, const temp_Rational& b) { return mpq_cmp(&a,&b)<0; }
inline bool operator> (const temp_Rational& a, const temp_Rational& b) { return mpq_cmp(&a,&b)>0; }
inline bool operator<= (const temp_Rational& a, const temp_Rational& b) { return mpq_cmp(&a,&b)<=0; }
inline bool operator>= (const temp_Rational& a, const temp_Rational& b) { return mpq_cmp(&a,&b)>=0; }

inline bool operator== (const Integer& a, const Rational& b) { return b==a; }
inline bool operator== (long a, const Rational& b) { return b==a; }
inline bool operator== (const Rational& a, int b) { return a==long(b); }
inline bool operator== (int a, const Rational& b) { return b==long(a); }

inline bool abs_equal(const Integer& a, const Rational& b) { return abs_equal(b,a); }
inline bool abs_equal(long a, const Rational& b) { return abs_equal(b,a); }
inline bool abs_equal(const Rational& a, int b) { return abs_equal(a,long(b)); }
inline bool abs_equal(int a, const Rational& b) { return abs_equal(b,long(a)); }

inline bool operator!= (const Rational& a, const Integer& b) { return !(a==b); }
inline bool operator!= (const Integer& a, const Rational& b) { return !(b==a); }
inline bool operator!= (const Rational& a, long b) { return !(a==b); }
inline bool operator!= (long a, const Rational& b) { return !(b==a); }
inline bool operator!= (const Rational& a, int b) { return !(a==b); }
inline bool operator!= (int a, const Rational& b) { return !(b==a); }

inline bool operator< (const Rational& a, const Integer& b) { return numerator(a) < b*denominator(a); }
inline bool operator< (const Rational& a, long b) { return numerator(a) < b*denominator(a); }
inline bool operator< (const Rational& a, int b) { return a < long(b); }

inline bool operator> (const Rational& a, const Integer& b) { return numerator(a) > b*denominator(a); }
inline bool operator> (const Rational& a, long b) { return numerator(a) > b*denominator(a); }
inline bool operator> (const Rational& a, int b) { return a > long(b); }

inline bool operator< (const Integer& a, const Rational& b) { return b>a; }
inline bool operator< (long a, const Rational& b) { return b>a; }
inline bool operator< (int a, const Rational& b) { return b>long(a); }
inline bool operator> (const Integer& a, const Rational& b) { return b<a; }
inline bool operator> (long a, const Rational& b) { return b<a; }
inline bool operator> (int a, const Rational& b) { return b<long(a); }

inline bool operator<= (const Rational& a, const Integer& b) { return !(a>b); }
inline bool operator<= (const Integer& a, const Rational& b) { return !(b<a); }

inline bool operator<= (const Rational& a, long b) { return !(a>b); }
inline bool operator<= (const Rational& a, int b) { return !(a>long(b)); }
inline bool operator<= (long a, const Rational& b) { return !(b<a); }
inline bool operator<= (int a, const Rational& b) { return !(b<long(a)); }

inline bool operator>= (const Rational& a, const Integer& b) { return !(a<b); }
inline bool operator>= (const Rational& a, long b) { return !(a<b); }
inline bool operator>= (const Rational& a, int b) { return !(a<long(b)); }
inline bool operator>= (const Integer& a, const Rational& b) { return !(b>a); }
inline bool operator>= (long a, const Rational& b) { return !(b>a); }
inline bool operator>= (int a, const Rational& b) { return !(b>long(a)); }

inline Integer& Integer::operator= (const Rational& r)
{
   if (! mpz_cmp_ui(mpq_denref(r.rep),1))
      mpz_set(rep,mpq_numref(r.rep));
   else
      mpz_tdiv_q(rep, mpq_numref(r.rep), mpq_denref(r.rep));
   return *this;
}

inline Integer::Integer(const Rational& r)
{
   if (! mpz_cmp_ui(mpq_denref(r.rep),1)) {
      mpz_init_set(rep,mpq_numref(r.rep));
   } else {
      mpz_init(rep);
      mpz_tdiv_q(rep, mpq_numref(r.rep), mpq_denref(r.rep));
   }
}

inline Rational& Rational::operator*= (const Integer& b)
{
   if (!*this) return *this;
   if (b) {
      Integer g=gcd(denominator(*this),b);
      if (g==1) {
	 mpz_mul(mpq_numref(rep), mpq_numref(rep), b.rep);
      } else {
	 mpz_divexact(mpq_denref(rep), mpq_denref(rep), g.rep);
	 mpz_divexact(g.rep, b.rep, g.rep);
	 mpz_mul(mpq_numref(rep), mpq_numref(rep), g.rep);
      }
   } else {
      *this=0;
   }
   return *this;
}

inline Rational& Rational::operator/= (const Integer& b)
{
   if (!b) zero_division();
   if (*this) {
      Integer g=gcd(numerator(*this),b);
      if (g==1) {
	 mpz_mul(mpq_denref(rep), mpq_denref(rep), b.rep);
      } else {
	 mpz_divexact(mpq_numref(rep), mpq_numref(rep), g.rep);
	 mpz_divexact(g.rep, b.rep, g.rep);
	 mpz_mul(mpq_denref(rep), mpq_denref(rep), g.rep);
      }
   }
   return *this;
}

inline Rational operator+ (const Rational& a, int b) { return a+long(b); }
inline Rational operator+ (const Integer& a, const Rational& b) { return b+a; }
inline Rational operator+ (long a, const Rational& b) { return b+a; }
inline Rational operator+ (int a, const Rational& b) { return b+long(a); }

inline Rational operator- (const Rational& a, int b) { return a-long(b); }
inline Rational operator- (int a, const Rational& b) { return long(a)-b; }

inline Rational operator* (const Rational& a, const Integer& b)
{
   if (!a || !b) return Rational();
   const Integer g=gcd(denominator(a),b);
   if (g==1)
      return Rational(mpz_mul, mpq_numref(a.rep), b.get_rep(), mpq_denref(a.rep));
   const Integer t=div_exact(b,g);
   return Rational(mpz_mul, mpq_numref(a.rep), t.get_rep(),
		   mpz_divexact, mpq_denref(a.rep), g.get_rep());
}

inline Rational operator* (const Rational& a, int b) { return a*long(b); }
inline Rational operator* (const Integer& a, const Rational& b) { return b*a; }
inline Rational operator* (long a, const Rational& b) { return b*a; }
inline Rational operator* (int a, const Rational& b) { return b*long(a); }

inline Rational operator/ (const Rational& a, int b) { return a/long(b); }
inline Rational operator/ (int a, const Rational& b) { return long(a)/b; }

inline Rational operator<< (const Rational& a, unsigned int k)
{
   return a << static_cast<unsigned long>(k);
}

inline Rational operator>> (const Rational& a, unsigned int k)
{
   return a >> static_cast<unsigned long>(k);
}

inline Rational operator<< (const Rational& a, long k)
{
   if (k<0) return a >> static_cast<unsigned long>(-k);
   return a << static_cast<unsigned long>(k);
}

inline Rational operator>> (const Rational& a, long k)
{
   if (k<0) return a << static_cast<unsigned long>(-k);
   return a >> static_cast<unsigned long>(k);
}

inline Rational operator<< (const Rational& a, int k) { return a << long(k); }
inline Rational operator>> (const Rational& a, int k) { return a >> long(k); }

template <typename Traits>
std::basic_ostream<char, Traits>& operator<< (std::basic_ostream<char, Traits>& os, const Rational& a)
{
   bool show_den=false;
   int s=numerator(a).strsize(os.flags());
   if (denominator(a)!=1) {
      show_den=true;
      s+=1+denominator(a).strsize(os.flags());
   }
   Integer::little_buffer buf(s);
   numerator(a).putstr(os.flags(), buf);
   if (show_den) {
      char *den_buf=buf+strlen(buf);
      *den_buf++ = '/';
      denominator(a).putstr(os.flags(), den_buf);
   }
   return os << static_cast<char*>(buf);
}

template <typename Traits>
std::basic_istream<char, Traits>& operator>> (std::basic_istream<char, Traits>& is, Rational& a)
{
   a.read(is);
   return is;
}

namespace std {

inline void swap(::Rational& r1, ::Rational& r2) { r1.swap(r2); }

inline Rational abs(const Rational& a)
{
   return Rational(mpz_abs, mpq_numref(a.rep), mpq_denref(a.rep));
}

template <>
struct numeric_limits<Rational> : numeric_limits<Integer> {
   static const bool is_integer=false;
};

} // end namespace std

namespace std_ext {

template <> struct hash<MP_RAT> : hash<MP_INT> {
protected:
   size_t _do(mpq_srcptr a) const
   {
      return hash<MP_INT>::_do(mpq_numref(a)) - hash<MP_INT>::_do(mpq_denref(a));
   }
public:
   size_t operator() (const MP_RAT& a) const { return _do(&a); }
};

template <> struct hash<Rational> : hash<MP_RAT>
{
   size_t operator() (const Rational& a) const { return _do(a.get_rep()); }
};

} // end namespace std_ext

#endif // _POLYMAKE_GMP_RATIONAL_H

// Local Variables:
// mode:C++
// End:
