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

#ident "$Project: polymake $$Id: gmp_integer.cc,v 1.5 1999/02/08 14:19:59 gawrilow Exp $"

#pragma implementation
#include "gmp_integer.h"

// fragments borrowed from read_int() of GNU libio
void read(istream& is, mpz_t dst, char delim, bool allow_sign)
{
   string text;
   char c;

   if (allow_sign) {
      is >> c;		// skip leading whitespaces, consume the sign or the first digit
      switch (c) {
      case '-':
	 text += c;	// NOBREAK
      case '+':
	 is >> c;	// there may be some spaces after the sign, too
      }
   } else {
      is.get(c);
   }
   if (is.eof()) {
      is.set(ios::failbit);
      return;
   }

   int base=0;		// per default mpz_set_str must detect the numeric base automatically
   switch (is.flags() & ios::basefield) {
   case ios::hex:
      base=16; break;
   case ios::oct:
      base=8;  break;
   case ios::dec:
      base=10; break;
   }

   // gather all characters up to the first whitespace, delimiter, or end of file
   int cnt=0;
   while (!(is.eof() || string::traits_type::is_del(c) || c==delim)) {
      text += c;
      is.get(c);
      ++cnt;
   }
   if (!is.eof()) is.unget();

   if (!cnt || mpz_set_str(dst, text.c_str(), base)) is.set(ios::failbit);
}

ostream& operator<< (ostream& os, const Integer& a)
{
   int base=10;
   const char *base_prefix="";
   switch (os.flags() & ios::basefield) {
   case ios::hex:
      base=16;
      base_prefix="0x";
      break;
   case ios::oct:
      base=8;
      base_prefix="0";
      break;
   }

   const unsigned size=mpz_sizeinbase(a.rep, base)+2;
   char text[size], *t=text;
   mpz_get_str(text, base, a.rep);

   if (os.flags() & ios::showbase) {
      if (*t != '0') {
	 if (*t == '-') {
	    os << '-';
	    ++t;
	 }
	 os << base_prefix;
      }
   }
   return os << t;
}

gmp_alloc_init::gmp_alloc_init()
{
   mp_set_memory_functions (alloc::allocate, alloc::reallocate, alloc::deallocate);
}
