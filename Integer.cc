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

#ident "$Project: polymake $$Id: Integer.cc,v 1.16 2002/12/11 14:08:35 gawrilow Exp $"

#if defined(__GNUC__) && __GNUC_MINOR__<3 && !defined(__APPLE__)
#pragma implementation
#endif

#include <cctype>
#include "Integer.h"
#ifdef __APPLE__
# include <cstdlib>
#else
# include <alloca.h>
#endif

// fragments borrowed from read_int() of GNU libio
void Integer::read(std::istream& is, mpz_ptr dst, bool allow_sign) {
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
      return;
   }

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
      mpz_set_str(dst, text.c_str(), base);
   } else {
      is.setstate(std::ios::failbit);
   }
}

size_t Integer::strsize(const std::ostream& os, mpz_srcptr src) {
   size_t s=1+(mpz_sgn(src)<0);	// terminating '\0' and possible sign
   int base;
   switch (os.flags() & (std::ios::basefield | std::ios::showbase)) {
   case int(std::ios::hex) | int(std::ios::showbase):
      s+=2;
   case std::ios::hex:
      base=16;
      break;
   case int(std::ios::oct) | int(std::ios::showbase):
      s+=1;
   case std::ios::oct:
      base=8;
      break;
   default:
      base=10;
   }
   return s+mpz_sizeinbase(src, base);
}

void Integer::putstr(const std::ostream& os, char* buf, mpz_srcptr src) {
   int base;
   switch (os.flags() & (std::ios::basefield | std::ios::showbase)) {
   case int(std::ios::hex) | int(std::ios::showbase):
      mpz_get_str(buf+2, 16, src);
      if (mpz_sgn(src)<0) *buf++='-';
      *buf++='0';
      *buf='x';
      return;
   case int(std::ios::oct) | int(std::ios::showbase):
      mpz_get_str(buf+1, 8, src);
      if (mpz_sgn(src)<0) *buf++='-';
      *buf='0';
      return;
   case std::ios::hex:
      base=16;
      break;
   case std::ios::oct:
      base=8;
      break;
   default:
      base=10;
   }
   mpz_get_str(buf, base, src);
}

std::string Integer::to_string(int base) const {
   std::string s(mpz_sizeinbase(rep,base)+2, '\0');
   mpz_get_str(const_cast<char*>(s.data()),base,rep);	// a nasty cast, I know...
   s.resize(s.find('\0'));
   return s;
}

std::ostream& operator<< (std::ostream& os, const Integer& a) {
   int s=Integer::strsize(os,a.rep);
   char *buf=(char*)alloca(s);
   Integer::putstr(os,buf,a.rep);
   return os << buf;
}
