/* Copyright (c) 1997-2006
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

#ident "$Project: polymake $$Id: Integer.cc 7315 2006-04-02 21:37:53Z gawrilow $"

#include "Integer.h"

// fragments borrowed from read_int() of GNU libio
size_t Integer::strsize(const std::ios::fmtflags flags) const
{
   size_t s=1+(mpz_sgn(rep)<0);	// terminating '\0' and possible sign
   int base;
   switch (flags & (std::ios::basefield | std::ios::showbase)) {
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
   return s+mpz_sizeinbase(rep, base);
}

void Integer::putstr(std::ios::fmtflags flags, char* buf) const
{
   int base;
   switch (flags & (std::ios::basefield | std::ios::showbase)) {
   case int(std::ios::hex) | int(std::ios::showbase):
      mpz_get_str(buf+2, 16, rep);
      if (mpz_sgn(rep)<0) *buf++='-';
      *buf++='0';
      *buf='x';
      return;
   case int(std::ios::oct) | int(std::ios::showbase):
      mpz_get_str(buf+1, 8, rep);
      if (mpz_sgn(rep)<0) *buf++='-';
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
   mpz_get_str(buf, base, rep);
}

std::string Integer::to_string(int base) const
{
   std::string s(mpz_sizeinbase(rep,base)+2, '\0');
   mpz_get_str(const_cast<char*>(s.data()), base, rep);	// a nasty cast, I know...
   s.resize(s.find('\0'));
   return s;
}

char Integer::little_buffer::little[256];
