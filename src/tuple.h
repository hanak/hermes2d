// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __HERMES2D_TUPLE_H
#define __HERMES2D_TUPLE_H

#include "common.h"

/// A vector of values. Used to pass a variable number of parameters type-safe and a little bit comfortably.
template<typename T>
class Tuple: public std::vector<T> {
public:
  //template<typename ForwardIterator>
  //Tuple(ForwardIterator first, ForwardIterator last) {
  //  reserve(last - first);
  //  while(first != last) {
  //    push_back(*first);
  //    first++;
  //  }
  //};
  Tuple(const T& a) { push_back(a); };
  Tuple(const T& a, const T& b) { push_back(a); push_back(b); };
  Tuple(const T& a, const T& b, const T& c) { push_back(a); push_back(b); push_back(c); };
  Tuple(const T& a, const T& b, const T& c, const T& d) { push_back(a); push_back(b); push_back(c); push_back(d); };
  Tuple(const T& a, const T& b, const T& c, const T& d, const T& e) { push_back(a); push_back(b); push_back(c); push_back(d); push_back(e); };
  Tuple(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f) { push_back(a); push_back(b); push_back(c); push_back(d); push_back(e); push_back(f); };
  Tuple(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g) { push_back(a); push_back(b); push_back(c); push_back(d); push_back(e); push_back(f); push_back(g); };
  Tuple(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h) { push_back(a); push_back(b); push_back(c); push_back(d); push_back(e); push_back(f); push_back(g); push_back(h); };
  Tuple(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h, const T& i) { push_back(a); push_back(b); push_back(c); push_back(d); push_back(e); push_back(f); push_back(g); push_back(h); push_back(i); };
  Tuple(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h, const T& i, const T& j) { push_back(a); push_back(b); push_back(c); push_back(d); push_back(e); push_back(f); push_back(g); push_back(h); push_back(i); push_back(j); };
};

#endif
