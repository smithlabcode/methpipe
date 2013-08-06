/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2012 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "clip_adaptor_from_reads.hpp"

using std::string;
using std::min;

// We are using a conservative approach to clip that adaptors in which
// the adaptor sequence is only required to match some at some initial
// portion, and then the rest of the read is not examined.

const size_t head_length = 14;
const size_t sufficient_head_match = 11;
const size_t min_overlap = 5;

size_t
similarity(const string &s, const size_t pos, const string &adaptor) {
  const size_t lim = min(min(s.length() - pos, adaptor.length()), head_length);
  size_t count = 0;
  for (size_t i = 0; i < lim; ++i)
    count += (s[pos + i] == adaptor[i]);
  return count;
}

size_t 
clip_adaptor_from_read(const string &adaptor, string &s) {
  size_t lim1 = s.length() - head_length + 1;
  for (size_t i = 0; i < lim1; ++i)
    if (similarity(s, i, adaptor) >= sufficient_head_match) {
      fill(s.begin() + i, s.end(), 'N');
      return s.length() - i;
    }
  const size_t lim2 = s.length() - min_overlap + 1;
  for (size_t i = lim1; i < lim2; ++i)
    if (similarity(s, i, adaptor) >= s.length() - i - 1) {
      fill(s.begin() + i, s.end(), 'N');
      return s.length() - i;
    }
  return 0;
}
