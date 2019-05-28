/*    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith and Fang Fang
 *
 *    Authors: Fang Fang and Andrew D. Smith
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

#include <limits>
#include <string>
#include <stdexcept>

#include "Epiread.hpp"

using std::vector;
using std::string;

size_t
adjust_read_offsets(vector<epiread> &reads) {
  size_t first_read_offset = std::numeric_limits<size_t>::max();
  for (size_t i = 0; i < reads.size(); ++i)
    first_read_offset = std::min(reads[i].pos, first_read_offset);
  for (size_t i = 0; i < reads.size(); ++i)
    reads[i].pos -= first_read_offset;
  return first_read_offset;
}


size_t
get_n_cpgs(const vector<epiread> &reads) {
  size_t n_cpgs = 0;
  for (size_t i = 0; i < reads.size(); ++i)
    n_cpgs = std::max(n_cpgs, reads[i].end());
  return n_cpgs;
}

std::istream&
operator>>(std::istream &in, epiread &er) {
  string buffer;
  if (getline(in, buffer)) {
    std::istringstream is(buffer);
    if (!(is >> er.chr >> er.pos >> er.seq))
      throw std::runtime_error("malformed epiread line:\n" + buffer);
  }
  return in;
}


std::ostream&
operator<<(std::ostream &out, const epiread &er) {
  return out << er.chr << '\t' << er.pos << '\t' << er.seq;
}
