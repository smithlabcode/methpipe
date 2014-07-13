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

#ifndef EPIREAD
#define EPIREAD

#include <string>
#include <vector>
#include "smithlab_utils.hpp"

struct epiread {
  std::string chr;
  size_t pos;
  std::string seq;
  epiread() {}
  epiread(const size_t p, const std::string &s) : pos(p), seq(s) {}
  epiread(const std::string &c, const size_t p, const std::string &s)
      : chr(c), pos(p), seq(s) {}

  bool operator<(const epiread &other) const {
    return (chr < other.chr || (chr == other.chr && pos < other.pos));
  }
  size_t end() const {return pos + seq.length();}
  size_t length() const {return seq.length();}
};

std::istream& operator>>(std::istream &in, epiread &er);
std::ostream& operator<<(std::ostream &out, const epiread &er);

size_t
adjust_read_offsets(std::vector<epiread> &reads);

size_t
get_n_cpgs(const std::vector<epiread> &reads);

#endif
