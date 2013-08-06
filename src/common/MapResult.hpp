/*
 *    Part of RMAP software
 *
 *    Copyright (C) 2009 University of Southern California and
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

#ifndef MAP_RESULT_HPP
#define MAP_RESULT_HPP

#include <numeric>
#include <string>
#include <limits>
#include <vector>

struct MapResult {
  MapResult(size_t ste = std::numeric_limits<size_t>::max(),
	    size_t chr = 0,
	    bool str = true) : site(ste), chrom(chr), strand(str) {}
  void set(size_t ste, size_t chr, bool str);
  unsigned site   : 32;
  unsigned chrom  : 21;
  unsigned strand : 1;
  bool operator<(const MapResult& rhs) const {
    return (chrom < rhs.chrom || (chrom == rhs.chrom && site < rhs.site));
  }
  bool operator==(const MapResult& rhs) const {
    return site == rhs.site && chrom == rhs.chrom;
  }
  std::string tostring() const {
    return smithlab::toa(site) + "\t" + smithlab::toa(chrom) + "\t"
      + smithlab::toa(strand);
  }
  
  static std::vector<size_t> chrom_sizes;
  static std::vector<std::string> chrom_names;
};

std::vector<size_t> MapResult::chrom_sizes;
std::vector<std::string> MapResult::chrom_names;

class MultiMapResult {
public:
  bool empty() const {return mr.empty();}
  void sort() {std::sort(mr.begin(), mr.end());}
  void clear() {std::vector<MapResult>().swap(mr);}
  void swap(MultiMapResult &rhs) {
    mr.swap(rhs.mr);
    std::swap(score, rhs.score);
  }
  void add(size_t scr, size_t chr, size_t ste, bool str) {
    const size_t mr_size = mr.size();
    if (mr_size == 0 || scr < score) {
      mr.resize(0);
      mr.push_back(MapResult(ste, chr, str));
      score = scr;
    }
    else if (scr == score && mr_size < twice_max_count_lim)
      mr.push_back(MapResult(ste, chr, str));
  }
  size_t score;
  std::vector<MapResult> mr;
  bool ambiguous() const {return mr.size() > max_count;}
  void collapse() {
    std::sort(mr.begin(), mr.end());
    mr.erase(std::unique(mr.begin(), mr.end()), mr.end());
    if (mr.size() > max_count + 1)
      mr.erase(mr.begin() + max_count + 1, mr.end());
  }
  static size_t max_count;
  static size_t twice_max_count_lim;
  static void init(const size_t mc) {
    max_count = mc;
    twice_max_count_lim = 2*mc + 1;
  }
};

size_t MultiMapResult::max_count;
size_t MultiMapResult::twice_max_count_lim;

#endif
