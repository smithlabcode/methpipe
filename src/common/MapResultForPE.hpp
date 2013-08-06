/*
 *    Part of RMAP software
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

#ifndef MAP_RESULT_HPP
#define MAP_RESULT_HPP

#include <numeric>
#include <string>
#include <limits>
#include <vector>


struct MapResult {
  MapResult(size_t ste = std::numeric_limits<size_t>::max(),
	    size_t chr = 0,
	    size_t scr = 1023,
	    bool str = true) : site(ste), chrom(chr), score(scr), strand(str) {}
  unsigned site   : 32;
  unsigned chrom  : 21;
  unsigned score  : 10;
  unsigned strand : 1;
  bool operator<(const MapResult& rhs) const {
    return (chrom < rhs.chrom || (chrom == rhs.chrom && site < rhs.site));
  }
  bool operator==(const MapResult& rhs) const {
    return site == rhs.site && chrom == rhs.chrom;
  }
  std::string tostring() const {
    return smithlab::toa(site) + "\t" + smithlab::toa(chrom) + "\t"
      + smithlab::toa(strand) + "\t" + smithlab::toa(score);
  }
  
  static std::vector<size_t> chrom_sizes;
  static std::vector<std::string> chrom_names;
};

std::vector<size_t> MapResult::chrom_sizes;
std::vector<std::string> MapResult::chrom_names;


struct mr_score_less {
  bool operator()(const MapResult &x, const MapResult &y) const 
  {return x.score < y.score;}
};


class MultiMapResult {
public:
  bool empty() const {return mr.empty();}
  void sort() {std::sort(mr.begin(), mr.end());}
  void clear() {std::vector<MapResult>().swap(mr);}
  void add(size_t scr, size_t chr, size_t ste, bool str) {
    if (mr.size() < twice_max_count || scr < mr.back().score) {
      mr.push_back(MapResult(ste, chr, scr, str)); 
      push_heap(mr.begin(), mr.end(), mr_sl);
      if (mr.size() > twice_max_count) {
	pop_heap(mr.begin(), mr.end(), mr_sl);
	mr.pop_back();
      }
    }
  }
  std::vector<MapResult> mr;
  
  void collapse() {
    std::sort(mr.begin(), mr.end());
    mr.erase(std::unique(mr.begin(), mr.end()), mr.end());
    std::sort(mr.begin(), mr.end(), mr_sl);
    if (mr.size() > max_count)
      mr.erase(mr.begin() + max_count, mr.end());
    std::make_heap(mr.begin(), mr.end(), mr_sl);
  }
  static size_t max_count;
  static size_t twice_max_count;
  static mr_score_less mr_sl;
  static void init(const size_t mc) {
    max_count = mc;
    twice_max_count = 2*mc;
  }
};

size_t MultiMapResult::max_count;
size_t MultiMapResult::twice_max_count;
mr_score_less MultiMapResult::mr_sl;

void swap(MultiMapResult &a, MultiMapResult &b) {
  swap(a.mr, b.mr);
}

#endif
