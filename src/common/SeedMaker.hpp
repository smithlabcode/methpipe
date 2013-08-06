/*
 *    Part of RMAP software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
 *                       University of Southern California and
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

#ifndef SEED_MAKER_HPP
#define SEED_MAKER_HPP

#include <string>
#include <vector>

class SeedMaker {
public:
  SeedMaker(const size_t r, const size_t s, const size_t w, 
	    const size_t b, const size_t f, const size_t g);
  size_t operator()(size_t read_word, size_t idx) const;
  std::string tostring() const;
  void get_seed_profiles(std::vector<size_t> &profiles) const;
  
  static void
  get_seed_set(const size_t rw, const size_t ns, 
	       const size_t sw, const size_t md, std::vector<size_t> &p);
  static void
  first_last_seeds(const size_t rw, const size_t ns, 
 		   const size_t sw, std::vector<size_t> &p);
  static void
  last_seeds(const size_t read_width, const size_t n_seeds, 
	     const size_t seed_weight, std::vector<size_t> &profs);
  static void
  first_seeds(const size_t read_width, const size_t n_seeds, 
	      const size_t seed_weight, std::vector<size_t> &profs);
  static void
  last_two_seeds(const size_t read_width, const size_t n_seeds, 
		 const size_t seed_weight, std::vector<size_t> &profs);
  static size_t make_read_word(const std::string &s);
  static void update_read_word(const size_t base, size_t &key);
  static void update_bad_bases(const size_t base, size_t &bads);
  static bool valid_seed(const size_t the_seed, const std::string &read);

  static size_t apply_seed(const size_t shift, const size_t sd, const size_t rw);
  static size_t get_heavy_shift(const size_t width, const size_t sd);
  
  static const size_t max_seed_part; // = 32ul;
  
private:
  size_t read_width;
  size_t n_seeds;
  size_t seed_weight;
  size_t block_size;
  size_t shift_size;
  size_t gap_size;
  size_t seed_mask;
};


inline std::ostream& 
operator<<(std::ostream& s, const SeedMaker& maker) {
  return s << maker.tostring();}

inline void
SeedMaker::update_read_word(const size_t base, size_t &key) {
  key <<= static_cast<size_t>(2);
  key += base;
}

inline void
SeedMaker::update_bad_bases(const size_t base, size_t &bads) {
  bads = (bads << static_cast<size_t>(2)) + 
    static_cast<size_t>(3*static_cast<size_t>(base == 4));
}

#endif
