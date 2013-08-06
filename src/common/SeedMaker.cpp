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

#include "SeedMaker.hpp"
#include "rmap_utils.hpp"
#include "smithlab_utils.hpp"

#include <cmath>

using std::string;
using std::vector;

const size_t SeedMaker::max_seed_part = 32ul;
////////////////////////////////////////////////////////////////////////
// MAIN FUNCTIONS TO GET THE SEED SETS

void
SeedMaker::get_seed_set(const size_t read_width, const size_t n_seeds, 
			const size_t seed_weight, const size_t max_depth,
			vector<size_t> &all_profs) {
  
  size_t depth = 0;
  for (size_t i = 1; i <= seed_weight && depth < max_depth; i *= 2) {
    ++depth;
    const size_t shift = static_cast<size_t>(std::ceil((1.0*read_width/n_seeds)/i));
    const SeedMaker sm(read_width, n_seeds, seed_weight,
		       seed_weight/i, shift,
		       (read_width - seed_weight)/i);
    vector<size_t> profs;
    sm.get_seed_profiles(profs);
    
    if (all_profs.empty())
      all_profs.swap(profs);
    else {
      vector<size_t> tmp_profs;
      for (size_t j = 0; j < all_profs.size(); ++j)
	for (size_t i = 0; i < profs.size(); ++i)
	  tmp_profs.push_back((profs[i] | all_profs[j]));
      all_profs.swap(tmp_profs);
    }
  }
}


void
SeedMaker::first_last_seeds(const size_t read_width, const size_t n_seeds, 
			    const size_t seed_weight, vector<size_t> &profs) {
  
  const size_t shift = static_cast<size_t>(
      std::ceil(static_cast<double>(read_width - seed_weight)/(n_seeds - 1)));
  const SeedMaker sm_first(read_width, n_seeds,  seed_weight, seed_weight,
			   shift, read_width);
  const SeedMaker sm_last(read_width, n_seeds, seed_weight, 1, 1, n_seeds);
  
  vector<size_t> first_profs, last_profs;
  sm_first.get_seed_profiles(first_profs);
  sm_last.get_seed_profiles(last_profs);

  for (size_t i = 0; i < first_profs.size(); ++i)
    for (size_t j = 0; j < last_profs.size(); ++j)
      profs.push_back((first_profs[i] | last_profs[j]));

  vector<size_t> non_redundant_profs;
  const size_t mask = (2ul << (2*read_width - 1)) - 1;
  for (size_t i = 0; i < profs.size(); ++i) {
    bool found_containing = false;
    for (size_t j = 0; j < non_redundant_profs.size() && !found_containing; ++j)
      if ((non_redundant_profs[j] | profs[i]) == non_redundant_profs[j])
	found_containing = true;
    if (!found_containing)
      non_redundant_profs.push_back(profs[i] & mask);
  }
  profs.swap(non_redundant_profs);
}


void
SeedMaker::last_seeds(const size_t read_width, const size_t n_seeds, 
		      const size_t seed_weight, vector<size_t> &profs) {
  const SeedMaker sm_last(read_width, n_seeds, seed_weight, 1, 1, n_seeds);
  sm_last.get_seed_profiles(profs);
}


void
SeedMaker::first_seeds(const size_t read_width, const size_t n_seeds, 
		       const size_t seed_weight, vector<size_t> &profs) {
  const size_t shift = static_cast<size_t>(
      std::ceil(static_cast<double>(read_width - seed_weight)/(n_seeds - 1)));
  const SeedMaker sm_first(read_width, n_seeds,  seed_weight, seed_weight,
			   shift, read_width);
  sm_first.get_seed_profiles(profs);
}


void
SeedMaker::last_two_seeds(const size_t read_width, const size_t n_seeds, 
			  const size_t seed_weight, vector<size_t> &profs) {
  
  const SeedMaker sm_last(read_width, n_seeds, seed_weight, 1, 1, n_seeds);
  const SeedMaker sm_first(read_width, n_seeds, seed_weight, 2, 2, 2*n_seeds);
  
  vector<size_t> first_profs, last_profs;
  sm_first.get_seed_profiles(first_profs);
  sm_last.get_seed_profiles(last_profs);
  
  for (size_t i = 0; i < first_profs.size(); ++i)
    for (size_t j = 0; j < last_profs.size(); ++j)
      profs.push_back((first_profs[i] | last_profs[j]));

  vector<size_t> non_redundant_profs;
  for (size_t i = 0; i < profs.size(); ++i) {
    bool found_containing = false;
    for (size_t j = 0; j < non_redundant_profs.size() && !found_containing; ++j)
      if ((non_redundant_profs[j] | profs[i]) == non_redundant_profs[j])
	found_containing = true;
    if (!found_containing)
      non_redundant_profs.push_back(profs[i]);
  }
  profs.swap(non_redundant_profs);
}


////////////////////////////////////////////////////////////////////////
// MAIN FUNCTIONS TO GET THE SEED SETS

SeedMaker::SeedMaker(const size_t r, const size_t s, const size_t w, 
		     const size_t b, const size_t f, const size_t g) :
  read_width(r), n_seeds(s), seed_weight(w), block_size(b), 
  shift_size(f), gap_size(g) {
  seed_mask = 0;
  const size_t n_blocks = static_cast<size_t>(std::ceil(1.0*seed_weight/block_size));
  for (size_t i = 0; i < n_blocks; ++i) {
    for (size_t j = 0; j < gap_size - block_size; ++j)
      seed_mask = (seed_mask << 2);
    for (size_t j = 0; j < block_size; ++j)
      seed_mask = (seed_mask << 2) + 3;
  }
  size_t full_mask = 0;
  for (size_t i = 0; i < 2*read_width; ++i) {
    full_mask <<= 1;
    full_mask |= 1ul;
  }
  seed_mask &= full_mask;
}

void
SeedMaker::get_seed_profiles(vector<size_t> &profiles) const {
  for (size_t i = 0; i < n_seeds; ++i)
    profiles.push_back((seed_mask << (2*i*shift_size)));
}

string
SeedMaker::tostring() const {
  std::ostringstream ss;
  ss << "READ WIDTH:  " << read_width  << "\n" 
     << "N SEEDS:     " << n_seeds     << "\n" 
     << "SEED INC:    " << block_size  << "\n" 
     << "SEED WEIGHT: " << seed_weight << "\n";
  vector<size_t> profiles;
  get_seed_profiles(profiles);
  for (size_t i = 0; i < profiles.size(); ++i)
    ss << bits2string_masked(rmap_bits::all_ones, profiles[i]) << "\n";
  return ss.str();
}

////////////////////////////////////////////////////////////////////////
// STATIC

size_t
SeedMaker::make_read_word(const string &s) {
  size_t multiplier = 1, index = 0;
  size_t n = s.length();
  do { 
    --n;
    index += base2int(s[n])*multiplier;
    multiplier *= smithlab::alphabet_size; 
  } while (n > 0);
  return index;
}


size_t
SeedMaker::apply_seed(const size_t shift, const size_t sd, const size_t rw) {
  return ((sd & rw) >> 2*shift) | ((sd & rw) << (64ul - 2*shift));
}


bool
SeedMaker::valid_seed(const size_t the_seed, const string &read) {
  return true;
}

size_t
SeedMaker::get_heavy_shift(const size_t width, const size_t sd) {
  size_t best_count = 0, best_shift = 0;
  for (size_t i = 0; i < 32 - width; ++i) {
    size_t count = 0;
    for (size_t j = 0; j < width; ++j)
      count += (((3ul << 2*(i + j)) & sd) != 0);
    if(count > best_count) {
      best_count = count;
      best_shift = i;
    }
  }
  return best_shift;
}
