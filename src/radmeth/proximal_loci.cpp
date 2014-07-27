/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
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
 */
 
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>

#include "proximal_loci.hpp"

using std::copy;  using std::vector;
using std::back_inserter; using std::cout;
using std::endl;

bool ProximalLoci::get(vector<PvalLocus> &neighbors) {
  
  if (next_pos_ == loci_.end())
    return false;
  
  vector<PvalLocus>::const_iterator cur_pos = next_pos_;
  
  neighbors.clear();
  
  neighbors.push_back(*cur_pos);
  
  if ( cur_pos != loci_.begin() ) {
  
    vector<PvalLocus>::const_iterator up_pos = cur_pos;
  
    bool too_far = false;
  
    do {
      --up_pos;
    
      size_t up_dist = cur_pos->pos - (up_pos->pos + 1);
      if(up_dist <= max_distance_ && cur_pos->chrom_ind == up_pos->chrom_ind) {
          neighbors.push_back(*up_pos);
      } else
        too_far = true;
          
    } while (!too_far && up_pos != loci_.begin());
  }
  
  std::reverse(neighbors.begin(), neighbors.end());
    
  if (cur_pos != loci_.end() - 1) {
   
    bool too_far = false;
    
    vector<PvalLocus>::const_iterator down_pos = cur_pos;
 
    do {
      
      ++down_pos;
      
      size_t down_dist = down_pos->pos - (cur_pos->pos + 1);
            
      if( down_dist <= max_distance_ 
          && down_pos->chrom_ind == cur_pos->chrom_ind) {
          neighbors.push_back(*down_pos);
      } else
        too_far = true;
      
    } while (!too_far && down_pos != loci_.end() - 1);
  }
  
  ++next_pos_;
  return true;
}
