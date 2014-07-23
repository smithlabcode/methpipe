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
 
#include <sstream>

#include "bin_for_distance.hpp"


BinForDistance::BinForDistance(std::string spec_string) {
  std::replace(spec_string.begin(), spec_string.end(), ':', ' ');
    
  std::stringstream ss(spec_string);
  ss >> min_dist_ >> max_dist_ >> bin_size_;
    
  num_bins_ = (max_dist_ - min_dist_) / bin_size_;
  invalid_bin_ = num_bins_ + 1;
}
      
size_t BinForDistance::which_bin(size_t value) const {
  if (value < min_dist_)
    return invalid_bin_;
      
  const size_t bin = (value - min_dist_) / bin_size_;
    
  //Bin numbering is 0 based.
  if (bin >= num_bins_)
    return invalid_bin_;
      
  return bin;
}
