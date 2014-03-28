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

#ifndef BIN_FOR_DISTANCE_HPP_
#define BIN_FOR_DISTANCE_HPP_

#include <iostream>
#include <vector>
#include <algorithm>

class BinForDistance {
public:
  BinForDistance(std::string spec_string);
      
  size_t which_bin(size_t value) const;
  
  size_t num_bins() const {return num_bins_;}
  size_t invalid_bin() const {return invalid_bin_;}
  size_t max_dist() const {return max_dist_;}
  
private:
  size_t min_dist_;
  size_t max_dist_;
  size_t bin_size_;
  size_t num_bins_;
  size_t invalid_bin_;
};

#endif //BIN_FOR_DISTANCE_
