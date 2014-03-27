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
