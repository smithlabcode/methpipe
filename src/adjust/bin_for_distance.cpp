
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
