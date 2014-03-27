
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>

#include "proximal_loci.hpp"

using std::copy;  using std::vector;
using std::back_inserter; using std::cout;
using std::endl;

bool ProximalLoci::get(vector<LocusIterator> &neighbors) {
  
  if (next_pos_ == loci_.end())
    return false;
  
  vector<LocusIterator>::const_iterator cur_pos = next_pos_;
  
  neighbors.clear();
  
  neighbors.push_back(*cur_pos);
  
  if ( cur_pos != loci_.begin() ) {
  
    vector<LocusIterator>::const_iterator up_pos = cur_pos;
  
    bool too_far = false;
  
    do {
      --up_pos;
    
      size_t up_dist = (*cur_pos)->begin() - (*up_pos)->end();
      if(up_dist <= max_distance_ 
          && (*cur_pos)->chrom() == (*up_pos)->chrom()) {
          neighbors.push_back(*up_pos);
      } else
        too_far = true;
          
    } while (!too_far && up_pos != loci_.begin());
  }
  
  std::reverse(neighbors.begin(), neighbors.end());
    
  if (cur_pos != loci_.end() - 1) {
   
    bool too_far = false;
    
    vector<LocusIterator>::const_iterator down_pos = cur_pos;
 
    do {
      
      ++down_pos;
      
      size_t down_dist = (*down_pos)->begin() - (*cur_pos)->end();
            
      if(down_dist <= max_distance_ 
          && (*down_pos)->chrom() == (*cur_pos)->chrom()) {
          neighbors.push_back(*down_pos);
      } else
        too_far = true;
      
    } while (!too_far && down_pos != loci_.end() - 1);
  }
  
  ++next_pos_;
  return true;
}
