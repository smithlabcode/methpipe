
#include <algorithm>
#include <sstream>
#include <iostream>

#include "fdr.hpp"

#include "locus.hpp"

using std::vector;

static bool 
lt_locus_iterator_score(LocusIterator r1, LocusIterator r2) {
  return r1->score() < r2->score();
}

void
fdr(vector<LocusIterator> &loci_iterators) {
        
      std::sort(loci_iterators.begin(), loci_iterators.end(), 
                lt_locus_iterator_score);
    
      for (size_t ind = 0; ind < loci_iterators.size(); ++ind) {
        const double current_score = loci_iterators[ind]->score();
        
        //Save current score.
        std::stringstream ss;
        ss << loci_iterators[ind]->name() << ":" << current_score;
        loci_iterators[ind]->set_name(ss.str());
        
        //Assign a new one.
        const double new_score = loci_iterators.size()*current_score/(ind + 1);
        loci_iterators[ind]->set_score(new_score);
      }
    
      for (vector<LocusIterator>::reverse_iterator 
          it = loci_iterators.rbegin() + 1; it != loci_iterators.rend(); ++it) {
        
        const LocusIterator &prev_region_it = *(it - 1);
        LocusIterator &cur_region_it = *(it);
        
        cur_region_it->set_score(
              std::min(prev_region_it->score(), cur_region_it->score())
            );
      }
      
      for (vector<LocusIterator>::iterator it = loci_iterators.begin();
            it != loci_iterators.end(); ++it) {
        LocusIterator &cur_region_it = *(it);
        if (cur_region_it->score() > 1.0)
          cur_region_it->set_score(1.0);
      }
}
