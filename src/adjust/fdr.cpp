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
 
#include <algorithm>
#include <sstream>
#include <iostream>

#include "fdr.hpp"

#include "pvallocus.hpp"

using std::vector;

static bool 
lt_locus_pval(const PvalLocus &r1, const PvalLocus &r2) {
  return r1.combined_pval < r2.combined_pval;
}

static bool
ls_locus_position(const PvalLocus &r1, const PvalLocus &r2) {
  if (r1.chrom_ind < r2.chrom_ind)
    return true;

  if (r1.chrom_ind == r2.chrom_ind && r1.pos < r2.pos)
    return true;

  return false;
}

void
fdr(vector<PvalLocus> &loci) {
        
      std::sort(loci.begin(), loci.end(), lt_locus_pval);
    
      for (size_t ind = 0; ind < loci.size(); ++ind) {
        const double current_score = loci[ind].combined_pval;
        
        //Save current score.
        //std::stringstream ss;
        //ss << loci_iterators[ind]->name() << ":" << current_score;
        //loci_iterators[ind]->set_name(ss.str());
        
        //Assign a new one.
        const double corrected_pval = loci.size()*current_score/(ind + 1);
        loci[ind].corrected_pval = corrected_pval;
      }
    
      for (vector<PvalLocus>::reverse_iterator 
            it = loci.rbegin() + 1; it != loci.rend(); ++it) {
        
        const PvalLocus &prev_locus = *(it - 1);
        PvalLocus &cur_locus = *(it);
        
        cur_locus.corrected_pval =
              std::min(prev_locus.corrected_pval, cur_locus.corrected_pval);
      }
    
      for (vector<PvalLocus>::iterator it = loci.begin(); 
            it != loci.end(); ++it) {
        PvalLocus &cur_locus = *(it);
        if (cur_locus.corrected_pval > 1.0)
          cur_locus.corrected_pval = 1.0;
      }

      // Restore original order
      std::sort(loci.begin(), loci.end(), ls_locus_position);
}
