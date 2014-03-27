
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>

#include "locus.hpp"
#include "merge.hpp"

using std::vector; using std::string;
using std::istringstream; using std::ostringstream;
using std::stringstream; using std::ostream; 
using std::istream;

static void
extract_scores(const string &name, double &lod, double &mindiff, double &pval) {
  istringstream name_stream(name);
  string token;
  getline(name_stream, token, ':');
  
  getline(name_stream, token, ':');
  lod = atof(token.c_str());
  
  getline(name_stream, token, ':');
  mindiff = atof(token.c_str());
  
  getline(name_stream, token, ':');
  pval = atof(token.c_str());
}

static bool
read_next_significant_cpg(istream &cpg_stream, Locus &cpg, double cutoff, 
                          bool &skipped_any) {
  
  skipped_any = false;
  string cpg_encoding;
  
  while (getline(cpg_stream, cpg_encoding)) {
    Locus cur_cpg(cpg_encoding);
    if (0 <= cur_cpg.score() && cur_cpg.score() < cutoff) {
      cpg = cur_cpg;
      return true;
    }
    skipped_any = true;
  }
  
  return false;
}

static void
update_dmr_name(Locus &dmr, double lod_sum, double mindiff_sum) {
  ostringstream name_stream;
  name_stream << "dmr:" 
              << lod_sum / dmr.score() << ":" 
              << mindiff_sum / dmr.score();
  dmr.set_name(name_stream.str());
}

void
merge(istream &cpg_stream, ostream &dmr_stream, double cutoff) {
  
  bool skipped_last_cpg;
  Locus dmr;
  
  if (!read_next_significant_cpg(cpg_stream, dmr, cutoff, skipped_last_cpg))
    return;
  
  dmr.set_score(1);
  double lod_sum, mindiff_sum, unadjusted_pval;
  extract_scores(dmr.name(), lod_sum, mindiff_sum, unadjusted_pval);
  
  if (unadjusted_pval >= cutoff) {
    lod_sum = 0;
    mindiff_sum = 0;
  }
  
  Locus cpg;
  while(read_next_significant_cpg(cpg_stream, cpg, cutoff, skipped_last_cpg)) {
    double lod, mindiff;
    
    extract_scores(cpg.name(), lod, mindiff, unadjusted_pval);
    
    if (skipped_last_cpg || cpg.chrom() != dmr.chrom()) {
      update_dmr_name(dmr, lod_sum, mindiff_sum);
      
      if (lod_sum != 0) {
        dmr_stream << dmr << std::endl;
      }
      
      dmr = cpg;
      dmr.set_score(1);
      
      if (unadjusted_pval >= cutoff) {
        lod_sum = 0;
        mindiff_sum = 0;
      } else {
        lod_sum = lod;
        mindiff_sum = mindiff;
      }
      
    } else {
      dmr.set_end(cpg.end());
      dmr.set_score(dmr.score() + 1);
      
      if (unadjusted_pval < cutoff) {
        lod_sum += lod;
        mindiff_sum += mindiff;
      }
    }
    
  }
  
  update_dmr_name(dmr, lod_sum, mindiff_sum);
  dmr_stream << dmr << std::endl;
  
}
