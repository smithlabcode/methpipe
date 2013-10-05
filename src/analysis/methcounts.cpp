/*    methcounts: a program for counting the methylated and
 *    unmethylated reads mapping over each CpG or C
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith and Elena Harris and Song Qiang
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

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <utility>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "QualityScore.hpp"

#include "bsutils.hpp"
#include "FileIterator.hpp"
#include "MappedRead.hpp"

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::ofstream;

// GLOBGAL VARIABLE TO CONTROL OUTPUT FORMAT
static bool USE_ALT_OUTPUT = true;

struct Compare : public std::binary_function<
  std::pair<string, size_t>, std::pair<string, size_t>, bool> {
  bool operator()(const std::pair<string, size_t> &a,
                  const std::pair<string, size_t> &b) 
  {return (a < b);}
};

struct MethStat {

  MethStat() : 
    total_sites(0), total_covered(0),
    max_cov(0), sum_cov(0), sum_cov_Cs(0) {}
  
  string tostring() const;
  
  void collect(const size_t meth_count, const size_t total) {
    total_sites++;
    if (total > 0) {
      total_covered++;
      max_cov = max(max_cov, total);
      sum_cov += total;
      sum_cov_Cs += meth_count;
    }
  }
  
  size_t total_sites;
  size_t total_covered;
  size_t max_cov;
  size_t sum_cov;
  size_t sum_cov_Cs;
};


string
MethStat::tostring() const {
  std::ostringstream out;
  
  out << "SITES:\t" << total_sites << endl
      << "SITES COVERED:\t" << total_covered << endl
      << "FRACTION:\t" << static_cast<double>(total_covered)/total_sites << endl;
  
  const double overall_cov = 
    static_cast<double>(sum_cov)/max(static_cast<size_t>(1), total_sites);
  const double covered_cov = 
    static_cast<double>(sum_cov)/max(static_cast<size_t>(1), total_covered);
  out << "MAX COVERAGE:\t" << max_cov << endl
      << "MEAN COVERAGE:\t" << overall_cov << endl
      << "MEAN (WHEN > 0):\t" << covered_cov << endl;
  
  const double meth_level = 
    static_cast<double>(sum_cov_Cs)/max(static_cast<size_t>(1), sum_cov);
  out << "MEAN METHYLATION:\t" << meth_level;
  return out.str();
}


std::ostream& 
operator<<(std::ostream& the_stream, const MethStat& ms) {
  return the_stream << ms.tostring();
}


static string
cytosine_type_tag(const string &chrom_seq, size_t chrom_pos, char strand) {
  if (strand == '+') {
    const size_t chr_len = chrom_seq.length();
    const size_t next_pos = chrom_pos + 1;
    if (next_pos == chr_len) 
      return "CHH";
    if (is_guanine(chrom_seq[next_pos])) 
      return "CG";
    const size_t next_next_pos = next_pos + 1;
    if (next_next_pos == chr_len) return "CHH";
    return is_guanine(chrom_seq[next_next_pos]) ? "CHG" : "CHH";
  }
  else {
    if (chrom_pos == 0) 
      return "CHH";
    const size_t next_pos = chrom_pos - 1;
    if (is_cytosine(chrom_seq[next_pos])) 
      return "CG";
    if (next_pos == 0) 
      return "CHH";
    const size_t next_next_pos = next_pos - 1;
    return is_cytosine(chrom_seq[next_next_pos]) ? "CHG" : "CHH";
  }
}


/***********************************************************************
 * FUNCTIONS BELOW ARE FOR FASTQRecord OJBECTS AND *NOT* USING THE
 * QUALITY INFORAMTION
 */

static void
add_contribution_cpg(const size_t offset, const MappedRead &r,
		     size_t &meth, size_t &unmeth) {
  if (r.r.pos_strand() && 
      // (offset - r.r.get_start() < s.first.length()) &&
      (r.r.get_start() <= offset)) {
    const size_t position = offset - r.r.get_start();
    if(position >= r.seq.length())
      return;//throw SMITHLABException("ERROR: Reads must be sorted by chromosome and end position.");
    //    assert(position < r.seq.length());
    if (is_cytosine(r.seq[position])) ++meth;
    if (is_thymine(r.seq[position])) ++unmeth;
  }
  if (r.r.neg_strand() && 
      // (r.r.get_start() <= offset + 1) && 
      (offset - r.r.get_start() + 2 <= r.seq.length())) {
    const size_t position = (r.seq.length() - 1) - 
      ((offset + 1) - r.r.get_start());
    if(position >= r.seq.length())
      return;//throw SMITHLABException("ERROR: Reads must be sorted by chromosome and end position.");

//    assert(position < r.seq.length());
    if (is_cytosine(r.seq[position])) ++meth;
    if (is_thymine(r.seq[position])) ++unmeth;
  }
}

static void
add_contribution_c(const size_t offset, const MappedRead &r,
		   size_t &meth, size_t &unmeth) {
  if (r.r.pos_strand()) {
    const size_t position = offset - r.r.get_start();
    //    assert(position < r.seq.length());
    if (position >= r.seq.length())
      return;//throw SMITHLABException("ERROR: Reads must be sorted by chromosome and end position.");

    if (is_cytosine(r.seq[position])) ++meth;
    if (is_thymine(r.seq[position])) ++unmeth;
  }
}

static void
add_contribution_g(const size_t offset, const MappedRead &r,
		   size_t &meth, size_t &unmeth) {
  if (r.r.neg_strand()) {
    const size_t position = (r.seq.length() - 1) - (offset - r.r.get_start());
    //    assert(position < r.seq.length());
    if (position >= r.seq.length())
      return;//throw SMITHLABException("ERROR: Reads must be sorted by chromosome and end position.");
    
    if (is_cytosine(r.seq[position])) ++meth;
    if (is_thymine(r.seq[position])) ++unmeth;
  }
}


/***********************************************************************
 * FUNCTIONS BELOW ARE FOR FASTQRecord OJBECTS AND USING THE QUALITY
 * INFORAMTION
 */

struct QualityChecker {
  QualityChecker(const FASTQScoreType &score_format, const double c) {
    cutoff = (FASTQScoreIsPhred(score_format) ?
	      phred_to_quality_character(error_probability_to_phred(c)) :
	      solexa_to_quality_character(error_probability_to_solexa(c)));
  }
  bool operator()(const MappedRead &r, const size_t position) const {
    return r.scr[position] >= cutoff;
  }
  char cutoff;
};

static void
add_contribution_cpg(const QualityChecker &qc,
		     const size_t offset, const MappedRead &r,
		     size_t &meth, size_t &unmeth) {
  if (r.r.pos_strand() && (r.r.get_start() <= offset)) {
    const size_t position = offset - r.r.get_start();
    //    assert(position < r.seq.length());
    if(position >= r.seq.length())
      return;//throw SMITHLABException("ERROR: Reads must be sorted by chromosome and end position.");

    if (qc(r, position)) {
      if (is_cytosine(r.seq[position])) ++meth;
      if (is_thymine(r.seq[position])) ++unmeth;
    }
  }
  if (r.r.neg_strand() && (offset - r.r.get_start() + 2 <= r.seq.length())) {
    const size_t position = (r.seq.length() - 1) - 
      ((offset + 1) - r.r.get_start());
    //    assert(position < r.seq.length());
    if(position >= r.seq.length())
      return;//throw SMITHLABException("ERROR: Reads must be sorted by chromosome and end position.");

    if (qc(r, position)) {
      if (is_cytosine(r.seq[position])) ++meth;
      if (is_thymine(r.seq[position])) ++unmeth;
    }
  }
}

static void
add_contribution_c(const QualityChecker &qc,
		   const size_t offset, const MappedRead &r,
		   size_t &meth, size_t &unmeth) {
  if (r.r.pos_strand()) {
    const size_t position = offset - r.r.get_start();
    //    assert(position < r.seq.length());
    if(position >= r.seq.length())
      return;//throw SMITHLABException("ERROR: Reads must be sorted by chromosome and end position.");

    if (qc(r, position)) {
      if (is_cytosine(r.seq[position])) ++meth;
      if (is_thymine(r.seq[position])) ++unmeth;
    }
  }
}

static void
add_contribution_g(const QualityChecker &qc,
		   const size_t offset, const MappedRead &r,
		   size_t &meth, size_t &unmeth) {
  if (r.r.neg_strand()) {
    const size_t position = (r.seq.length() - 1) - (offset - r.r.get_start());
    //    assert(position < r.seq.length());
    if(position >= r.seq.length())
      return;//throw SMITHLABException("ERROR: Reads must be sorted by chromosome and end position.");

    if (qc(r, position)) {
      if (is_cytosine(r.seq[position])) ++meth;
      if (is_thymine(r.seq[position])) ++unmeth;
    }
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static bool
precedes(const MappedRead &r, const size_t offset) {
  return r.r.get_end() <= offset;
}


static bool
succeeds(const MappedRead &r, const size_t offset) {
  return r.r.get_start() > offset;
}


static void
advance(const size_t first, const size_t last,
	const GenomicRegion &chrom_region, 
	FileIterator<MappedRead> &regions,
        const size_t max_length) {

  while (regions.last_is_good() && 
	 chrom_region.same_chrom(regions.get_last()->r) &&
	 !succeeds(*regions.get_last(), last+max_length)) {
    regions.increment_last();
  }

  //   if (regions.last_is_good() != reads.last_is_good())
  //     throw SMITHLABException("read and map files seem out of sync");
  while (regions.first_is_good() && 
	 chrom_region.same_chrom(regions.get_first()->r) &&
	 precedes(*regions.get_first(), first)) {
    regions.increment_first();
  }

  //   if (regions.first_is_good() != reads.first_is_good())
  //     throw SMITHLABException("read and map files seem out of sync");
}


static void
scan_chromosome_cpg(const QualityChecker &qc,
		    const string &chrom,
		    const GenomicRegion &chrom_region,
		    const double max_mismatches,
		    FileIterator<MappedRead> &regions,
		    std::ostream &out, MethStat &meth_stat_collector,
                    const size_t max_length) {
  const string chrom_name(chrom_region.get_chrom());
  size_t i = 0;
  for ( i = 0; i < chrom.length() - 1 && regions.first_is_good(); ++i) {
    if (is_cpg(chrom, i)) {
      /* need the "+1" below because of the 'G' in CpG */
      advance(i, i + 1, chrom_region, regions, max_length);
      size_t meth_count = 0, unmeth_count = 0;
      for (vector<MappedRead>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j)
	if (j->r.get_score() <= max_mismatches)
	  add_contribution_cpg(qc, i, *j, meth_count, unmeth_count);
      const size_t total = meth_count + unmeth_count;
      meth_stat_collector.collect(meth_count, total);
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t+\tCpG\t" 
            << meth_count/max(1.0, static_cast<double>(total)) << "\t"
            << total << endl;
      else
        out << chrom_name << "\t" 
            << i << "\t" << i + 1 << "\tCpG:" << total << "\t" 
            << meth_count/max(1.0, static_cast<double>(total)) << "\t+\n";
    }//if
  }//for
  for(;  i < chrom.length() - 1; ++i){
    if (is_cpg(chrom, i)) {
      const size_t total = 0;
      meth_stat_collector.collect(total, total);
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t+\tCpG\t" 
            << total << "\t" << total << endl;
      else
        out << chrom_name << "\t" << i << "\t" << i + 1 << "\tCpG:"
            << total << "\t" << total << "\t+\n";
    }//if
  }//for

}//scan_chromosome


static void
scan_chromosome(const QualityChecker &qc,
		const string &chrom, const GenomicRegion &chrom_region,
		const double max_mismatches,
		FileIterator<MappedRead> &regions, 
		std::ostream &out, MethStat &meth_stat_collector,
                const size_t max_length) {

  const string chrom_name(chrom_region.get_chrom());
  size_t i = 0; 
  size_t chrom_len = chrom.length();
  for (i = 0; i < chrom_len - 1 && regions.first_is_good(); ++i) {
    advance(i, i, chrom_region, regions, max_length);
    if (is_cytosine(chrom[i])) { // && !is_guanine(chrom[i + 1])) {
      size_t meth_count = 0, unmeth_count = 0;
      for (vector<MappedRead>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j)
	if (j->r.get_score() <= max_mismatches)
	  add_contribution_c(qc, i, *j, meth_count, unmeth_count);
      const size_t total = meth_count + unmeth_count;
      meth_stat_collector.collect(meth_count, total);
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t+\t"  
            << cytosine_type_tag(chrom, i, '+') << "\t"
            << meth_count/max(1.0, static_cast<double>(total)) << "\t"
            << total << endl;
      else
        out << chrom_name << "\t" << i << "\t" << i + 1 << "\t" 
            << cytosine_type_tag(chrom, i, '+') << ":"
            << total << "\t" 
            << meth_count/max(1.0, static_cast<double>(total)) << "\t+\n";
    }
    if (is_guanine(chrom[i])) { // && !is_cytosine(chrom[i - 1])) {
      size_t meth_count = 0, unmeth_count = 0;
      for (vector<MappedRead>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j)
	if (j->r.get_score() <= max_mismatches)
	  add_contribution_g(qc, i, *j, meth_count, unmeth_count);
      const size_t total = meth_count + unmeth_count;
      meth_stat_collector.collect(meth_count, total);
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t-\t" 
            << cytosine_type_tag(chrom, i, '-') << "\t" 
            << meth_count/max(1.0, static_cast<double>(total)) << "\t"
            << total << endl;
      else
        out << chrom_name << "\t" << i << "\t" << i + 1 << "\t" 
            << cytosine_type_tag(chrom, i, '-') << ":" 
            << total << "\t" 
            << meth_count/max(1.0, static_cast<double>(total)) << "\t-\n";
    }
  }//for
  for (; i < chrom_len - 1 ; ++i){
    const size_t total = 0;
    if (is_cytosine(chrom[i])) {
      meth_stat_collector.collect(total, total);
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t+\t"  
            << cytosine_type_tag(chrom, i, '+') << "\t"
            << total << "\t" << total << endl;
      else
        out << chrom_name << "\t" << i << "\t" << i + 1 << "\t" 
            << cytosine_type_tag(chrom, i, '+') << ":"
            << total << "\t" << total << "\t+\n";
    }
    if (is_guanine(chrom[i])) {
      meth_stat_collector.collect(total, total);
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t+\t"  
            << cytosine_type_tag(chrom, i, '-') << "\t"
            << total << "\t" << total << endl;
      else
        out << chrom_name << "\t" << i << "\t" << i + 1 << "\t" 
            << cytosine_type_tag(chrom, i, '-') << ":"
            << total << "\t" << total << "\t+\n";
    }
  }//for
}


static void
scan_chromosome_cpg(const string &chrom,
		    const GenomicRegion &chrom_region,
		    const double max_mismatches,
		    FileIterator<MappedRead> &regions, 
		    std::ostream &out, MethStat &meth_stat_collector,
                    const size_t max_length) {
  const string chrom_name(chrom_region.get_chrom());
  size_t i = 0;
  for (i = 0; i < chrom.length() - 1 && regions.first_is_good(); ++i) {
    if (is_cpg(chrom, i)) {
      /* need the "+1" below because of the 'G' in CpG */
      advance(i, i + 1, chrom_region, regions, max_length);
      size_t meth_count = 0, unmeth_count = 0;
      for (vector<MappedRead>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j)
	if (j->r.get_score() <= max_mismatches)
	  add_contribution_cpg(i, *j, meth_count, unmeth_count);
      const size_t total = meth_count + unmeth_count;
      meth_stat_collector.collect(meth_count, total);
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t+\tCpG\t" 
            << meth_count/max(1.0, static_cast<double>(total)) << "\t"
            << total << endl;
      else
        out << chrom_name << "\t" 
            << i << "\t" << i + 1 << "\tCpG:" << total << "\t" 
            << meth_count/max(1.0, static_cast<double>(total)) << "\t+\n";
    }
  }

  for(;  i < chrom.length() - 1; ++i){
    if (is_cpg(chrom, i)) {
      const size_t total = 0;
      meth_stat_collector.collect(total, total);
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t+\tCpG\t" 
            << total << "\t" << total << endl;
      else
        out << chrom_name << "\t" << i << "\t" << i + 1 << "\tCpG:"
            << total << "\t" << total << "\t+\n";
    }
  }
}


static void
scan_chromosome(const string &chrom, const GenomicRegion &chrom_region,
		const double max_mismatches,
		FileIterator<MappedRead> &regions, 
		std::ostream &out, MethStat &meth_stat_collector,
                const size_t max_length) {
  
  const string chrom_name(chrom_region.get_chrom());
  size_t i = 0; 
  size_t chrom_len = chrom.length();
  for (i = 0; i < chrom_len - 1 && regions.first_is_good(); ++i) {
    advance(i, i, chrom_region, regions, max_length);
    if (is_cytosine(chrom[i])) { // && !is_guanine(chrom[i + 1])) {
      size_t meth_count = 0, unmeth_count = 0;
      for (vector<MappedRead>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j)
	if (j->r.get_score() <= max_mismatches)
	  add_contribution_c(i, *j, meth_count, unmeth_count);
      const size_t total = meth_count + unmeth_count;
      meth_stat_collector.collect(meth_count, total);
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t+\t"  
            << cytosine_type_tag(chrom, i, '+') << "\t"
            << meth_count/max(1.0, static_cast<double>(total)) << "\t"
            << total << endl;
      else
        out << chrom_name << "\t" << i << "\t" << i + 1 << "\t" 
            << cytosine_type_tag(chrom, i, '+') << ":"
            << total << "\t" 
            << meth_count/max(1.0, static_cast<double>(total)) << "\t+\n";
    }
    if (is_guanine(chrom[i])) { // && !is_cytosine(chrom[i - 1])) {
      size_t meth_count = 0, unmeth_count = 0;
      for (vector<MappedRead>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j)
	if (j->r.get_score() <= max_mismatches)
	  add_contribution_g(i, *j, meth_count, unmeth_count);
      const size_t total = meth_count + unmeth_count;
      meth_stat_collector.collect(meth_count, total);
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t-\t" 
            << cytosine_type_tag(chrom, i, '-') << "\t" 
            << meth_count/max(1.0, static_cast<double>(total)) << "\t"
            << total << endl;
      else
        out << chrom_name << "\t" << i << "\t" << i + 1 << "\t" 
            << cytosine_type_tag(chrom, i, '-') << ":" 
            << total << "\t" 
            << meth_count/max(1.0, static_cast<double>(total)) << "\t-\n";
    }
  }

  for (; i < chrom_len - 1 ; ++i){
    const size_t total = 0;
    if (is_cytosine(chrom[i])) {
      meth_stat_collector.collect(total, total);      
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t+\t"  
            << cytosine_type_tag(chrom, i, '+') << "\t"
            << total << "\t" << total << endl;
      else
        out << chrom_name << "\t" << i << "\t" << i + 1 << "\t" 
            << cytosine_type_tag(chrom, i, '+') << ":"
            << total << "\t" << total << "\t+\n";
    }
    if (is_guanine(chrom[i])) {
      meth_stat_collector.collect(total, total);       
      if (USE_ALT_OUTPUT)
        out << chrom_name << "\t" << i << "\t+\t"  
            << cytosine_type_tag(chrom, i, '-') << "\t"
            << total << "\t" << total << endl;
      else
        out << chrom_name << "\t" << i << "\t" << i + 1 << "\t" 
            << cytosine_type_tag(chrom, i, '-') << ":"
            << total << "\t" << total << "\t+\n";
    }
  }
}


static void
identify_chromosomes(const bool VERBOSE, const string chrom_file,
		     const string fasta_suffix, vector<string> &chrom_files) {
  if (VERBOSE)
    cerr << "[IDENTIFYING CHROMS] ";
  if (isdir(chrom_file.c_str())) 
    read_dir(chrom_file, fasta_suffix, chrom_files);
  else chrom_files.push_back(chrom_file);
  if (VERBOSE) {
    cerr << "[DONE]" << endl 
	 << "chromosome files found (approx size):" << endl;
    for (vector<string>::const_iterator i = chrom_files.begin();
	 i != chrom_files.end(); ++i)
      cerr << *i << " (" << roundf(get_filesize(*i)/1e06) << "Mbp)" << endl;
    cerr << endl;
  }
}


static void
advance_chromosome(const GenomicRegion &chrom_region, 
		   FileIterator<MappedRead> &regions) {

  while (regions.last_is_good() && 
	 (regions.get_last()->r < chrom_region)) {
    assert(regions.last_is_good());
    regions.increment_last();
  }
  while (regions.first_is_good() && 
	 (regions.get_first()->r < chrom_region))
    regions.increment_first();
}


static void
fix_chrom_names(vector<string> &chrom_names)
{
  // make sure the chrom names don't have spaces

  for (size_t i = 0; i < chrom_names.size(); ++i) {
    const size_t chr_name_end = chrom_names[i].find_first_of(" \t");
    if (chr_name_end != string::npos)
      chrom_names[i].erase(chr_name_end);
  }
}


static 
void
scan_chroms(const bool VERBOSE, const bool PROCESS_NON_CPGS,
	    const QualityChecker &qc,
	    const double max_mismatches,
	    const string &outfile, const vector<string> &chrom_files, 
	    FileIterator<MappedRead> &regions, MethStat &meth_stat_collector,
            const size_t max_length) {
  
  // Get the output stream
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  for (size_t i = 0; i < chrom_files.size(); ++i) {
    const string fn(strip_path(chrom_files[i]));
    if (VERBOSE)
      cerr << "[LOADING CHROM FILE=" << fn << "]" << endl;
    vector<string> chrom_names, chroms;
    read_fasta_file(chrom_files[i].c_str(), chrom_names, chroms);
    fix_chrom_names(chrom_names);
    if (chrom_names.size() > 1) {
      // make sure chrosomes comes in lexical order as how input reads
      // is sorted
      vector<std::pair<string, size_t> > chrom_idx;
      for (size_t j = 0; j < chrom_names.size(); ++j)
        chrom_idx.push_back(std::make_pair(chrom_names[j], j));
      std::sort(chrom_idx.begin(), chrom_idx.end(), Compare());
      
      vector<string> chrom_names_tmp(chrom_names.size()),
        chroms_tmp(chrom_names.size());
      for (size_t j = 0; j < chrom_names.size(); ++j)
      {
        std::swap(chrom_names_tmp[j], chrom_names[chrom_idx[j].second]);
        std::swap(chroms_tmp[j], chroms[chrom_idx[j].second]);
      }
      std::swap(chrom_names_tmp, chrom_names);
      std::swap(chroms_tmp, chroms);
    }
    for (size_t j = 0; j < chroms.size(); ++j) {
      if (VERBOSE) cerr << "[SCANNING=" << chrom_names[j] << "]";
      //TODO: WHAT HAPPENS IF A CHROM IS MISSING??
      const GenomicRegion chrom_region(chrom_names[j], 0, 0);
      advance_chromosome(chrom_region, regions);
      if (PROCESS_NON_CPGS)
	scan_chromosome(qc, chroms[j], chrom_region, max_mismatches,
			regions, out, meth_stat_collector, max_length);
      else scan_chromosome_cpg(qc, chroms[j], chrom_region, max_mismatches,
			       regions, out, meth_stat_collector, max_length);
      if (VERBOSE) cerr << " [DONE]" << endl;
    }
  }
}


static void
scan_chroms(const bool VERBOSE, const bool PROCESS_NON_CPGS,
	    const double max_mismatches,
	    const string &outfile, const vector<string> &chrom_files, 
	    FileIterator<MappedRead> &regions, MethStat &meth_stat_collector,
            const size_t max_length) {

  // Get the output stream
  std::ofstream of;
  if (!outfile.empty()) of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
  
  for (size_t i = 0; i < chrom_files.size(); ++i) {
    const string fn(strip_path(chrom_files[i]));
    if (VERBOSE)
      cerr << "[LOADING CHROM FILE=" << fn << "]" << endl;
    vector<string> chrom_names, chroms;
    read_fasta_file(chrom_files[i].c_str(), chrom_names, chroms);
    fix_chrom_names(chrom_names);
    if (chrom_names.size() > 1) {
      // make sure chrosomes comes in lexical order as how input reads
      // is sorted
      vector<std::pair<string, size_t> > chrom_idx;
      for (size_t j = 0; j < chrom_names.size(); ++j)
        chrom_idx.push_back(std::make_pair(chrom_names[j], j));
      std::sort(chrom_idx.begin(), chrom_idx.end(), Compare());
      
      vector<string> chrom_names_tmp(chrom_names.size()),
        chroms_tmp(chrom_names.size());
      for (size_t j = 0; j < chrom_names.size(); ++j)
      {
        std::swap(chrom_names_tmp[j], chrom_names[chrom_idx[j].second]);
        std::swap(chroms_tmp[j], chroms[chrom_idx[j].second]);
      }
      std::swap(chrom_names_tmp, chrom_names);
      std::swap(chroms_tmp, chroms);
    }
    for (size_t j = 0; j < chroms.size(); ++j) {
      if (VERBOSE) cerr << "[SCANNING=" << chrom_names[j] << "]";
      //TODO: WHAT HAPPENS IF A CHROM IS MISSING??
      const GenomicRegion chrom_region(chrom_names[j], 0, 0);
      advance_chromosome(chrom_region, regions);
      if (PROCESS_NON_CPGS)
	scan_chromosome(chroms[j], chrom_region, max_mismatches,
			regions, out, meth_stat_collector, max_length);
      else scan_chromosome_cpg(chroms[j], chrom_region, max_mismatches,
			       regions, out, meth_stat_collector, max_length);
      if (VERBOSE) cerr << " [DONE]" << endl;
    }
  }
}



int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    bool PROCESS_NON_CPGS = false;
    
    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";
    
    size_t BUFFER_SIZE = 100000;
    size_t max_length = 10000;
    string out_stat;
    double max_mismatches = std::numeric_limits<double>::max();

    double cutoff = -std::numeric_limits<double>::max();
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Program to count "
			   "methylated/unmethylated bases in "
			   "reads mapping over each CpG or C",
			   "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("chrom", 'c', "file or dir of chroms (FASTA format; .fa suffix)",
		      true , chrom_file);
    opt_parse.add_opt("non", 'N', "process non-CpG cytosines", 
		      false , PROCESS_NON_CPGS);
    //!!!!!! OPTION IS HIDDEN BECAUSE USERS DON'T NEED TO CHANGE IT...
    //     opt_parse.add_opt("buffer", 'B', "buffer size (in records, not bytes)", 
    // 		      false , BUFFER_SIZE);
    //     opt_parse.add_opt("suffix", 's', "suffix of FASTA files "
    // 		      "(assumes -c indicates dir)", 
    // 		      false , fasta_suffix);
    //     opt_parse.add_opt("cutoff", 'C', "cutoff for high-quality bases (assumes fastq reads)", 
    // 		      false , cutoff);
    opt_parse.add_opt("max", 'M', "max mismatches (can be fractional)", 
		      false , max_mismatches);
    opt_parse.add_opt("max_length", 'L', "The maximal read length of the input file (default 10000)",
                      false , max_length);
    opt_parse.add_opt("output_stat", 'S', "Name of output file with statistics",
                      false , out_stat);
    opt_parse.add_opt("alt_out", 'A', "use alternative output format",
                      false , USE_ALT_OUTPUT);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "MAX MISMATCHES=" << max_mismatches << endl;

    FASTQScoreType score_format = FASTQ_Solexa;
    
    vector<string> chrom_files;
    identify_chromosomes(VERBOSE, chrom_file, fasta_suffix, chrom_files);
    sort(chrom_files.begin(), chrom_files.end());
    
    MethStat meth_stat_collector;
    
    FileIterator<MappedRead> regions(mapped_reads_file, BUFFER_SIZE);
    if (cutoff != -std::numeric_limits<double>::max()) {
      const QualityChecker qc(score_format, cutoff);
      scan_chroms(VERBOSE, PROCESS_NON_CPGS, qc, max_mismatches, 
		  outfile, chrom_files, regions, meth_stat_collector, max_length);
    }
    else scan_chroms(VERBOSE, PROCESS_NON_CPGS, max_mismatches, 
		     outfile, chrom_files, regions, meth_stat_collector, max_length);
    
    if (VERBOSE || !out_stat.empty()) {
      std::ofstream of;
      if (!out_stat.empty()) of.open(out_stat.c_str());
      std::ostream out(out_stat.empty() ? cerr.rdbuf() : of.rdbuf());
      out << meth_stat_collector << endl;
    }
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
