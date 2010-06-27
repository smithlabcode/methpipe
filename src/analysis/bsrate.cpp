/*    bsrate: a program for determining the rate of bisulfite
 *    conversion in a bisulfite sequencing experiment, with too many
 *    reads to load into memory at once.
 *
 *    Copyright (C) 2009 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
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

// TODO: check that the mapped locations and the read sequences are
// equal in number.


#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "QualityScore.hpp"

#include "bsutils.hpp"
#include "FileIterator.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;


/***********************************************************************
 * FUNCTIONS BELOW ARE FOR FASTQRecord OJBECTS AND USING THE QUALITY
 * INFORAMTION
 */

struct QualityChecker {
  QualityChecker(const FASTQScoreType &score_format, const double c) {
    is_phred = FASTQScoreIsPhred(score_format);
    cutoff = (is_phred) ?
      phred_to_quality_character(error_probability_to_phred(c)) :
      solexa_to_quality_character(error_probability_to_solexa(c));
  }
  bool operator()(const FASTQRecord &s, const size_t position) const {
    return s.second[position] >= cutoff;
  }
  double error_probability(const FASTQRecord &s, const size_t position) const {
    return (is_phred) ?
      phred_to_error_probability(quality_character_to_phred(s.second[position])) :
      solexa_to_error_probability(quality_character_to_solexa(s.second[position]));
  }
  char cutoff;
  bool is_phred;
};

static void
add_contribution_c(const size_t offset, const GenomicRegion &r,
		   const string &s, 
		   vector<size_t> &unconv, vector<size_t> &conv,
		   vector<size_t> &err, vector<double> &qual) {
  if (r.pos_strand()) {
    const size_t position = offset - r.get_start();
    assert(position < s.length());
    if (is_cytosine(s[position])) ++unconv[position];
    else if (is_thymine(s[position])) ++conv[position];
    else ++err[position];
  }
}


static void
add_contribution_g(const size_t offset, const GenomicRegion &r,
		   const string &s, 
		   vector<size_t> &unconv, vector<size_t> &conv,
		   vector<size_t> &err, vector<double> &qual) {
  if (r.neg_strand()) {
    const size_t position = s.length() - (offset - r.get_start() + 1);
    assert(position < s.length());
    if (is_cytosine(s[position])) ++unconv[position];
    else if (is_thymine(s[position])) ++conv[position];
    else ++err[position];
  }
}


static void
add_contribution_c(const QualityChecker &qc,
		   const size_t offset, const GenomicRegion &r,
		   const FASTQRecord &s, 
		   vector<size_t> &unconv, vector<size_t> &conv,
		   vector<size_t> &err, vector<double> &qual) {
  if (r.pos_strand()) {
    const size_t position = offset - r.get_start();
    assert(position < s.first.length());
    if (qc(s, position)) {
      if (is_cytosine(s.first[position])) ++unconv[position];
      else if (is_thymine(s.first[position])) ++conv[position];
      else ++err[position];
      qual[position] += qc.error_probability(s, position);
    }
  }
}


static void
add_contribution_g(const QualityChecker &qc,
		   const size_t offset, const GenomicRegion &r,
		   const FASTQRecord &s, 
		   vector<size_t> &unconv, vector<size_t> &conv,
		   vector<size_t> &err, vector<double> &qual) {
  if (r.neg_strand()) {
    const size_t position = s.first.length() - (offset - r.get_start() + 1);
    assert(position < s.first.length());
    if (qc(s, position)) {
      if (is_cytosine(s.first[position])) ++unconv[position];
      else if (is_thymine(s.first[position])) ++conv[position];
      else ++err[position];
      qual[position] += qc.error_probability(s, position);
    }
  }
}

static void
add_contribution_c(const size_t offset, const GenomicRegion &r,
		   const FASTQRecord &s, 
		   vector<size_t> &unconv, vector<size_t> &conv,
		   vector<size_t> &err, vector<double> &qual) {
  if (r.pos_strand()) {
    const size_t position = offset - r.get_start();
    assert(position < s.first.length());
    if (is_cytosine(s.first[position])) ++unconv[position];
    else if (is_thymine(s.first[position])) ++conv[position];
    else ++err[position];
  }
}


static void
add_contribution_g(const size_t offset, const GenomicRegion &r,
		   const FASTQRecord &s, 
		   vector<size_t> &unconv, vector<size_t> &conv,
		   vector<size_t> &err, vector<double> &qual) {
  if (r.neg_strand()) {
    const size_t position = s.first.length() - (offset - r.get_start() + 1);
    assert(position < s.first.length());
    if (is_cytosine(s.first[position])) ++unconv[position];
    else if (is_thymine(s.first[position])) ++conv[position];
    else ++err[position];
  }
}


static bool
precedes(const GenomicRegion &r, const size_t offset) {
  return r.get_end() <= offset;
}


static bool
succeeds(const GenomicRegion &r, const size_t offset) {
  return r.get_start() > offset;
}


template <class T> void
advance(const size_t first, const size_t last,
	const GenomicRegion &chrom_region, 
	FileIterator<GenomicRegion> &regions,
	FileIterator<T> &reads) {
  while (regions.last_is_good() && 
	 chrom_region.same_chrom(*regions.get_last()) &&
	 !succeeds(*regions.get_last(), last)) {
    regions.increment_last();
    reads.increment_last();
  }
  //   if (regions.last_is_good() != reads.last_is_good())
  //     throw RMAPException("read and map files seem out of sync");
  while (regions.first_is_good() && 
	 chrom_region.same_chrom(*regions.get_first()) &&
	 precedes(*regions.get_first(), first)) {
    regions.increment_first();
    reads.increment_first();
  }
  //   if (regions.first_is_good() != reads.first_is_good())
  //     throw RMAPException("read and map files seem out of sync");
}


template <class T> void
scan_chromosome(const string &chrom, const GenomicRegion &chrom_region,
		const double max_mismatches,
		FileIterator<GenomicRegion> &regions, 
		FileIterator<T> &reads,
		vector<size_t> &unconv_count_pos, vector<size_t> &conv_count_pos,
		vector<size_t> &unconv_count_neg, vector<size_t> &conv_count_neg,
		vector<size_t> &err_pos, vector<size_t> &err_neg,
		vector<double> &qual_pos, vector<double> &qual_neg) {
  const string chrom_name(chrom_region.get_chrom());
  
  for (size_t i = 0; i < chrom.length() - 1 && regions.first_is_good(); ++i) {
    advance(i, i, chrom_region, regions, reads);
    if (is_cytosine(chrom[i]) && !is_guanine(chrom[i + 1])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k) {
	if (j->get_score() <= max_mismatches)
	  add_contribution_c(i, *j, *k, 
			     unconv_count_pos, conv_count_pos, 
			     err_pos, qual_pos);
      }
    }
    if (is_guanine(chrom[i]) && !is_cytosine(chrom[i - 1])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k)
	if (j->get_score() <= max_mismatches)
	  add_contribution_g(i, *j, *k, 
			     unconv_count_neg, conv_count_neg, 
			     err_neg, qual_neg);
    }
  }
}


template <class T> void
scan_chromosome_cpg(const string &chrom, const GenomicRegion &chrom_region,
		    const double max_mismatches,
		    FileIterator<GenomicRegion> &regions, 
		    FileIterator<T> &reads,
		    vector<size_t> &unconv_count_pos, vector<size_t> &conv_count_pos,
		    vector<size_t> &unconv_count_neg, vector<size_t> &conv_count_neg,
		    vector<size_t> &err_pos, vector<size_t> &err_neg,
		    vector<double> &qual_pos, vector<double> &qual_neg) {
  const string chrom_name(chrom_region.get_chrom());
  
  for (size_t i = 0; i < chrom.length() - 1 && regions.first_is_good(); ++i) {
    advance(i, i, chrom_region, regions, reads);
    if (is_cytosine(chrom[i])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k) {
	if (j->get_score() <= max_mismatches)
	  add_contribution_c(i, *j, *k, 
			     unconv_count_pos, conv_count_pos, 
			     err_pos, qual_pos);
      }
    }
    if (is_guanine(chrom[i])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k)
	if (j->get_score() <= max_mismatches)
	  add_contribution_g(i, *j, *k, 
			     unconv_count_neg, conv_count_neg, 
			     err_neg, qual_neg);
    }
  }
}

template <class T> void
scan_chromosome(const string &chrom, const GenomicRegion &chrom_region,
		const QualityChecker &qc,
		const double max_mismatches,
		FileIterator<GenomicRegion> &regions, 
		FileIterator<T> &reads,
		vector<size_t> &unconv_count_pos, vector<size_t> &conv_count_pos,
		vector<size_t> &unconv_count_neg, vector<size_t> &conv_count_neg,
		vector<size_t> &err_pos, vector<size_t> &err_neg,
		vector<double> &qual_pos, vector<double> &qual_neg) {
  const string chrom_name(chrom_region.get_chrom());
  
  for (size_t i = 0; i < chrom.length() - 1 && regions.first_is_good(); ++i) {
    advance(i, i, chrom_region, regions, reads);
    if (is_cytosine(chrom[i]) && !is_guanine(chrom[i + 1])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k) {
	if (j->get_score() <= max_mismatches)
	  add_contribution_c(qc, i, *j, *k, 
			     unconv_count_pos, conv_count_pos, 
			     err_pos, qual_pos);
      }
    }
    if (is_guanine(chrom[i]) && !is_cytosine(chrom[i - 1])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k)
	if (j->get_score() <= max_mismatches)
	  add_contribution_g(qc, i, *j, *k, 
			     unconv_count_neg, conv_count_neg, 
			     err_neg, qual_neg);
    }
  }
}


template <class T> void
scan_chromosome_cpg(const string &chrom, const GenomicRegion &chrom_region,
		    const QualityChecker &qc,
		    const double max_mismatches,
		    FileIterator<GenomicRegion> &regions, 
		    FileIterator<T> &reads,
		    vector<size_t> &unconv_count_pos, vector<size_t> &conv_count_pos,
		    vector<size_t> &unconv_count_neg, vector<size_t> &conv_count_neg,
		    vector<size_t> &err_pos, vector<size_t> &err_neg,
		    vector<double> &qual_pos, vector<double> &qual_neg) {
  const string chrom_name(chrom_region.get_chrom());
  
  for (size_t i = 0; i < chrom.length() - 1 && regions.first_is_good(); ++i) {
    advance(i, i, chrom_region, regions, reads);
    if (is_cytosine(chrom[i])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k) {
	if (j->get_score() <= max_mismatches)
	  add_contribution_c(qc, i, *j, *k, 
			     unconv_count_pos, conv_count_pos, 
			     err_pos, qual_pos);
      }
    }
    if (is_guanine(chrom[i])) {
      typename vector<T>::const_iterator k(reads.get_first());
      for (vector<GenomicRegion>::const_iterator j(regions.get_first());
	   j != regions.get_last(); ++j, ++k)
	if (j->get_score() <= max_mismatches)
	  add_contribution_g(qc, i, *j, *k, 
			     unconv_count_neg, conv_count_neg, 
			     err_neg, qual_neg);
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


template <class T> void
advance_chromosome(const GenomicRegion &chrom_region, 
		   FileIterator<GenomicRegion> &regions, 
		   FileIterator<T> &reads) {
  while (regions.last_is_good() && 
	 (*regions.get_last() < chrom_region)) {
    assert(regions.last_is_good());
    regions.increment_last();
    reads.increment_last();
  }
  while (regions.first_is_good() && 
	 (*regions.get_first() < chrom_region)) {
    regions.increment_first();
    reads.increment_first();
  }
}


template <class T> 
void
scan_chroms(const bool VERBOSE, const bool PROCESS_CPGS,
	    const double max_mismatches,
	    const string &outfile, const vector<string> &chrom_files, 
	    FileIterator<GenomicRegion> &regions,
	    FileIterator<T> &reads) {

  const size_t read_len = regions.get_first()->get_width();
  
  vector<size_t> unconv_count_pos(read_len, 0ul);
  vector<size_t> conv_count_pos(read_len, 0ul);
  vector<size_t> unconv_count_neg(read_len, 0ul);
  vector<size_t> conv_count_neg(read_len, 0ul);
  vector<size_t> err_pos(read_len, 0ul), err_neg(read_len, 0ul);
  vector<double> qual_pos(read_len, 0.0), qual_neg(read_len, 0.0);
  
  for (size_t i = 0; i < chrom_files.size(); ++i) {
    const string fn(strip_path_and_suffix(chrom_files[i]));
    if (VERBOSE)
      cerr << "[LOADING CHROM FILE=" << fn << "]";
    vector<string> chrom_names, chroms;
    read_fasta_file(chrom_files[i].c_str(), chrom_names, chroms);
    for (size_t j = 0; j < chroms.size(); ++j) {
      if (VERBOSE) cerr << "[SCANNING=" << chrom_names[j] << "]";
      //TODO: WHAT HAPPENS IF A CHROM IS MISSING??
      const GenomicRegion chrom_region(chrom_names[j], 0, 0);
      advance_chromosome(chrom_region, regions, reads);
      if (PROCESS_CPGS)
	scan_chromosome_cpg(chroms[j], chrom_region, max_mismatches,
			    regions, reads, 
			    unconv_count_pos, conv_count_pos,
			    unconv_count_neg, conv_count_neg,
			    err_pos, err_neg, qual_pos, qual_neg);
      else scan_chromosome(chroms[j], chrom_region, max_mismatches, 
			   regions, reads, 
			   unconv_count_pos, conv_count_pos,
			   unconv_count_neg, conv_count_neg,
			   err_pos, err_neg, qual_pos, qual_neg);
    }
    if (VERBOSE) cerr << " [DONE]" << endl;
  }
  
  std::ostream *out = (outfile.empty()) ? &cout : new std::ofstream(outfile.c_str());
  *out << "BASE" << '\t'
       << "PTOT" << '\t'
       << "PCONV" << '\t'
       << "PRATE" << '\t'
       << "NTOT" << '\t'
       << "NCONV" << '\t'
       << "NRATE" << '\t'
       << "BTHTOT" << '\t'
       << "BTHCONV" << '\t'
       << "BTHRATE" << '\t'
       << "ERR" << '\t'
       << "ALL" << '\t'
       << "ERRRATE" << '\t'
       << "QUAL" << endl;
  static const size_t precision_val = 5;
  for (size_t i = 0; i < read_len; ++i) {
    const double total_pos = unconv_count_pos[i] + conv_count_pos[i];
    const double total_neg = unconv_count_neg[i] + conv_count_neg[i];
    const double total = total_pos + total_neg;
    *out << i << "\t";

    out->precision(precision_val);
    out->width(precision_val + 2);
    *out << static_cast<size_t>(total_pos) << "\t" 
	 << conv_count_pos[i] << "\t" 
	 << conv_count_pos[i]/total_pos << "\t";

    out->precision(precision_val);
    out->width(precision_val + 2);
    *out << static_cast<size_t>(total_neg) << "\t" 
	 << conv_count_neg[i] << "\t" 
	 << conv_count_neg[i]/total_neg << "\t";
    
    out->precision(precision_val);
    out->width(precision_val + 2);
    *out << static_cast<size_t>(total) << "\t" 
	 << conv_count_pos[i] + conv_count_neg[i] << "\t" 
	 << (conv_count_pos[i] + conv_count_neg[i])/total << '\t';
    
    out->precision(precision_val);
    out->width(precision_val + 2);
    *out << err_pos[i] + err_neg[i] << "\t" 
	 << static_cast<size_t>(total + err_pos[i] + err_neg[i]) << "\t" 
	 << (err_pos[i] + err_neg[i])/(err_pos[i] + err_neg[i] + total) << '\t';
    
    out->precision(precision_val);
    out->width(precision_val + 2);
    *out << (qual_pos[i] + qual_neg[i])/total << endl;
  }
  if (out != &cout) delete out;
}

template <class T> 
void
scan_chroms(const bool VERBOSE, const bool PROCESS_CPGS,
	    const QualityChecker &qc,
	    const double max_mismatches,
	    const string &outfile, const vector<string> &chrom_files, 
	    FileIterator<GenomicRegion> &regions,
	    FileIterator<T> &reads) {

  const size_t read_len = regions.get_first()->get_width();

  vector<size_t> unconv_count_pos(read_len, 0ul);
  vector<size_t> conv_count_pos(read_len, 0ul);
  vector<size_t> unconv_count_neg(read_len, 0ul);
  vector<size_t> conv_count_neg(read_len, 0ul);
  vector<size_t> err_pos(read_len, 0ul), err_neg(read_len, 0ul);
  vector<double> qual_pos(read_len, 0.0), qual_neg(read_len, 0.0);
  
  for (size_t i = 0; i < chrom_files.size(); ++i) {
    const string fn(strip_path_and_suffix(chrom_files[i]));
    if (VERBOSE)
      cerr << "[LOADING CHROM FILE=" << fn << "]";
    vector<string> chrom_names, chroms;
    read_fasta_file(chrom_files[i].c_str(), chrom_names, chroms);
    for (size_t j = 0; j < chroms.size(); ++j) {
      if (VERBOSE) cerr << "[SCANNING=" << chrom_names[j] << "]";
      //TODO: WHAT HAPPENS IF A CHROM IS MISSING??
      const GenomicRegion chrom_region(chrom_names[j], 0, 0);
      advance_chromosome(chrom_region, regions, reads);
      if (PROCESS_CPGS)
	scan_chromosome_cpg(chroms[j], chrom_region, qc, max_mismatches,
			    regions, reads, 
			    unconv_count_pos, conv_count_pos,
			    unconv_count_neg, conv_count_neg,
			    err_pos, err_neg, qual_pos, qual_neg);
      else scan_chromosome(chroms[j], chrom_region, qc, max_mismatches, 
			   regions, reads, 
			   unconv_count_pos, conv_count_pos,
			   unconv_count_neg, conv_count_neg,
			   err_pos, err_neg, qual_pos, qual_neg);
    }
    if (VERBOSE) cerr << " [DONE]" << endl;
  }
  
  std::ostream *out = (outfile.empty()) ? &cout : new std::ofstream(outfile.c_str());
  *out << "BASE" << '\t'
       << "PTOT" << '\t'
       << "PCONV" << '\t'
       << "PRATE" << '\t'
       << "NTOT" << '\t'
       << "NCoNV" << '\t'
       << "NRATE" << '\t'
       << "BTHTOT" << '\t'
       << "BTHCONV" << '\t'
       << "BTHRATE" << '\t'
       << "ERR" << '\t'
       << "ALL" << '\t'
       << "ERRRATE" << endl;

  static const size_t precision_val = 5;
  for (size_t i = 0; i < read_len; ++i) {
    const double total_pos = unconv_count_pos[i] + conv_count_pos[i];
    const double total_neg = unconv_count_neg[i] + conv_count_neg[i];
    const double total = total_pos + total_neg;
    *out << i << "\t";

    out->precision(precision_val);
    out->width(precision_val + 2);
    *out << static_cast<size_t>(total_pos) << "\t" 
	 << conv_count_pos[i] << "\t" << conv_count_pos[i]/total_pos << "\t";
    
    out->precision(precision_val);
    out->width(precision_val + 2);
    *out << static_cast<size_t>(total_neg) << "\t" 
	 << conv_count_neg[i] << "\t" << conv_count_neg[i]/total_neg << "\t";
    
    out->precision(precision_val);
    out->width(precision_val + 2);
    *out << static_cast<size_t>(total) << "\t" 
	 << conv_count_pos[i] + conv_count_neg[i] << "\t" 
	 << (conv_count_pos[i] + conv_count_neg[i])/total << '\t';
    
    out->precision(precision_val);
    out->width(precision_val + 2);
    *out << static_cast<size_t>(err_pos[i] + err_neg[i]) << "\t" 
	 << static_cast<size_t>(total_pos + total_neg + 
				err_pos[i] + err_neg[i]) << "\t" 
	 << (err_pos[i] + err_neg[i])/(err_pos[i] + err_neg[i] + total) << '\t';
    
    out->precision(precision_val);
    out->width(precision_val + 2);
    *out << qual_pos[i]/total_pos << "\t" 
	 << qual_neg[i]/total_neg << "\t" 
	 << (qual_pos[i] + qual_neg[i])/total << endl;
  }
  if (out != &cout) delete out;
}


int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    bool PROCESS_CPGS = false;
    
    string mapped_file;
    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";

    double cutoff = -std::numeric_limits<double>::max();
    
    size_t BUFFER_SIZE = 100000;
    double max_mismatches = std::numeric_limits<double>::max();
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "a program for determining the "
			   "rate of bisulfite conversion in a "
			   "bisulfite sequencing experiment, with "
			   "too many reads to load into memory at once.",
			   "<FASTA_READS>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("chrom", 'c', "FASTA file or dir containing chromosome(s)", 
		      true , chrom_file);
    opt_parse.add_opt("suffix", 's', "suffix of FASTA files "
		      "(assumes -c indicates dir)", 
		      false , fasta_suffix);
    opt_parse.add_opt("mapped", 'm', "file of mapped locations", 
		      true , mapped_file);
    opt_parse.add_opt("all", 'A', "process all Cs", 
		      false , PROCESS_CPGS);
    opt_parse.add_opt("buffer", 'B', "buffer size (in records, not bytes)", 
		      false , BUFFER_SIZE);
    opt_parse.add_opt("cutoff", 'C', "cutoff for high-quality bases (assumes fastq reads)", 
		      false , cutoff);
    opt_parse.add_opt("max", 'M', "max mismatches (can be fractional)", 
		      false , max_mismatches);
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
    const string reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "MAX MISMATCHES=" << max_mismatches << endl;

    const bool FASTQ = is_fastq(reads_file);
    if (VERBOSE)
      cerr << "READS FILE FORMAT: " << ((FASTQ) ? "FASTQ" : "FASTA") << endl;

    FASTQScoreType score_format = (FASTQ) ?
      fastq_score_type(reads_file) : FASTQ_Solexa;

    if (VERBOSE && FASTQ)
      cerr << "SCORE FORMAT: " 
	   << ((FASTQScoreIsPhred(score_format)) ? "Phred" : "Solexa") << endl;
    
    vector<string> chrom_files;
    identify_chromosomes(VERBOSE, chrom_file, fasta_suffix, chrom_files);
    sort(chrom_files.begin(), chrom_files.end());
    
    FileIterator<GenomicRegion> regions(mapped_file, BUFFER_SIZE);
    if (FASTQ) {
      FileIterator<FASTQRecord> reads(reads_file, BUFFER_SIZE);
      if (cutoff != -std::numeric_limits<double>::max()) {
	const QualityChecker qc(score_format, cutoff);
	scan_chroms(VERBOSE, PROCESS_CPGS, qc, max_mismatches,
		    outfile, chrom_files, regions, reads);
      }
      else scan_chroms(VERBOSE, PROCESS_CPGS, max_mismatches,
		       outfile, chrom_files, regions, reads);
    }
    else {
      FileIterator<string> reads(reads_file, BUFFER_SIZE);
      scan_chroms(VERBOSE, PROCESS_CPGS, max_mismatches, 
		  outfile, chrom_files, regions, reads);
    }
  }
  catch (const RMAPException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
