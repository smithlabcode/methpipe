/* Copyright (C) 2013 University of Southern California
 *                    Andrew D Smith and Jenny Qu
 *
 * Author: Jenny Qu and Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include <string>
#include <vector>
#include <iostream>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"


using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::unordered_map;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}



static void
build_coordinate_converter(const unordered_map<string, string> &chrom_files,
                           const string &chrom,
                           unordered_map<size_t, size_t> &cpg_lookup) {

  const unordered_map<string, string>::const_iterator
    file_name(chrom_files.find(chrom));
  if (file_name == chrom_files.end())
    throw SMITHLABException("chrom file not found for chrom: " + chrom);

  vector<string> dummy, chrom_seq;
  read_fasta_file(file_name->second.c_str(), dummy, chrom_seq);
  if (chrom_seq.size() > 1)
    throw SMITHLABException("multiple chroms/file: " + file_name->second);

  const size_t lim = chrom_seq.front().length() - 1;
  size_t cpg_count = 0;
  for (size_t i = 0; i < lim; ++i)
    if (is_cpg(chrom_seq.front(), i)) {
      cpg_lookup[cpg_count] = i;
      ++cpg_count;
    }
}



static size_t
convert_coordinates(const unordered_map<size_t, size_t> &cpgs,
                    const size_t position)  {
  const unordered_map<size_t, size_t>::const_iterator
    i(cpgs.find(position));
  if (i == cpgs.end())
    throw SMITHLABException("could not convert:\n" + toa(position));
  return i->second;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////
/////// CODE FOR WORKING WITH EPIREADS BELOW HERE
///////



struct epiread {
  bool operator<(const epiread &other) const {
    return (chr < other.chr || (chr == other.chr && pos < other.pos));
  }
  size_t get_end() const {return pos + seq.length();}
  size_t length() const {return seq.length();}
  bool flip_states();

  string chr;
  size_t pos;
  string seq;
};



bool
epiread::flip_states() {
  const size_t meth_states_count =
    std::count(seq.begin(), seq.end(), 'C');
  if (meth_states_count < 0.5*length()) {
    for (size_t i = 0; i < seq.length(); ++i)
      seq[i] = (seq[i] == 'T') ? 'C' : ((seq[i] == 'C') ? 'T' : seq[i]);
    return true;
  }
  return false;
}



static std::istream&
operator>>(std::istream &in, epiread &er) {
  string buffer;
  if (getline(in, buffer)) {
    std::istringstream is(buffer);
    if (!(is >> er.chr >> er.pos >> er.seq))
      throw SMITHLABException("malformed epiread line:\n" + buffer);
  }
  return in;
}



std::ostream&
operator<<(std::ostream &out, const epiread &er) {
  return out << er.chr << '\t' << er.pos << '\t' << er.seq;
}

static bool
check_sorted(const vector<epiread> &epireads){
  for (size_t i = 1; i < epireads.size(); ++i)
    if (epireads[i] < epireads[i - 1])
      return false;
  return true;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////
/////// CODE FOR ACTUAL ENTROPY CALCULATIONS BELOW HERE
///////



static double
compute_prob_read_has_state(const vector<double> &site_probs,
                            const size_t start_cpg, const size_t end_cpg,
                            const size_t state, const epiread &er) {
  double prob = 1.0;
  size_t cpg_pos = start_cpg;
  const size_t lim = end_cpg - start_cpg;
  for (size_t i = 0; i < lim; ++i) {
    const char curr_state = ((state & (1ul << i)) > 0) ? 'C' : 'T';
    const size_t er_idx = (cpg_pos >= er.pos) ? cpg_pos - er.pos :
      std::numeric_limits<size_t>::max();
    if (er_idx < er.length() && er.seq[er_idx] != 'N')
      prob *= static_cast<double>(curr_state == er.seq[er_idx]);
    else {
      if (curr_state == 'C')
        prob *= site_probs[cpg_pos];
      else prob *= (1.0 - site_probs[cpg_pos]);
    }
    ++cpg_pos;
  }
  return prob;
}



static double
compute_entropy_for_window(const vector<double> &site_probs,
                           const vector<epiread> &epireads,
                           const size_t start_idx,
                           const size_t end_idx,
                           const size_t start_cpg,
                           const size_t end_cpg,
                           size_t &reads_in_window) {

  const size_t n_states = 1ul << (end_cpg - start_cpg);

  double entropy = 0.0;
  for (size_t i = 0; i < n_states; ++i) {

    double state_prob = 0.0;
    reads_in_window = 0;
    for (size_t j = start_idx; j < end_idx; ++j)
      if (epireads[j].get_end() > start_cpg && epireads[j].pos < end_cpg) {
        state_prob += compute_prob_read_has_state(site_probs, start_cpg,
                                                  end_cpg, i, epireads[j]);
        ++reads_in_window;
      }
    if (reads_in_window > 0) {
      state_prob /= reads_in_window;
      entropy += (state_prob > 0.0) ? state_prob*log2(state_prob) : 0.0;
    }
  }
  return entropy;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////
/////// CODE FOR SLIDING THE WINDOW ALONG THE CHROMOSOME BELOW HERE
///////



/* This function just basically computes the same thing as methcounts
   output, so that unobserved states can be imputed  */
static void
compute_site_probs(const size_t n_cpgs, const vector<epiread> &epireads,
                   vector<double> &site_probs) {

  site_probs = vector<double>(n_cpgs);
  vector<size_t> totals(n_cpgs);

  for (size_t i = 0; i < epireads.size(); ++i) {
    const size_t len = epireads[i].length();
    size_t idx = epireads[i].pos;
    for (size_t j = 0; j < len; ++j, ++idx) {
      site_probs[idx] += (epireads[i].seq[j] == 'C');
      totals[idx] += (epireads[i].seq[j] != 'N');
    }
  }
  for (size_t i = 0; i < site_probs.size(); ++i)
    if (totals[i] > 0.0)
      site_probs[i] /= totals[i];
}



static void
move_start_index(const size_t max_epiread_len,
                 const vector<epiread> &epireads,
                 const size_t start_cpg, size_t &idx) {
  while (idx < epireads.size() &&
         epireads[idx].get_end() + max_epiread_len <= start_cpg)
    ++idx;
}



static void
move_end_index(const vector<epiread> &epireads,
               const size_t start_cpg, const size_t cpg_window, size_t &idx) {
  while (idx < epireads.size() && epireads[idx].pos < start_cpg + cpg_window)
    ++idx;
}



static void
process_chrom(const bool VERBOSE, const size_t cpg_window,
              const vector<epiread> &epireads,
              const unordered_map<size_t, size_t> &cpg_lookup,
              std::ostream &out) {

  const string chrom(epireads.front().chr);
  if (!check_sorted(epireads))
    throw SMITHLABException("epireads not sorted in chrom: " + chrom);

  const size_t n_cpgs = cpg_lookup.size();
  if (VERBOSE)
    cerr << "processing chrom: " << chrom
         << " (cpgs = " << n_cpgs << ")" << endl;

  vector<double> site_probs;
  compute_site_probs(n_cpgs, epireads, site_probs);

  size_t max_epiread_len = 0;
  for (size_t i = 0; i < epireads.size(); ++i)
    max_epiread_len = std::max(max_epiread_len, epireads[i].length());

  size_t start_cpg = 0;
  size_t start_idx = 0, end_idx = 0;
  while (start_cpg + cpg_window < n_cpgs) {

    move_start_index(max_epiread_len, epireads, start_cpg, start_idx);
    move_end_index(epireads, start_cpg, cpg_window, end_idx);

    size_t reads_used = 0;
    const double entropy =
      compute_entropy_for_window(site_probs, epireads, start_idx,
                                 end_idx, start_cpg, start_cpg + cpg_window,
                                 reads_used);

    out << chrom << '\t'
        << convert_coordinates(cpg_lookup, start_cpg + cpg_window/2) << '\t'
        << "+\tCpG\t"
        << entropy << '\t'
        << reads_used << endl;

    ++start_cpg;
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

int
main(int argc, const char **argv) {

  try {

    static const string fasta_suffix = "fa";

    bool VERBOSE = false;
    bool FLIP_MAJORITY_STATE = false;

    size_t cpg_window = 4;
    string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "compute methylation entropy in sliding window",
                           "<chroms> <epireads-file>");
    opt_parse.add_opt("window", 'w', "number of CpGs in sliding window",
                      false, cpg_window);
    opt_parse.add_opt("flip", 'F', "flip read majority state to meth",
                      false, FLIP_MAJORITY_STATE);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false, outfile);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string chroms_dir = leftover_args.front();
    const string epi_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ifstream in(epi_file.c_str());
    if (!in)
      throw SMITHLABException("cannot open input file: " + epi_file);

    unordered_map<string, string> chrom_files;
    identify_chromosomes(chroms_dir, fasta_suffix, chrom_files);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    vector<epiread> epireads;
    epiread tmp_er;
    while (in >> tmp_er) {
      if (!epireads.empty() && tmp_er.chr != epireads.back().chr) {
        unordered_map<size_t, size_t> cpg_lookup;
        build_coordinate_converter(chrom_files,
                                   epireads.back().chr, cpg_lookup);
        process_chrom(VERBOSE, cpg_window, epireads, cpg_lookup, out);
        epireads.clear();
      }
      if (FLIP_MAJORITY_STATE)
        tmp_er.flip_states();
      epireads.push_back(tmp_er);
    }
    if (!epireads.empty()) {
      unordered_map<size_t, size_t> cpg_lookup;
      build_coordinate_converter(chrom_files, epireads.back().chr, cpg_lookup);
      process_chrom(VERBOSE, cpg_window, epireads, cpg_lookup, out);
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
