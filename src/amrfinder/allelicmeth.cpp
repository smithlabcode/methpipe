/* allelicmeth:
 *
 * Copyright (C) 2014-2022 University of Southern California,
 *                         Andrew D. Smith, and Benjamin E. Decato
 *
 * Authors: Andrew D. Smith and Benjamin E. Decato
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <utility>
#include <stdexcept>

#include <cmath>
#include <sstream>
#include <gsl/gsl_sf_gamma.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

#include "Epiread.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::istringstream;
using std::unordered_map;
using std::ostringstream;
using std::max;
using std::min;
using std::runtime_error;

static inline double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


// p(k) =  C(n1, k) C(n2, t - k) / C(n1 + n2, t)
static double
log_hyper_g(const size_t k, const size_t n1, const size_t n2, const size_t t) {
  return (gsl_sf_lnfact(n1) - gsl_sf_lnfact(k) - gsl_sf_lnfact(n1 - k) +
          gsl_sf_lnfact(n2) - gsl_sf_lnfact(t - k) - gsl_sf_lnfact(n2 - (t - k)) -
          (gsl_sf_lnfact(n1 + n2) - gsl_sf_lnfact(t) - gsl_sf_lnfact(n1 + n2 - t)));
}


static double
fishers_exact(size_t a, size_t b, size_t c, size_t d) {
  const size_t m = a + c; // sum of first column
  const size_t n = b + d; // sum of second column
  const size_t k = a + b; // sum of first row
  const double observed = log_hyper_g(a, m, n, k);
  double p = 0.0;
  for (size_t i = (n > k ? 0ul : k - n); i <= std::min(k, m); ++i) {
    const double curr = log_hyper_g(i, m, n, k);
    if (curr <= observed)
      p = log_sum_log(p, curr);
  }
  return exp(p);
}


static size_t
state_pair_to_index(const string &s, const size_t idx) {
  assert(idx < s.length() - 1);
  const char a = s[idx];
  if (a == 'C') {
    const char b = s[idx+1];
    if (b == 'C') return 0;
    if (b == 'T') return 1;
    return 4;
  }
  if (a == 'T') {
    const char b = s[idx+1];
    if (b == 'C') return 2;
    if (b == 'T') return 3;
    return 4;
  }
  return 4;
}


template <class T>
struct PairStateCounter {
  T CC;
  T CT;
  T TC;
  T TT;

  double score() const {
    return (CC*TT > CT*TC) ?
      fishers_exact(CC, CT, TC, TT) : fishers_exact(CT, CC, TT, TC);
  }
  double total() const {return CC + CT + TC + TT;}

  string tostring() const {
    return toa(CC) + '\t' + toa(CT) + '\t' + toa(TC) + '\t' + toa(TT);
  }

  void increment(const size_t state) {
    if (state == 0) ++CC;
    else if (state == 1) ++CT;
    else if (state == 2) ++TC;
    else if (state == 3) ++TT;
  }
};


template <class T> void
fit_states(const epiread &er, vector<PairStateCounter<T> > &counts) {
  for (size_t i = 0; i < er.length() - 1; ++i) {
    const size_t pos = er.pos + i;
    assert(pos < counts.size());
    const size_t curr_state = state_pair_to_index(er.seq, i);
    counts[pos].increment(curr_state);
  }
}


static void
get_chrom_sizes(const bool VERBOSE, const string &epi_file,
                unordered_map<string, size_t> &chrom_sizes) {
  std::ifstream in(epi_file);
  if (!in)
    throw runtime_error("cannot open input file: " + epi_file);

  string chrom;
  epiread er;
  while (in >> er) {
    if (chrom_sizes.find(er.chr) == chrom_sizes.end())
      chrom_sizes[er.chr] = 0;
    chrom_sizes[er.chr] = std::max(chrom_sizes[er.chr], er.pos + er.length());
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////
//////////////  CODE FOR CONVERTING BETWEEN CPG AND BASE PAIR
//////////////  COORDINATES BELOW HERE
//////////////
inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}


static void
collect_cpgs(const string &s, unordered_map<size_t, size_t> &cpgs) {
  const size_t lim = s.length() - 1;
  size_t cpg_count = 0;
  for (size_t i = 0; i < lim; ++i)
    if (is_cpg(s, i)) {
      cpgs[cpg_count] = i;
      ++cpg_count;
    }
}


static void
convert_coordinates(const unordered_map<size_t, size_t> &cpgs,
                    GenomicRegion &region)  {

  auto start_itr(cpgs.find(region.get_start()));
  if (start_itr == end(cpgs))
    throw runtime_error("could not convert:\n" + region.tostring());

  region.set_start(start_itr->second);
  region.set_end(start_itr->second + 1);
}


static void
get_chrom(const bool VERBOSE, const GenomicRegion &r,
          const unordered_map<string, string>& chrom_files,
          GenomicRegion &chrom_region,  string &chrom) {
  const unordered_map<string, string>::const_iterator
    fn(chrom_files.find(r.get_chrom()));
  if (fn == chrom_files.end())
    throw runtime_error("could not find chrom: " + r.get_chrom());
  chrom.clear();
  read_fasta_file(fn->second, r.get_chrom(), chrom);
  if (chrom.empty())
    throw runtime_error("could not find chrom: " + r.get_chrom());
  else {
    chrom_region.set_chrom(r.get_chrom());
  }
}


/* The "sites" in the convert_coordinates function represent pairs of
    consecutive CpG sites */
static void
convert_coordinates(const bool VERBOSE, const string chroms_dir,
                    const string fasta_suffix, vector<GenomicRegion> &sites) {

  unordered_map<string, string> chrom_files;
  identify_and_read_chromosomes(chroms_dir, fasta_suffix, chrom_files);
  if (VERBOSE)
    cerr << "CHROMS:\t" << chrom_files.size() << endl;

  unordered_map<size_t, size_t> cpgs;
  string chrom;
  GenomicRegion chrom_region("chr0", 0, 0);
  for (size_t i = 0; i < sites.size(); ++i) {
    if (!sites[i].same_chrom(chrom_region)) {
      get_chrom(VERBOSE, sites[i], chrom_files, chrom_region, chrom);
      collect_cpgs(chrom, cpgs);
      if (VERBOSE)
        cerr << "CONVERTING: " << chrom_region.get_chrom() << endl;
    }
    convert_coordinates(cpgs, sites[i]);
  }
}


static void
add_cytosine(const string &chrom_name, const size_t start_cpg,
             vector<PairStateCounter<unsigned short> > &counts,
             vector<GenomicRegion> &cytosines) {
  const size_t end_cpg = start_cpg + 1;
  ostringstream s;
  s << counts[start_cpg].score() << "\t"
    << counts[start_cpg].total() << "\t"
    << counts[start_cpg].tostring();
  const string name(s.str());
  cytosines.push_back(GenomicRegion(chrom_name, start_cpg, end_cpg,
                                    name, 0, '+'));
}


static size_t
process_chrom(const bool VERBOSE, const string &chrom_name,
              const vector<epiread> &epireads, vector<GenomicRegion> &cytosines,
              vector<PairStateCounter<unsigned short> > &counts) {
  size_t max_epiread_len = 0;
  for (size_t i = 0; i < epireads.size(); ++i)
    max_epiread_len = std::max(max_epiread_len, epireads[i].length());
  const size_t chrom_cpgs = get_n_cpgs(epireads);
  if (VERBOSE)
    cerr << "PROCESSING: " << chrom_name << " "
         << "[reads: " << epireads.size() << "] "
         << "[cpgs: " << chrom_cpgs << "]" << endl;

  for (size_t i = 0; i < epireads.size(); ++i) {
    fit_states(epireads[i], counts);
  }
  for (size_t i = 0; i < counts.size(); ++i) {
    add_cytosine(chrom_name, i, counts, cytosines);
  }
  return 0;
}


int
main(int argc, const char **argv) {

  try {

    static const string description =
      "computes probability of allele-specific \
       methylation at each tuple of CpGs";

    static const string fasta_suffix = "fa";
    bool VERBOSE = false;

    string outfile;
    string chroms_dir;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description, "<epireads>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("chrom", 'c', "genome sequence file/directory",
                      true, chroms_dir);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string epi_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    unordered_map<string, size_t> chrom_sizes;
    get_chrom_sizes(VERBOSE, epi_file, chrom_sizes);

    if (VERBOSE)
      cerr << "CHROMS: " << chrom_sizes.size() << endl;

    std::ifstream in(epi_file);
    if (!in)
      throw runtime_error("cannot open input file: " + epi_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    string chrom;
    epiread er;
    vector<epiread> epireads;
    vector<GenomicRegion> cytosines;
    vector<PairStateCounter<unsigned short> > counts;
    while (in >> er) {
      if (er.chr != chrom && !chrom.empty()) {
        counts = vector<PairStateCounter<unsigned short> >(chrom_sizes[chrom] - 1);
        process_chrom(VERBOSE, chrom, epireads, cytosines, counts);
        // convert coordinates
        convert_coordinates(VERBOSE, chroms_dir, fasta_suffix, cytosines);
        epireads.clear();
        for( size_t i = 0; i < cytosines.size()-1; ++i ) {
          out << cytosines[i].get_chrom() << "\t"
              << cytosines[i].get_start() << "\t+\tCpG\t"
              << cytosines[i].get_name() << endl;
        }
        cytosines.clear();
      }
      epireads.push_back(er);
      chrom.swap(er.chr);
    }
    if (!chrom.empty()) {
      counts = vector<PairStateCounter<unsigned short> >(chrom_sizes[chrom]);
      process_chrom(VERBOSE, chrom, epireads, cytosines, counts);
      convert_coordinates(VERBOSE, chroms_dir, fasta_suffix, cytosines);
      for (size_t i = 0; i < cytosines.size() - 1; ++i) {
        out << cytosines[i].get_chrom() << "\t"
            << cytosines[i].get_start() << "\t+\tCpG\t"
            << cytosines[i].get_name() << endl;
      }
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
