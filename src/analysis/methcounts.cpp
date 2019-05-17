/*    methcounts: a program for counting the methylated and
 *    unmethylated reads mapping over each CpG or C
 *
 *    Copyright (C) 2011-2014 University of Southern California and
 *                            Andrew D. Smith
 *
 *    Authors: Andrew D. Smith and Song Qiang
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
#include <fstream>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"
#include "MethpipeFiles.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unordered_map;
using std::runtime_error;



/* The three functions below here should probably be moved into
   bsutils.hpp. I am not sure if the DDG function is needed, but it
   seems like if one considers strand, and the CHH is not symmetric,
   then one needs this. Also, Qiang should be consulted on this
   because he spent much time thinking about it in the context of
   plants. */
bool
is_chh(const std::string &s, size_t i) {
  return (i < (s.length() - 2)) &&
    is_cytosine(s[i]) &&
    !is_guanine(s[i + 1]) &&
    !is_guanine(s[i + 2]);
}


bool
is_ddg(const std::string &s, size_t i) {
  return (i < (s.length() - 2)) &&
    !is_cytosine(s[i]) &&
    !is_cytosine(s[i + 1]) &&
    is_guanine(s[i + 2]);
}


bool
is_c_at_g(const std::string &s, size_t i) {
  return (i < (s.length() - 2)) &&
    is_cytosine(s[i]) &&
    !is_cytosine(s[i + 1]) &&
    !is_guanine(s[i + 1]) &&
    is_guanine(s[i + 2]);
}


/* Right now the CountSet objects below are much larger than they need
   to be, for the things we are computing. However, it's not clear
   that the minimum information would really put the memory
   requirement of the program into a more reasonable range, so keeping
   all the information seems reasonable. */

template <class count_type>
struct CountSet {

  string tostring() const {
    std::ostringstream oss;
    oss << pA << '\t' << pC << '\t' << pG << '\t' << pT << '\t'
        << nA << '\t' << nC << '\t' << nG << '\t' << nT << '\t' << N;
    return oss.str();
  }
  void add_count_pos(const char x) {
    if (x == 'T') ++pT; // conditions ordered for efficiency
    else if (x == 'C') ++pC;
    else if (x == 'G') ++pG;
    else if (x == 'A') ++pA;
    else ++N;
  }
  void add_count_neg(const char x) {
    if (x == 'T') ++nT; // conditions ordered for efficiency
    else if (x == 'C') ++nC;
    else if (x == 'G') ++nG;
    else if (x == 'A') ++nA;
    else ++N;
  }
  count_type pos_total() const {return pA + pC + pG + pT;}
  count_type neg_total() const {return nA + nC + nG + nT;}

  count_type unconverted_cytosine() const {return pC;}
  count_type converted_cytosine() const {return pT;}
  count_type unconverted_guanine() const {return nC;}
  count_type converted_guanine() const {return nT;}

  // using "int" here because it is smaller
  count_type pA, pC, pG, pT;
  count_type nA, nC, nG, nT;
  count_type N;
};


/* The "tag" returned by this function should be exclusive, so that
 * the order of checking conditions doesn't matter. There is also a
 * bit of a hack in that the unsigned "pos" could wrap, but this still
 * works as long as the chromosome size is not the maximum size of a
 * size_t.
 */
string
get_methylation_context_tag_from_genome(const string &s, const size_t pos) {
  if (is_cytosine(s[pos])) {
    if (is_cpg(s, pos)) return "CpG";
    else if (is_chh(s, pos)) return "CHH";
    else if (is_c_at_g(s, pos)) return "CXG";
    else return "CCG";
  }
  if (is_guanine(s[pos])) {
    if (is_cpg(s, pos - 1)) return "CpG";
    else if (is_ddg(s, pos - 2)) return "CHH";
    else if (is_c_at_g(s, pos - 2)) return "CXG";
    else return "CCG";
  }
  return "N";
}


template <class count_type>
static void
count_states_pos(const size_t chrom_len, const MappedRead &r,
                 vector<CountSet<count_type> > &counts) {
  const size_t width = r.r.get_width();

  size_t position = r.r.get_start();
  assert(position < chrom_len);
  for (size_t i = 0; i < width; ++i, ++position)
    if (position < chrom_len)
      counts[position].add_count_pos(r.seq[i]);
}


template <class count_type>
static void
count_states_neg(const size_t chrom_len, const MappedRead &r,
                 vector<CountSet<count_type> > &counts) {
  const size_t width = r.r.get_width();

  size_t position = r.r.get_start() + width - 1;
  assert(r.r.get_start() < chrom_len);
  for (size_t i = 0; i < width; ++i, --position)
    if (position < chrom_len)
      counts[position].add_count_neg(r.seq[i]);
}


/* This "has_mutated" function looks on the opposite strand to see
 * if the apparent conversion from C->T was actually already in the
 * DNA because of a mutation or SNP.
 */
template <class count_type>
static bool
has_mutated(const char base, const CountSet<count_type> &cs) {
  static const double MUTATION_DEFINING_FRACTION = 0.5;
  return is_cytosine(base) ?
    (cs.nG < MUTATION_DEFINING_FRACTION*(cs.neg_total())) :
    (cs.pG < MUTATION_DEFINING_FRACTION*(cs.pos_total()));
}

inline static bool
is_cpg_site(const string &s, const size_t pos) {
  return (is_cytosine(s[pos]) ? is_guanine(s[pos+1]) :
          (is_guanine(s[pos]) ?
           (pos > 0 && is_cytosine(s[pos - 1])) : false));
}

template <class count_type>
static void
write_output(std::ostream &out,
             const string &chrom_name, const string &chrom,
             const vector<CountSet<count_type> > &counts,
             bool CPG_ONLY) {

  for (size_t i = 0; i < counts.size(); ++i) {
    const char base = chrom[i];
    if (is_cytosine(base) || is_guanine(base)) {
      const double unconverted = is_cytosine(base) ?
        counts[i].unconverted_cytosine() : counts[i].unconverted_guanine();
      const double converted = is_cytosine(base) ?
        counts[i].converted_cytosine() : counts[i].converted_guanine();
      const double meth = unconverted/(converted + unconverted);
      const string tag = get_methylation_context_tag_from_genome(chrom, i) +
        (has_mutated(base, counts[i]) ? "x" : "");
      if (!CPG_ONLY || is_cpg_site(chrom, i)) {
        methpipe::write_site(out, chrom_name, i,
                             (is_cytosine(base) ? "+" : "-"),
                             tag, meth, converted + unconverted);
      }
    }
  }
}


inline size_t
get_chrom_id(const string &chrom_name,
             const unordered_map<string, size_t> &cl) {

  const unordered_map<string, size_t>::const_iterator
    the_chrom(cl.find(chrom_name));
  if (the_chrom == end(cl))
    throw runtime_error("could not find chrom: " + chrom_name);

  return the_chrom->second;
}

int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool CPG_ONLY = false;

    string chrom_file;
    string outfile;
    string fasta_suffix = "fa";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "get methylation levels from "
                           "mapped WGBS reads", "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("chrom", 'c', "file or dir of chroms (FASTA format; "
                      ".fa suffix)", true , chrom_file);
    opt_parse.add_opt("suffix", 's', "suffix of FASTA files "
                      "(assumes -c specifies dir)", false , fasta_suffix);
    opt_parse.add_opt("cpg-only", 'n', "print only CpG context cytosines",
                      false, CPG_ONLY);
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
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (!outfile.empty() && !is_valid_output_file(outfile))
      throw runtime_error("bad output file: " + outfile);

    vector<string> chrom_files;
    if (isdir(chrom_file.c_str()))
      read_dir(chrom_file, fasta_suffix, chrom_files);
    else chrom_files.push_back(chrom_file);

    if (VERBOSE)
      cerr << "n_chrom_files: " << chrom_files.size() << endl;

    vector<string> chroms, names;
    unordered_map<string, size_t> chrom_lookup;
    vector<size_t> chrom_sizes;
    for (auto i(begin(chrom_files)); i != end(chrom_files); ++i) {

      const size_t chrom_counter = chroms.size();
      vector<string> tmp_chroms, tmp_names;
      read_fasta_file(*i, tmp_names, tmp_chroms);

      chroms.resize(chrom_counter + tmp_chroms.size());
      names.resize(chrom_counter + tmp_chroms.size());
      chrom_sizes.resize(chrom_counter + tmp_chroms.size());

      for (size_t j = 0; j < tmp_chroms.size(); ++j) {
        const size_t k = chrom_counter + j;
        names[k].swap(tmp_names[j]);
        chroms[k].swap(tmp_chroms[j]);
        chrom_sizes[k] = chroms[k].size();
        chrom_lookup[names[k]] = k;
      }
    }
    vector<string> chrom_order(names);
    sort(begin(chrom_order), end(chrom_order));

    if (VERBOSE)
      cerr << "n_chroms: " << chroms.size() << endl;

    std::ifstream in(mapped_reads_file.c_str());
    if (!in)
      throw runtime_error("cannot open file: " + mapped_reads_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    size_t chrom_id = 0;
    size_t chrom_size = 0;
    GenomicRegion chrom_region; // holds chrom name for fast comparisons
    MappedRead mr;

    // this is where all the counts are accumulated
    vector<CountSet<unsigned short> > counts;

    size_t j = 0; // current chromosome

    while (in >> mr) {

      // if chrom changes, output previous results, get new one
      if (counts.empty() || !mr.r.same_chrom(chrom_region)) {

        // make sure all reads from same chrom are contiguous in the file
        if (mr.r.get_chrom() < chrom_region.get_chrom())
          throw runtime_error("chroms out of order: " + mapped_reads_file);

        if (!counts.empty()) {// should be true after first iteration
          write_output(out, chrom_order[j], chroms[chrom_id], counts, CPG_ONLY);
          ++j;
        }
        while (j < chrom_order.size() && chrom_order[j] != mr.r.get_chrom()) {
          if (VERBOSE)
            cerr << "NO_READS:\t" << chrom_order[j] << endl;
          chrom_id = get_chrom_id(chrom_order[j], chrom_lookup);
          chrom_size = chrom_sizes[chrom_id];
          chrom_region.set_chrom(chrom_order[j]);
          counts.clear();
          counts.resize(chrom_size);
          write_output(out, chrom_order[j], chroms[chrom_id], counts, CPG_ONLY);
          ++j;
        }

        if (j == chrom_order.size() || chrom_order[j] != mr.r.get_chrom())
          throw runtime_error("problem with chrom order in mapped reads");

        // move to the current chromosome
        chrom_id = get_chrom_id(chrom_order[j], chrom_lookup);
        chrom_size = chrom_sizes[chrom_id];
        chrom_region.set_chrom(chrom_order[j]);
        // reset the counts
        counts.clear();
        counts.resize(chrom_size);

        if (VERBOSE)
          cerr << "PROCESSING:\t" << chrom_order[j] << endl;
      }

      // do the work for this mapped read, depending on strand
      if (mr.r.pos_strand())
        count_states_pos(chrom_size, mr, counts);
      else count_states_neg(chrom_size, mr, counts);

    }
    if (!counts.empty()) {// should be true after first iteration
      write_output(out, chrom_order[j], chroms[chrom_id], counts, CPG_ONLY);
      ++j;
    }
    while (j < chrom_order.size()) {
      if (VERBOSE)
        cerr << "NO_READS:\t" << chrom_order[j] << endl;
      chrom_id = get_chrom_id(chrom_order[j], chrom_lookup);
      chrom_size = chrom_sizes[chrom_id];
      chrom_region.set_chrom(chrom_order[j]);
      counts.clear();
      counts.resize(chrom_size);
      write_output(out, chrom_order[j], chroms[chrom_id], counts, CPG_ONLY);
      ++j;
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
