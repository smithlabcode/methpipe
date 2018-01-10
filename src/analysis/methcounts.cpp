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
using std::max;
using std::accumulate;
using std::unordered_map;

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
get_methylation_context_tag(const string &s, const size_t pos) {
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
count_states_pos(const string &chrom, const MappedRead &r,
                 vector<CountSet<count_type> > &counts) {
  const size_t width = r.r.get_width();

  size_t position = r.r.get_start();
  assert(position + width <= chrom.length());
  for (size_t i = 0; i < width; ++i, ++position)
    counts[position].add_count_pos(r.seq[i]);
}


template <class count_type>
static void
count_states_neg(const string &chrom, const MappedRead &r,
                 vector<CountSet<count_type> > &counts) {
  const size_t width = r.r.get_width();

  size_t position = r.r.get_start() + width - 1;
  assert(position < chrom.length());
  for (size_t i = 0; i < width; ++i, --position)
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
      const string tag = get_methylation_context_tag(chrom, i) +
        (has_mutated(base, counts[i]) ? "x" : "");
      if (!CPG_ONLY || (!tag.compare("CpG")||!tag.compare("CpGx"))) {
        methpipe::write_site(out, chrom_name, i,
                             (is_cytosine(base) ? "+" : "-"),
                             tag, meth, converted + unconverted);
      }
    }
  }
}


typedef unordered_map<string, string> chrom_file_map;
static void
get_chrom(const MappedRead &mr,
          const chrom_file_map &chrom_files,
          GenomicRegion &chrom_region, string &chrom) {

  const chrom_file_map::const_iterator fn(chrom_files.find(mr.r.get_chrom()));
  if (fn == chrom_files.end())
    throw SMITHLABException("could not find chrom: " + mr.r.get_chrom());

  chrom.clear();
  read_fasta_file(fn->second, mr.r.get_chrom(), chrom);
  if (chrom.empty())
    throw SMITHLABException("could not find chrom: " + mr.r.get_chrom());

  chrom_region.set_chrom(mr.r.get_chrom());
}

int
main(int argc, const char **argv) {

  try {

    static const double gigs_per_base =
      (1.0 + sizeof(CountSet<unsigned short>))/(1073741824.0);

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
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested() || leftover_args.size() != 1) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (!outfile.empty() && !is_valid_output_file(outfile))
      throw SMITHLABException("bad output file: " + outfile);

    chrom_file_map chrom_files;
    identify_and_read_chromosomes(chrom_file, fasta_suffix, chrom_files);
    if (VERBOSE)
      cerr << "CHROMS_FOUND=" << chrom_files.size() << endl;

    std::ifstream in(mapped_reads_file.c_str());
    if (!in)
      throw SMITHLABException("cannot open file: " + mapped_reads_file);

    // this is where all the counts are accumulated
    vector<CountSet<unsigned short> > counts;

    string chrom; // holds the current chromosome being processed
    GenomicRegion chrom_region; // holds chrom name for fast comparisons

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    MappedRead mr;
    while (in >> mr) {

      // if chrom changes, output previous results, get new one
      if (chrom.empty() || !mr.r.same_chrom(chrom_region)) {

        // make sure all reads from same chrom are contiguous in the file
        if (mr.r.get_chrom() < chrom_region.get_chrom())
          throw SMITHLABException("chroms out of order: " + mapped_reads_file);

        if (!counts.empty()) // if we have results, output them
          write_output(out, chrom_region.get_chrom(), chrom, counts, CPG_ONLY);

        // load the new chromosome and reset the counts
        get_chrom(mr, chrom_files, chrom_region, chrom);
        if (VERBOSE)
          cerr << "PROCESSING:\t" << chrom_region.get_chrom() << '\t'
               << "(REQD MEM = "
               << std::setprecision(3)
               << chrom.length()*gigs_per_base << "GB)" << endl;

        counts.clear();
        counts.resize(chrom.size());
      }

      // do the work for this mapped read, depending on strand
      if (mr.r.pos_strand())
        count_states_pos(chrom, mr, counts);
      else count_states_neg(chrom, mr, counts);
    }
    // ALWAYS output the chromosome, even if all sites are uncovered.
    write_output(out, chrom_region.get_chrom(), chrom, counts, CPG_ONLY);
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
