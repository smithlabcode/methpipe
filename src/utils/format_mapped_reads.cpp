/* format_mapped_reads: a program to ensure SAM and BAM format reads
 * are conforming to expectations of methpipe software
 *
 * Copyright (C) 2020 University of Southern California and
 *                    Andrew D. Smith
 *
 * Authors: Andrew Smith
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
 */

/* The output of this program should include mapped reads that are
 * T-rich, which might require reverse-complementing sequences, and
 * switching their strand, along with any tag that indicates T-rich
 * vs. A-rich.
 *
 * Focusing only on single end reads, for each supported read mapping
 * tool, we require a means of determining whether or not the read is
 * A-rich and then changing that format to indicate T-rich.
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <stdexcept>
#include <cmath>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "htslib_wrapper.hpp"
#include "sam_record.hpp"
#include "cigar_utils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::max;
using std::min;
using std::runtime_error;
using std::unordered_map;
using std::swap;
using std::to_string;

// ADS: do we need to check that both mates are on the same strand? Or
// that they are on opposite strands?

static bool
abismal_is_a_rich(const sam_rec &aln) {
  auto the_cv_tag = find_if(begin(aln.tags), end(aln.tags),
                            [](const string &t) {
                              return t.compare (0, 3, "CV:") == 0;
                            });
  if (the_cv_tag != end(aln.tags))
    return the_cv_tag->back() == 'A';
  return false;
}

static bool
is_a_rich(const sam_rec &aln) {
  return abismal_is_a_rich(aln);
}

static bool
is_mapped(const sam_rec &aln) {
  return !check_flag(aln, samflags::read_unmapped);
}

static bool
is_primary(const sam_rec &aln) {
  return !check_flag(aln, samflags::secondary_aln);
}

static bool
is_mapped_single_end(const sam_rec &aln) {
  return is_mapped(aln) &&
    (!check_flag(aln, samflags::read_paired) ||
     check_flag(aln, samflags::mate_unmapped) ||
     aln.tlen == 0);
  // ||
  //    aln.rnext == "=");
}

inline bool
is_rc(const sam_rec &aln) {
  return check_flag(aln, samflags::read_rc);
}

static void
flip_strand(sam_rec &aln) {
  if (check_flag(aln, samflags::read_rc))
    unset_flag(aln, samflags::read_rc);
  else
    set_flag(aln, samflags::read_rc);
}

static void
flip_conversion(sam_rec &aln) {
  flip_strand(aln); // set strand to opposite of current value
  revcomp_inplace(aln.seq); // reverse complement sequence
  std::reverse(begin(aln.qual), end(aln.qual)); // and quality scores

  // ADS: assuming abismal here
  auto the_cv_tag = find_if(begin(aln.tags), end(aln.tags),
                            [](const string &t) {
                              return t.compare (0, 3, "CV:") == 0;
                            });
  if (the_cv_tag != end(aln.tags)) {
    if (abismal_is_a_rich(aln))
      the_cv_tag->back() = 'T';
    else
      the_cv_tag->back() = 'A';
  }
}

static uint32_t
get_r_end(const sam_rec &sr) {
  return sr.pos + cigar_rseq_ops(sr.cigar);
}


static size_t
merge_mates(const size_t suffix_len, const size_t range,
            const sam_rec &one, const sam_rec &two, sam_rec &merged) {

  // cerr << one << '\t' << is_rc(one) << endl
  //      << two << '\t' << is_rc(two) << endl;

  assert(is_rc(one) == false && is_rc(two) == true);

  // ADS: not sure this can be consistent across mappers
  assert(is_a_rich(one) != is_a_rich(two));

  merged = one;

  // arithmetic easier using base 0 so subtracting 1 from pos
  const int one_s = one.pos - 1;
  const int one_e = one_s + cigar_rseq_ops(one.cigar);
  const int two_s = two.pos - 1;
  const int two_e = two_s + cigar_rseq_ops(two.cigar);
  assert(one_s >= 0 && two_s >= 0);

  const int spacer = two_s - one_e;
  if (spacer >= 0) {
    /* fragments longer enough that there is space between them: this
     * size of the spacer ("_") is determined based on the reference
     * positions of the two ends, and here we assume "one" maps to
     * positive genome strand.
     *
     * left                                                         right
     * one_s                    one_e      two_s                    two_e
     * [------------end1------------]______[------------end2------------]
     */
    merged.seq = one.seq + revcomp(two.seq);
    // ADS: need to take care of soft clipping in between;
    merged.cigar = one.cigar + to_string(spacer) + "N" + two.cigar;
    merged.qual = (one.qual == "*" ? one.qual : one.qual + revcomp(two.qual));
  }
  else {
    const int head = two_s - one_s;
    if (head >= 0) {
      /* fragment longer than or equal to the read length, but shorter
       * than twice the read length: this is determined by obtaining
       * the size of the "head" in the diagram below: the portion of
       * end1 that is not within [=]. If the read maps to the positive
       * strand, this depends on the reference start of end2 minus the
       * reference start of end1. For negative strand, this is
       * reference start of end1 minus reference start of end2.
       *
       * <======= head ================>
       *
       * left                                             right
       * one_s              two_s      one_e              two_e
       * [------------end1------[======]------end2------------]
       */
      truncate_cigar_r(merged.cigar, head);
      const uint32_t one_seq_len = cigar_qseq_ops(merged.cigar);
      merged.cigar += two.cigar;
      merge_equal_neighbor_cigar_ops(merged.cigar);
      // ADS: need to take care of soft clipping in between;
      merged.seq.resize(one_seq_len);
      merged.seq += revcomp(two.seq);
      if (merged.qual != "*") {
        merged.qual.resize(one_seq_len);
        merged.qual += two.qual;
      }
    }
    else {
      /* dovetail fragments shorter than read length: this is
       * identified if the above conditions are not satisfied, but
       * there is still some overlap. The overlap will be at the 5'
       * ends of reads, which in theory shouldn't happen unless the
       * two ends are covering identical genomic intervals.
       *
       * left                                       right
       * two_s            one_s    two_e            one_e
       * [--end2----------[============]----------end1--]
       */
      const int overlap = two_e - one_s;
      if (overlap >= 0) {
        truncate_cigar_r(merged.cigar, overlap);
        const uint32_t overlap_qlen = cigar_qseq_ops(merged.cigar);
        merged.seq.resize(overlap_qlen);
        if (merged.qual != "*")
          merged.qual.resize(overlap_qlen);
      }
    }
  }
  merged.rnext = "*";
  merged.pnext = 0;
  merged.tlen = 0;

  return two_e - one_s;
}

/********Above are functions for merging pair-end reads********/

// ADS: there is a bug somewhere when a value of 0 is given for
// suffix_len
static string
remove_suffix(const string &x, const size_t suffix_len) {
  return x.size() > suffix_len ? x.substr(0, x.size() - suffix_len) : x;
}

static bool
precedes_by_more_than(const sam_rec &a, const sam_rec &b,
                      const size_t max_frag_len) {
  return (a.rname < b.rname ||
          (a.rname == b.rname && (get_r_end(a) + max_frag_len < b.pos)));
}


inline bool
bsmap_get_rc(const string &strand_tag) {
  return strand_tag.size() >= 5 && strand_tag[5] == '-';
}

inline bool
bsmap_get_a_rich(const string &richness_tag) {
  return richness_tag.size() >= 6 && richness_tag[6] == '-';
}

static void
standardize_format(const string &input_format, sam_rec &aln) {

  if (input_format == "abismal") return;

  if (input_format == "bsmap") {
    auto z_tag_itr = find_if(begin(aln.tags), end(aln.tags),
                             [](const string &t) {
                               return t.compare (0, 3, "ZS:") == 0;
                             });
    if (z_tag_itr == end(aln.tags))
      throw runtime_error("record appears to be invalid for bsmap");
    const string z_tag = *z_tag_itr;
    aln.tags.erase(z_tag_itr);

    if (bsmap_get_rc(z_tag))
      set_flag(aln, samflags::read_rc);
    else
      unset_flag(aln, samflags::read_rc);

    const bool ar = bsmap_get_a_rich(z_tag);
    if (ar)
      flip_strand(aln);
    aln.add_tag(ar ? "CV:A:A" : "CV:A:T");
  }

  // doesn't depend on mapper
  aln.qual = "*";
}


int
main(int argc, const char **argv) {

  try {

    string outfile;
    string input_format;
    int max_frag_len = 3000;
    size_t max_dangling = 5000;
    size_t suffix_len = 1;
    bool VERBOSE = false;

    const string description = "convert SAM/BAM mapped bs-seq reads "
      "to standard methpipe format";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<sam/bam-file>", 1);
    opt_parse.add_opt("format", 'f', "input format (abismal, bsmap, bismark)",
                      false, input_format);
    opt_parse.add_opt("output", 'o', "output file name",
                      false, outfile);
    opt_parse.add_opt("suff", 's', "read name suffix length (default: 1)",
                      false, suffix_len);
    opt_parse.add_opt("max-frag", 'L', "maximum allowed insert size",
                      false, max_frag_len);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
    opt_parse.set_show_defaults();
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc < 3 || opt_parse.help_requested()) {
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

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    if (VERBOSE)
      cerr << "[input file: " << mapped_reads_file << "]" << endl
           << "[output file: "
           << (outfile.empty() ? "stdout" : outfile) << "]" << endl;

    SAMReader sam_reader(mapped_reads_file);
    unordered_map<string, sam_rec> dangling_mates;

    out << sam_reader.get_header(); // includes newline

    size_t count = 0;
    sam_rec aln;

    while (sam_reader >> aln) {

      standardize_format(input_format, aln);

      if (is_mapped(aln) && is_primary(aln)) {
        if (is_mapped_single_end(aln)) {
          if (is_a_rich(aln))
            flip_conversion(aln);
          out << aln << '\n';
        }
        else { // is_mapped_paired(aln)

          const string read_name(remove_suffix(aln.qname, suffix_len));
          auto the_mate = dangling_mates.find(read_name);
          if (the_mate == end(dangling_mates)) {
            dangling_mates[read_name] = aln; // no mate seen yet
          }
          else { // found a mate

            // ADS: below is essentially a check for dovetail
            if (!is_rc(aln))
              swap(aln, the_mate->second);

            sam_rec merged;
            const int frag_len = merge_mates(suffix_len, max_frag_len,
                                             the_mate->second, aln, merged);

            if (is_a_rich(merged))
              flip_conversion(merged);

            if (frag_len <= max_frag_len)
              out << merged << '\n';
            else if (frag_len > 0)
              out << the_mate->second << '\n'
                  << aln << '\n';
            dangling_mates.erase(read_name);
          }

          // flush dangling_mates
          if (dangling_mates.size() > max_dangling) {
            unordered_map<string, sam_rec> to_keep;
            for (auto &&mates : dangling_mates)
              if (precedes_by_more_than(the_mate->second, aln, max_frag_len)) {
                if (is_a_rich(mates.second))
                  flip_conversion(mates.second);
                out << mates.second << endl;
              }
              else to_keep.insert(mates);
            swap(to_keep, dangling_mates);
          }
        }
        ++count;
      }
    }

    // flushing dangling_mates
    while (!dangling_mates.empty()) {
      if (is_a_rich(begin(dangling_mates)->second))
        flip_conversion(begin(dangling_mates)->second);
      out << begin(dangling_mates)->second << endl;
      dangling_mates.erase(begin(dangling_mates));
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
