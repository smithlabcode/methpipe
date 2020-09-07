/* format_reads: a program to ensure SAM and BAM format reads
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
#include <sstream>

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
using std::ostringstream;

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

bool
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


static bool
are_mates(const sam_rec &one, const sam_rec &two) {
  return ((one.rnext == "=" && two.rnext == "=") ||
          (one.rnext == two.qname)) &&
           one.pnext == two.pos &&
           two.pnext == one.pos;
}

static bool
are_opposite_strands(const sam_rec &one, const sam_rec &two) {
  return (check_flag(one, samflags::template_first) !=
          check_flag(two, samflags::template_first)) &&
         (check_flag(one, samflags::template_last) !=
          check_flag(two, samflags::template_last));
}

static size_t
merge_mates(const size_t range,
            const sam_rec &one, const sam_rec &two, sam_rec &merged) {

  if (!are_mates(one, two)) {
    return -std::numeric_limits<int>::max();
  }

  // assert(is_rc(one) == false && is_rc(two) == true);
  if (!are_opposite_strands(one, two)) {
    ostringstream fail;
    fail << "Reads below are not in opposite strands\n";
    fail << one << '\n';
    fail << two << '\n';
    throw runtime_error(fail.str());
  }

  // ADS: not sure this can be consistent across mappers
  // GS: true after standardization
  // assert(is_a_rich(one) != is_a_rich(two));
  if (is_a_rich(one) == is_a_rich(two)) {
    ostringstream fail;
    fail << "Reads below do not have the same richness:\n";
    fail << one << '\n';
    fail << two << '\n';
    throw runtime_error(fail.str());
  }

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
    merged.seq += revcomp(two.seq);
    // ADS: need to take care of soft clipping in between;
    merged.cigar += to_string(spacer) + "N";
    merged.cigar += two.cigar;
    // merged.qual = (one.qual == "*" ? one.qual : one.qual + revcomp(two.qual));
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
      // if (merged.qual != "*") {
      //   merged.qual.resize(one_seq_len);
      //   merged.qual += two.qual;
      // }
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
        // if (merged.qual != "*")
        //   merged.qual.resize(overlap_qlen);
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
remove_suff(const string &x, const size_t suffix_len) {
  return x.size() > suffix_len ? x.substr(0, x.size() - suffix_len) : x;
}

bool
bsmap_get_rc(const string &strand_tag) {
  return strand_tag.size() > 5 && strand_tag[5] == '-';
}

bool
bsmap_get_a_rich(const string &richness_tag) {
  return richness_tag.size() > 6 && richness_tag[6] == '-';
}

bool
bismark_get_a_rich(const string &richness_tag) {
  return richness_tag.size() > 5 && richness_tag[5] == 'G';
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

    aln.add_tag(bsmap_get_a_rich(z_tag) ? "CV:A:A" : "CV:A:T");
    if (is_rc(aln)) revcomp_inplace(aln.seq);
  }

  if (input_format == "bismark") {
    // remove everything after _ in read name
    aln.qname = aln.qname.substr(0, aln.qname.find_first_of("_"));
    auto xr_tag_itr = find_if(begin(aln.tags), end(aln.tags),
                             [](const string &t) {
                               return t.compare (0, 3, "XR:") == 0;
                             }),
         nm_tag_itr = find_if(begin(aln.tags), end(aln.tags),
                             [](const string &t) {
                               return t.compare (0, 3, "NM:") == 0;
                             });

    if (xr_tag_itr == end(aln.tags))
      throw runtime_error("record appears to be invalid for bismark");

    const string xr_tag = *xr_tag_itr,
                 nm_tag = *nm_tag_itr;

    aln.tags.clear();
    aln.add_tag(nm_tag);
    aln.add_tag(bismark_get_a_rich(xr_tag) ? "CV:A:A" : "CV:A:T");

    if (is_rc(aln)) revcomp_inplace(aln.seq);
  }

  // doesn't depend on mapper
  aln.qual = "*";
}


int
main(int argc, const char **argv) {

  try {

    string outfile;
    string input_format;
    int max_frag_len = 10000;
    size_t buff_size = 10000;
    size_t suff_len = 1;
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
                      false, suff_len);
    opt_parse.add_opt("max-frag", 'L', "maximum allowed insert size",
                      false, max_frag_len);
    opt_parse.add_opt("buf-size", 'B', "maximum buffer size",
                      false, buff_size);
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
    std::ifstream in(mapped_reads_file);
    if (!in)
      throw std::runtime_error("problem with input file: " + mapped_reads_file);

    out << sam_reader.get_header(); // includes newline

    size_t count_a = 0, count_b = 0;
    sam_rec aln;

    vector<sam_rec> buffer(buff_size); // keeps previous records
    unordered_map<string, uint32_t> mate_lookup; // allows them to be accessed

    while (sam_reader >> aln) {
      standardize_format(input_format, aln);

      const string read_name(remove_suff(aln.qname, suff_len));
      auto the_mate = mate_lookup.find(read_name);
      if (the_mate == end(mate_lookup)) { // add solo end to buffer
        const size_t real_idx = count_a % buff_size;
        buffer[real_idx] = std::move(aln);
        mate_lookup[read_name] = real_idx;
        ++count_a;
      }
      else { // found a mate; attempt to merge
        const size_t mate_idx = the_mate->second;

        if (!is_rc(aln)) // essentially check for dovetail
          swap(buffer[mate_idx], aln);

        sam_rec merged;

        const int frag_len =
          merge_mates(max_frag_len, buffer[mate_idx], aln, merged);

        if (frag_len > 0 && frag_len < max_frag_len) {
          swap(buffer[mate_idx], merged);
          // nothing to do for mate_lookup
        }
        else { // if (frag_len > 0)
          // leave "mate" alone, but add current read
          const size_t real_idx = count_a % buff_size;
          swap(buffer[real_idx], aln);
          // no need to index this one; mate already incompatible
          ++count_a;
        }
      }

      if ((count_a % buff_size) == (count_b % buff_size)) {
        const size_t real_idx = count_b % buff_size;
        mate_lookup.erase(remove_suff(buffer[real_idx].qname, suff_len));
        if (is_a_rich(buffer[real_idx]))
          flip_conversion(buffer[real_idx]);
        out << buffer[real_idx] << '\n';
        ++count_b;
      }
    }

    for (; count_b < count_a; ++count_b) { // no need to erase names
      const size_t real_idx = count_b % buff_size;
      if (is_a_rich(buffer[real_idx]))
        flip_conversion(buffer[real_idx]);
      out << buffer[real_idx] << '\n';
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
