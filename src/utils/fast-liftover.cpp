/* fast-liftover: lift over sites using index file
 *
 * Copyright (C) 2014 University of Southern California and
 *                    Andrew D. Smith
 *
 * Authors: Jenny Qu, Qiang Song
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */


/*
  Sample indexfile line:
  [T-chr] [T-start] [T-end] [S-chr]:[S-start]:[S-end]:[S-strand] [] [T-strand]
  chr21   26608683        26608684        chr1:3007015:3007016:-  0       +
*/

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "MethpipeFiles.hpp"
#include "MethpipeSite.hpp"

using std::string;
using std::ios_base;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unordered_map;
using std::runtime_error;

struct SimpleSite {
  string chrom;
  uint32_t pos;
  char strand;
  SimpleSite() {}
  SimpleSite(const string &c, const uint32_t p, const char s) :
    chrom(c), pos(p), strand(s) {}
};

void
flip_strand(SimpleSite &s) {
  if (s.strand == '-') {
    s.pos--;
    s.strand = '+';
  }
}

typedef
unordered_map<size_t, SimpleSite> liftover_index;

static void
read_index_file(const bool plus_strand, const string &index_file,
                unordered_map<string, liftover_index> &index) {

  std::ifstream in(index_file);
  if (!in)
    throw runtime_error("problem opening index file: " + index_file);

  size_t from_pos, to_pos;
  string from_chrom, to_chrom;
  string to_strand;
  MSite curr_site;
  while (in >> from_chrom >> from_pos >> to_chrom >> to_pos >> to_strand) {
    SimpleSite the_site(to_chrom, to_pos, to_strand[0]);
    if (plus_strand)
      flip_strand(the_site);
    index[from_chrom][from_pos] = the_site;
  }
}

static bool
lift_site(const unordered_map<string, liftover_index> &index,
          MSite &meth_site) {

  auto chrom_index = index.find(meth_site.chrom);
  if (chrom_index == end(index))
    return false;

  auto pos_index = chrom_index->second.find(meth_site.pos);
  if (pos_index == end(chrom_index->second))
    return false;

  meth_site.chrom = pos_index->second.chrom;
  meth_site.pos = pos_index->second.pos;
  meth_site.strand = pos_index->second.strand;
  return true;
}

int
main(int argc, const char **argv) {
  try {
    string indexfile;
    string tofile;
    string fromfile;
    string unlifted_file;

    bool VERBOSE = false;
    bool plus_strand = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "Fast liftOver-all cytosine-by strand" );
    opt_parse.add_opt("indexfile", 'i', "index file", true, indexfile);
    opt_parse.add_opt("from", 'f', "Original file", true, fromfile);
    opt_parse.add_opt("to", 't', "Output file liftovered", true, tofile);
    opt_parse.add_opt("unlifted", 'u', "(optional) File for unlifted sites",
                      false, unlifted_file);
    opt_parse.add_opt("plus-strand", 'p', "(optional) Report sites on + strand",
                      false, plus_strand);
    opt_parse.add_opt("verbose", 'v', "(optional) Print more information",
                      false, VERBOSE);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
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
    /****************** END COMMAND LINE OPTIONS *****************/

    unordered_map<string, liftover_index> index;
    if (VERBOSE)
      cerr << "[loading liftover index file " << indexfile << "]" << endl;
    read_index_file(plus_strand, indexfile, index);

    std::ifstream in(fromfile);
    if (!in)
      throw runtime_error("cannot open input file: " + fromfile);

    std::ofstream out(tofile);
    if (!out)
      throw runtime_error("cannot open output file: " + tofile);

    std::ofstream unlifted;
    if (!unlifted_file.empty())
      unlifted.open(unlifted_file.c_str());

    if (VERBOSE)
      cerr << "[lifting from: " << fromfile << " to: " << tofile << "]" << endl;

    MSite lifted, meth_site;
    while (in >> meth_site) {
      if (lift_site(index, meth_site))
        out << meth_site.tostring() << '\n';
      else if (!unlifted_file.empty())
        unlifted << meth_site.tostring() << '\n';
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
