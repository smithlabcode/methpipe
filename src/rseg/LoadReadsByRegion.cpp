/*
 * Copyright (C) 2011 University of Southern California
 *                    Andrew D Smith and Qiang Song
 * Author: Qiang Song and Andrew D. Smith
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include "LoadReadsByRegion.hpp"
#include "SortGenomicRegion.hpp"

#include "smithlab_utils.hpp"

#include "BAMFile.hpp"

#include <iomanip>
#include <fstream>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::min;
using std::max;

template <class T> bool
check_sorted(const std::vector<T> &regions, bool require_unique = false) {
  if (require_unique) {
    for (size_t i = 1; i < regions.size(); ++i)
      if (regions[i] <= regions[i - 1])
        return false;
  }
  else
    for (size_t i = 1; i < regions.size(); ++i)
      if (regions[i] < regions[i - 1])
        return false;
  return true;
}

/*************************************************
 * This function takes the names of three files (a reads file, a
 * chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegionBED(const bool VERBOSE,
                     const std::string &chroms_file,
                     const std::string &reads_file,
                     const std::string &deads_file,
                     const size_t bin_size,
                     std::vector<SimpleGenomicRegion> &bin_boundaries,
                     std::vector<double> &read_bins,
                     std::vector<double> &nondead_scales,
                     std::vector<size_t> &reset_points,
                     const size_t FRAGMENT_LEN,
                     const bool REMOVE_JACKPOT)
{
  // get the chroms
  if (VERBOSE)
    cout << "[LOADING_DATA] chromosomes" << endl;
  vector<SimpleGenomicRegion> chroms;
  ReadBEDFile(chroms_file, chroms);
  if (!check_sorted(chroms, true))
    SortGenomicRegion::sort_regions(chroms);

  // Create bins
  reset_points.push_back(0);
  for (size_t i = 0; i < chroms.size(); ++i) {
    const string chrom(chroms[i].get_chrom());
    for (size_t j = 0; j < chroms[i].get_width(); j += bin_size)
      bin_boundaries.push_back(SimpleGenomicRegion(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    cout << "[LOADING_DATA] reads"  << endl;
  read_bins.resize(bin_boundaries.size(), 0);
  std::ifstream in(reads_file.c_str());
  string line;
  size_t i = 0;
  GenomicRegion prev_gr;
  while (getline(in, line)) {
    GenomicRegion gr(line);
    if (REMOVE_JACKPOT
        && prev_gr.get_start() == gr.get_start()
        && prev_gr.get_strand() == gr.get_strand()
        && prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr < prev_gr) {
      cerr << "ERROR: reads not sorted in " << reads_file << endl;
      cerr << prev_gr << endl << gr << endl;
      exit(-1);
    }
    prev_gr = gr;

    const size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand()) {
      gr.set_start(gr.get_start() + half_len);
      gr.set_end(gr.get_start() + 1);
    }
    else {
      gr.set_end(gr.get_end() - half_len);
      gr.set_start(gr.get_end() - 1);
    }

    while (i < bin_boundaries.size() &&
           (bin_boundaries[i].get_chrom() < gr.get_chrom() ||
            (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
             bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0 &&
           (bin_boundaries[i].get_chrom() == gr.get_chrom() &&
            bin_boundaries[i].get_start() >= gr.get_end()))
      --i;

    if (i >= bin_boundaries.size() || !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins[i];
  }

  // load the dead zones
  if (VERBOSE)
    cerr << "[LOADING_DATA] deadzones"  << endl;
  nondead_scales.resize(bin_boundaries.size(), 1.0);
  in.close();
  in.open(deads_file.c_str());
  i = 0;
  while (getline(in, line)) {
    SimpleGenomicRegion gr(line);
    while (i < bin_boundaries.size() && !bin_boundaries[i].overlaps(gr)) ++i;
    while (i < bin_boundaries.size() && bin_boundaries[i].overlaps(gr)) {
      const double dead = min(bin_boundaries[i].get_end(), gr.get_end())
        - max(bin_boundaries[i].get_start(), gr.get_start());
      nondead_scales[i] -= dead / bin_size;
      ++i;
    }
    i = i > 0 ? i - 1 : 0;
  }
}

/*************************************************
 * This function takes the names of four files (reads file a, reads
 * file b, a chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegionBED(const bool VERBOSE,
                     const std::string &chroms_file,
                     const std::string &reads_file_a,
                     const std::string &reads_file_b,
                     const std::string &deads_file,
                     const size_t bin_size,
                     std::vector<SimpleGenomicRegion> &bin_boundaries,
                     std::vector<double> &read_bins_a,
                     std::vector<double> &read_bins_b,
                     std::vector<double> &nondead_scales,
                     std::vector<size_t> &reset_points,
                     const size_t FRAGMENT_LEN,
                     const bool REMOVE_JACKPOT)
{
  // get the chroms
  if (VERBOSE)
    cout << "[LOADING_DATA] chromosomes" << endl;
  vector<SimpleGenomicRegion> chroms;
  ReadBEDFile(chroms_file, chroms);
  if (!check_sorted(chroms, true))
    SortGenomicRegion::sort_regions(chroms);

  // Create bins
  reset_points.push_back(0);
  for (size_t i = 0; i < chroms.size(); ++i) {
    const string chrom(chroms[i].get_chrom());
    for (size_t j = 0; j < chroms[i].get_width(); j += bin_size)
      bin_boundaries.push_back(SimpleGenomicRegion(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    cout << "[LOADING_DATA] reads"  << endl;
  read_bins_a.resize(bin_boundaries.size(), 0);
  std::ifstream in(reads_file_a.c_str());
  string line;
  size_t i = 0;
  GenomicRegion prev_gr;
  while (getline(in, line)) {
    GenomicRegion gr(line);
    if (REMOVE_JACKPOT
        && prev_gr.get_start() == gr.get_start()
        && prev_gr.get_strand() == gr.get_strand()
        && prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr < prev_gr) {
      cerr << "ERROR: reads not sorted in " << reads_file_a << endl;
      cerr << prev_gr << endl << gr << endl;
      exit(-1);
    }
    prev_gr = gr;

    const size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand())
      {
        gr.set_start(gr.get_start() + half_len);
        gr.set_end(gr.get_start() + 1);
      }
    else
      {
        gr.set_end(gr.get_end() - half_len);
        gr.set_start(gr.get_end() - 1);
      }

    while (i < bin_boundaries.size()
           && (bin_boundaries[i].get_chrom() < gr.get_chrom()
               || (bin_boundaries[i].get_chrom() == gr.get_chrom()
                   && bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0
           && (bin_boundaries[i].get_chrom() == gr.get_chrom()
               && bin_boundaries[i].get_start() >= gr.get_end()))
      --i;

    if (i >= bin_boundaries.size() || !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins_a[i];
  }

  // Load the reads, tabulating counts in bins
  read_bins_b.resize(bin_boundaries.size(), 0);
  in.close();
  in.open(reads_file_b.c_str());
  i = 0;
  prev_gr = GenomicRegion();
  while (getline(in, line)) {
    GenomicRegion gr(line);
    if (REMOVE_JACKPOT
        && prev_gr.get_start() == gr.get_start()
        && prev_gr.get_strand() == gr.get_strand()
        && prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr < prev_gr) {
      cerr << "ERROR: reads not sorted in " << reads_file_b << endl;
      cerr << prev_gr << endl << gr << endl;
      exit(-1);
    }
    prev_gr = gr;

    const size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand())
      {
        gr.set_start(gr.get_start() + half_len);
        gr.set_end(gr.get_start() + 1);
      }
    else
      {
        gr.set_end(gr.get_end() - half_len);
        gr.set_start(gr.get_end() - 1);
      }

    while (i < bin_boundaries.size()
           && (bin_boundaries[i].get_chrom() < gr.get_chrom()
               || (bin_boundaries[i].get_chrom() == gr.get_chrom()
                   && bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0
           && (bin_boundaries[i].get_chrom() == gr.get_chrom()
               && bin_boundaries[i].get_start() >= gr.get_end()))
      --i;
    if (i >= bin_boundaries.size() || !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins_b[i];
  }

  // load the dead zones
  if (VERBOSE)
    cout << "[LOADING_DATA] deadzones"  << endl;
  nondead_scales.resize(bin_boundaries.size(), 1.0);
  in.close();
  in.open(deads_file.c_str());
  i = 0;
  while (getline(in, line)) {
    SimpleGenomicRegion gr(line);
    while (i < bin_boundaries.size() && !bin_boundaries[i].overlaps(gr)) ++i;
    while (i < bin_boundaries.size() && bin_boundaries[i].overlaps(gr)) {
      const double dead = min(bin_boundaries[i].get_end(), gr.get_end())
        - max(bin_boundaries[i].get_start(), gr.get_start());
      nondead_scales[i] -= dead / bin_size;
      ++i;
    }
    i = i > 0 ? i - 1 : 0;
  }
}

/*************************************************
 * This function takes the names of three files (a reads file, a
 * chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegionBAM(const bool VERBOSE,
                     const std::string &chroms_file,
                     const std::string &reads_file,
                     const std::string &deads_file,
                     const size_t bin_size,
                     std::vector<SimpleGenomicRegion> &bin_boundaries,
                     std::vector<double> &read_bins,
                     std::vector<double> &nondead_scales,
                     std::vector<size_t> &reset_points,
                     const size_t FRAGMENT_LEN,
                     const bool REMOVE_JACKPOT)
{
  // get the chroms
  if (VERBOSE)
    cout << "[LOADING_DATA] chromosomes" << endl;
  vector<SimpleGenomicRegion> chroms;
  ReadBEDFile(chroms_file, chroms);
  if (!check_sorted(chroms, true))
    SortGenomicRegion::sort_regions(chroms);

  // Create bins
  reset_points.push_back(0);
  for (size_t i = 0; i < chroms.size(); ++i) {
    const string chrom(chroms[i].get_chrom());
    for (size_t j = 0; j < chroms[i].get_width(); j += bin_size)
      bin_boundaries.push_back(SimpleGenomicRegion(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    cout << "[LOADING_DATA] reads"  << endl;
  read_bins.resize(bin_boundaries.size(), 0);
  BAMFile in(reads_file);
  size_t i = 0;
  GenomicRegion prev_gr, gr;
  while ((in >> gr).good()) {
    if (REMOVE_JACKPOT
        && prev_gr.get_start() == gr.get_start()
        && prev_gr.get_strand() == gr.get_strand()
        && prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr.get_chrom() < prev_gr.get_chrom()
        || gr.get_start() < prev_gr.get_start()
        || gr.get_end() < prev_gr.get_end()
        || gr.get_strand() < prev_gr.get_strand()) {
      cerr << "ERROR: reads not sorted in " << reads_file << endl;
      cerr << prev_gr << endl << gr << endl;
      exit(-1);
    }
    prev_gr = gr;

    const size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand())
      {
        gr.set_start(gr.get_start() + half_len);
        gr.set_end(gr.get_start() + 1);
      }
    else
      {
        gr.set_end(gr.get_end() - half_len);
        gr.set_start(gr.get_end() - 1);
      }

    while (i < bin_boundaries.size()
           && (bin_boundaries[i].get_chrom() < gr.get_chrom()
               || (bin_boundaries[i].get_chrom() == gr.get_chrom()
                   && bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0
           && (bin_boundaries[i].get_chrom() == gr.get_chrom()
               && bin_boundaries[i].get_start() >= gr.get_end()))
      --i;
    if (i >= bin_boundaries.size() || !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins[i];
  }

  // load the dead zones
  if (VERBOSE)
    cerr << "[LOADING_DATA] deadzones"  << endl;
  nondead_scales.resize(bin_boundaries.size(), 1.0);
  std::ifstream dead_in(deads_file.c_str());
  i = 0;
  string line;
  while (getline(dead_in, line)) {
    SimpleGenomicRegion gr(line);
    while (i < bin_boundaries.size() && !bin_boundaries[i].overlaps(gr)) ++i;
    while (i < bin_boundaries.size() && bin_boundaries[i].overlaps(gr)) {
      const double dead = min(bin_boundaries[i].get_end(), gr.get_end())
        - max(bin_boundaries[i].get_start(), gr.get_start());
      nondead_scales[i] -= dead / bin_size;
      ++i;
    }
    i = i > 0 ? i - 1 : 0;
  }
}

/*************************************************
 * This function takes the names of four files (reads file a, reads
 * file b, a chromosome file, and a dead zones file [possibly empty])
 */
static void
LoadReadsByRegionBAM(const bool VERBOSE,
                     const std::string &chroms_file,
                     const std::string &reads_file_a,
                     const std::string &reads_file_b,
                     const std::string &deads_file,
                     const size_t bin_size,
                     std::vector<SimpleGenomicRegion> &bin_boundaries,
                     std::vector<double> &read_bins_a,
                     std::vector<double> &read_bins_b,
                     std::vector<double> &nondead_scales,
                     std::vector<size_t> &reset_points,
                     const size_t FRAGMENT_LEN,
                     const bool REMOVE_JACKPOT)
{
  // get the chroms
  if (VERBOSE)
    cout << "[LOADING_DATA] chromosomes" << endl;
  vector<SimpleGenomicRegion> chroms;
  ReadBEDFile(chroms_file, chroms);
  if (!check_sorted(chroms, true))
    SortGenomicRegion::sort_regions(chroms);

  // Create bins
  reset_points.push_back(0);
  for (size_t i = 0; i < chroms.size(); ++i) {
    const string chrom(chroms[i].get_chrom());
    for (size_t j = 0; j < chroms[i].get_width(); j += bin_size)
      bin_boundaries.push_back(SimpleGenomicRegion(chrom, j, j + bin_size));
    if (bin_boundaries.size() > reset_points.back())
      reset_points.push_back(bin_boundaries.size());
  }

  // Load the reads, tabulating counts in bins
  if (VERBOSE)
    cout << "[LOADING_DATA] reads"  << endl;
  read_bins_a.resize(bin_boundaries.size(), 0);
  BAMFile in(reads_file_a);
  size_t i = 0;
  GenomicRegion prev_gr, gr;
  while ((in >> gr).good()) {
    if (REMOVE_JACKPOT
        && prev_gr.get_start() == gr.get_start()
        && prev_gr.get_strand() == gr.get_strand()
        && prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr.get_chrom() < prev_gr.get_chrom()
        || gr.get_start() < prev_gr.get_start()
        || gr.get_end() < prev_gr.get_end()
        || gr.get_strand() < prev_gr.get_strand()) {
      cerr << "ERROR: reads not sorted in " << reads_file_a << endl;
      cerr << prev_gr << endl << gr << endl;
      exit(-1);
    }
    prev_gr = gr;

    const size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand())
      {
        gr.set_start(gr.get_start() + half_len);
        gr.set_end(gr.get_start() + 1);
      }
    else
      {
        gr.set_end(gr.get_end() - half_len);
        gr.set_start(gr.get_end() - 1);
      }

    while (i < bin_boundaries.size()
           && (bin_boundaries[i].get_chrom() < gr.get_chrom()
               || (bin_boundaries[i].get_chrom() == gr.get_chrom()
                   && bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0
           && (bin_boundaries[i].get_chrom() == gr.get_chrom()
               && bin_boundaries[i].get_start() >= gr.get_end()))
      --i;
    if (i >= bin_boundaries.size() || !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins_a[i];
  }

  // Load the reads, tabulating counts in bins
  read_bins_b.resize(bin_boundaries.size(), 0);
  in.close();
  in.open(reads_file_b);
  i = 0;
  prev_gr = GenomicRegion();
  while ((in >> gr).good()) {
    if (REMOVE_JACKPOT
        && prev_gr.get_start() == gr.get_start()
        && prev_gr.get_strand() == gr.get_strand()
        && prev_gr.get_chrom() == gr.get_chrom())
      continue;
    if (gr.get_chrom() < prev_gr.get_chrom()
        || gr.get_start() < prev_gr.get_start()
        || gr.get_end() < prev_gr.get_end()
        || gr.get_strand() < prev_gr.get_strand()) {
      cerr << "ERROR: reads not sorted in " << reads_file_b << endl;
      cerr << prev_gr << endl << gr << endl;
      exit(-1);
    }
    prev_gr = gr;

    const size_t half_len = FRAGMENT_LEN / 2;
    if (gr.pos_strand())
      {
        gr.set_start(gr.get_start() + half_len);
        gr.set_end(gr.get_start() + 1);
      }
    else
      {
        gr.set_end(gr.get_end() - half_len);
        gr.set_start(gr.get_end() - 1);
      }

    while (i < bin_boundaries.size()
           && (bin_boundaries[i].get_chrom() < gr.get_chrom()
               || (bin_boundaries[i].get_chrom() == gr.get_chrom()
                   && bin_boundaries[i].get_end() <= gr.get_start())))
      ++i;

    // adjust for jump caused by negative strand reads
    while (i >= 0
           && (bin_boundaries[i].get_chrom() == gr.get_chrom()
               && bin_boundaries[i].get_start() >= gr.get_end()))
      --i;
    if (i >= bin_boundaries.size() || !bin_boundaries[i].contains(gr))
      continue;
    ++read_bins_b[i];
  }

  // load the dead zones
  if (VERBOSE)
    cout << "[LOADING_DATA] deadzones"  << endl;
  nondead_scales.resize(bin_boundaries.size(), 1.0);
  std::ifstream dead_in(deads_file.c_str());
  i = 0;
  string line;
  while (getline(dead_in, line)) {
    SimpleGenomicRegion gr(line);
    while (i < bin_boundaries.size() && !bin_boundaries[i].overlaps(gr)) ++i;
    while (i < bin_boundaries.size() && bin_boundaries[i].overlaps(gr)) {
      const double dead = min(bin_boundaries[i].get_end(), gr.get_end())
        - max(bin_boundaries[i].get_start(), gr.get_start());
      nondead_scales[i] -= dead / bin_size;
      ++i;
    }
    i = i > 0 ? i - 1 : 0;
  }
}


/////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
void
LoadReadsByRegion(const bool VERBOSE,
                  const std::string &chroms_file,
                  const std::string &reads_file,
                  const std::string &deads_file,
                  const size_t bin_size,
                  std::vector<SimpleGenomicRegion> &bin_boundaries,
                  std::vector<double> &read_bins,
                  std::vector<double> &nondead_scales,
                  std::vector<size_t> &reset_points,
                  const size_t FRAGMENT_LEN,
                  const bool BAM_FORMAT,
                  const bool REMOVE_JACKPOT) {
  if (BAM_FORMAT)
    LoadReadsByRegionBAM(VERBOSE, chroms_file, reads_file, deads_file,
                         bin_size, bin_boundaries, read_bins,
                         nondead_scales, reset_points,
                         FRAGMENT_LEN, REMOVE_JACKPOT);
  else
    LoadReadsByRegionBED(VERBOSE, chroms_file, reads_file, deads_file,
                         bin_size, bin_boundaries, read_bins,
                         nondead_scales, reset_points,
                         FRAGMENT_LEN, REMOVE_JACKPOT);
}


void
LoadReadsByRegion(const bool VERBOSE,
                  const std::string &chroms_file,
                  const std::string &reads_file_a,
                  const std::string &reads_file_b,
                  const std::string &deads_file,
                  const size_t bin_size,
                  std::vector<SimpleGenomicRegion> &bin_boundaries,
                  std::vector<double> &read_bins_a,
                  std::vector<double> &read_bins_b,
                  std::vector<double> &nondead_scales,
                  std::vector<size_t> &reset_points,
                  const size_t FRAGMENT_LEN,
                  const bool BAM_FORMAT,
                  const bool REMOVE_JACKPOT) {
  if (BAM_FORMAT)
    LoadReadsByRegionBAM(VERBOSE, chroms_file, reads_file_a, reads_file_b,
                         deads_file, bin_size, bin_boundaries,
                         read_bins_a, read_bins_b, nondead_scales,
                         reset_points, FRAGMENT_LEN, REMOVE_JACKPOT);
  else
    LoadReadsByRegionBED(VERBOSE, chroms_file, reads_file_a, reads_file_b,
                         deads_file, bin_size, bin_boundaries,
                         read_bins_a, read_bins_b, nondead_scales,
                         reset_points, FRAGMENT_LEN, REMOVE_JACKPOT);
}


