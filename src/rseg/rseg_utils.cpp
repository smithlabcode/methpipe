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

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <utility>
#include <numeric>

#include <cmath>
#include <cassert>
#include <cstdlib>

#include "GenomicRegion.hpp"
#include "rseg_utils.hpp"
#include "Distro.hpp"
#include "SplitDistro.hpp"
#include "numerical_utils.hpp"


using std::vector;
using std::string;
using std::cin;
using std::cerr;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::max;
using std::min;

void
pick_training_sample(const vector<double> &read_bins,
                     const vector<double> &read_bins_a,
                     const vector<double> &read_bins_b,
                     const vector<double> &scales,
                     const vector<size_t> &reset_points,
                     const size_t training_sample_size,
                     vector<double> &read_bins_sample,
                     vector<double> &read_bins_a_sample,
                     vector<double> &read_bins_b_sample,
                     vector<double> &scales_sample,
                     vector<size_t> &reset_points_sample)
{
  vector<size_t> reset_points_control;
  for (size_t i = 0; i < reset_points.size(); ++i)
    if (reset_points[i] < training_sample_size)
      {
	reset_points_sample.push_back(reset_points[i]);
      }
    else
      {
	reset_points_sample.push_back(training_sample_size);
	break;
      }

  const size_t sample_sz = reset_points_sample.back();

  std::copy(read_bins.begin(), read_bins.begin() + sample_sz,
	    std::back_inserter(read_bins_sample));
  std::copy(read_bins_a.begin(), read_bins_a.begin() + sample_sz,
	    std::back_inserter(read_bins_a_sample));
  std::copy(read_bins_b.begin(), read_bins_b.begin() + sample_sz,
	    std::back_inserter(read_bins_b_sample));
  std::copy(scales.begin(), scales.begin() + sample_sz,
	    std::back_inserter(scales_sample));
}

void
clear_training_sample(vector<double> &read_bins_sample,
		      vector<double> &read_bins_a_sample,
		      vector<double> &read_bins_b_sample,
		      vector<double> &scales_sample,
		      vector<size_t> &reset_points_sample)
{
    vector<double>().swap(read_bins_sample);
    vector<double>().swap(read_bins_a_sample);
    vector<double>().swap(read_bins_b_sample);
    vector<double>().swap(scales_sample);
    vector<size_t>().swap(reset_points_sample);
}

void
set_transitions(const size_t bin_size, const double fg_size,
		const vector<double> &mixing,
		const bool VERBOSE,
		vector<double> &start_trans,
		vector<vector<double> > &trans,
		vector<double> &end_trans)
{
  start_trans.resize(3, 0);
  start_trans[0] = mixing[0];
  start_trans[1] = mixing[1];
  start_trans[2] = mixing[2];

  end_trans.resize(3, 0);
  end_trans[0] = 1e-10;
  end_trans[1] = 1e-10;
  end_trans[2] = 1e-10;

  trans.resize(3, vector<double>(3, 0));

  double fg_bin_n = fg_size / bin_size;
  if (fg_bin_n <= 1)
    {
      cerr << "\n[Warning] RSEG may not work as expected. "
	   << "The expected differential domain size is smaller than bin size.\n" << endl;
      fg_bin_n = 2;
    }
  trans[0][0] = 1 - 1.0 / fg_bin_n;
  trans[0][1] = 1.0 / fg_bin_n * mixing[1] / (mixing[1] + mixing[2]);
  trans[0][2] = 1.0 / fg_bin_n * mixing[2] / (mixing[1] + mixing[2]);

  // assumming all fg and bg domain are sperated by middle state regions
  double mid_bin_n = mixing[1] * fg_bin_n / (2 * mixing[0]);
  if (mid_bin_n <= 1)
    {
      cerr << "\n[Warning] RSEG may not work as expected. "
	   << "The expected basal domain size is smaller than bin size.\n" << endl;
      mid_bin_n = 2;
    }
  trans[1][1] = 1 - 1.0 / mid_bin_n;
  trans[1][0] = 1.0 / mid_bin_n * mixing[0] / (mixing[0] + mixing[2]);
  trans[1][2] = 1.0 / mid_bin_n * mixing[2] / (mixing[0] + mixing[2]);

  double bg_bin_n = fg_bin_n;  // assuming symmetry
  trans[2][2] = 1 - 1.0 / bg_bin_n;
  trans[2][0] = 1.0 / bg_bin_n * mixing[0] / (mixing[0] + mixing[1]);
  trans[2][1] = 1.0 / bg_bin_n * mixing[1] / (mixing[0] + mixing[1]);

  assert(start_trans[0] > 0 && start_trans[1] > 0 && start_trans[2] > 0);
  assert(end_trans[0] > 0 && end_trans[1] > 0 && end_trans[2] > 0);

  assert(trans[0][0] > 0 && trans[0][1] > 0 && trans[0][2] > 0 &&
	 trans[1][0] > 0 && trans[1][1] > 0 && trans[1][2] > 0 &&
	 trans[2][0] > 0 && trans[2][1] > 0 && trans[2][2] > 0);
}

void
set_transitions(const size_t bin_size, const double fg_size,
		const double mixing, const bool VERBOSE,
		vector<double> &start_trans,
		vector<vector<double> > &trans,
		vector<double> &end_trans)
{
  start_trans.resize(2, 0);
  start_trans[0] = mixing;
  start_trans[1] = 1 - mixing;

  end_trans.resize(2, 0);
  end_trans[0] = 1e-10;
  end_trans[1] = 1e-10;

  trans.resize(2, vector<double>(2, 0));
  // foreground
  double fg_bin_n = fg_size / bin_size;
  if (fg_bin_n <= 1)
    {
      cerr << "[Warning] RSEG may not work as expected. "
	   << "The expected foreground domain size is smaller than bin size." << endl;
      fg_bin_n = 2;
    }
  trans[0][1] = 1.0 / fg_bin_n;
  trans[0][0] = 1 - trans[0][1];

  double bg_bin_n = fg_bin_n * (1 - mixing) / mixing;
  if (bg_bin_n <= 1)
    {
      cerr << "\n[Warning] RSEG may not work as expected."
	   << "The expected background domain size is smaller than bin size.\n" << endl;
      bg_bin_n = 2;
    }
  trans[1][0] = 1.0 / bg_bin_n;
  trans[1][1] = 1 - trans[1][0];

  assert(start_trans[0] > 0 && start_trans[1] > 0);
  assert(end_trans[0] > 0 && end_trans[1] > 0);
  assert(trans[0][1] > 0 && trans[0][0] > 0 &&
	 trans[1][0] > 0 && trans[1][1] > 0);
}

void
report_final_values(const vector<Distro> &distros,
		    const vector<double> &start_trans,
		    const vector<vector<double> > &trans,
		    const vector<double> &end_trans)
{
  cout << "FINAL ESTIMATES" << endl
       << "---------------" << endl;

  cout << "Emission distributions" << endl;
  for (size_t i = 0; i < distros.size(); ++i)
    cout << "State " << i << ":\t" << distros[i] << endl;

  cout << "Expected sizes" << endl;
  for (size_t i = 0; i < trans.size(); ++i)
    cout << "State " << i << ":\t" << 1 / (1 - trans[i][i]) << endl;

  cout << "Transition probabilities" << endl;
  for (size_t i = 0; i < trans.size(); ++i)
    {
      for (size_t j = 0; j < trans[i].size(); ++j)
	cout << trans[i][j] << "\t";
      cout << endl;
    }
}

void
report_final_values(const vector<SplitDistro> &distros,
		    const vector<double> &start_trans,
		    const vector<vector<double> > &trans,
		    const vector<double> &end_trans)
{
  cout << "FINAL ESTIMATES" << endl
       << "---------------" << endl;

  cout << "Emission distributions" << endl;
  for (size_t i = 0; i < distros.size(); ++i)
    cout << "State " << i << ":\t" << distros[i] << endl;

  cout << "Expected sizes" << endl;
  for (size_t i = 0; i < trans.size(); ++i)
    cout << "State " << i << ":\t" << 1 / (1 - trans[i][i]) << endl;

  cout << "Transition probabilities" << endl;
  for (size_t i = 0; i < trans.size(); ++i)
    {
      for (size_t j = 0; j < trans[i].size(); ++j)
	cout << trans[i][j] << "\t";
      cout << endl;
    }
}

// The slow version with system call
void
chk_and_mk_dirs(const string & path)
{
  int s = system(("mkdir -p " + path).c_str());
  if (s != 0)
    {
      cerr << "Cannot access or create the directory " + path << endl;
      exit(s);
    }
}


void
write_read_counts_by_bin(const vector< vector<SimpleGenomicRegion> > &bin_boundaries,
                         const vector<double> &read_bins,
                         const vector<double> &scales,
                         const vector<bool> &classes,
                         const string &file_name,
                         const bool VERBOSE)
{
  std::ofstream outf(file_name.c_str());

  if (!outf.good())
    {
      cerr << "Error: cannot open file " + file_name << endl;
      exit(-1);
    }

  size_t k = 0;
  for (size_t i = 0; i < bin_boundaries.size(); ++i)
    for (size_t j = 0; j < bin_boundaries[i].size(); ++j)
      {
	outf << bin_boundaries[i][j] << "\t"
	     << read_bins[k] << "\t"
	     << scales[k] << "\t"
	     << classes[k] << std::endl;
	++k;
      }

  if (!outf.good())
    {
      cerr << "Error: writing file " + file_name << endl;
      exit(-1);
    }

  outf.close();
  if (VERBOSE)
    cout << "Read count file: "
	 << file_name << std::endl;
}

void
write_read_counts_by_bin(const vector< vector<SimpleGenomicRegion> > &bin_boundaries,
                         const vector<double> &read_bins,
                         const vector<double> &read_bins_a,
                         const vector<double> &read_bins_b,
                         const vector<bool> &classes,
                         const string &file_name,
                         const bool VERBOSE)
{
  std::ofstream outf(file_name.c_str());

  if (!outf.good())
    {
      cerr << "Error: cannot open file " + file_name << endl;
      exit(-1);
    }

  size_t k = 0;
  for (size_t i = 0; i < bin_boundaries.size(); ++i)
    for (size_t j = 0; j < bin_boundaries[i].size(); ++j)
      {
	outf << bin_boundaries[i][j] << "\t"
	     << read_bins[k] << "\t"
	     << read_bins_a[k] << "\t"
	     << read_bins_b[k] << "\t"
	     << classes[k] << std::endl;
	++k;
      }

  if (!outf.good())
    {
      cerr << "Error: writing file " + file_name << endl;
      exit(-1);
    }

  outf.close();
  if (VERBOSE)
    cout << "Read count file: "
	 << file_name << std::endl;
}

void
write_read_counts_by_bin(const vector< vector<SimpleGenomicRegion> > &bin_boundaries,
                         const vector<double> &read_bins,
                         const vector<double> &read_bins_a,
                         const vector<double> &read_bins_b,
                         const vector<size_t> &classes,
                         const string &file_name,
                         const bool VERBOSE)
{
  std::ofstream outf(file_name.c_str());

  if (!outf.good())
    {
      cerr << "Error: cannot open file " + file_name << endl;
      exit(-1);
    }

  size_t k = 0;
  for (size_t i = 0; i < bin_boundaries.size(); ++i)
    for (size_t j = 0; j < bin_boundaries[i].size(); ++j)
      {
	outf << bin_boundaries[i][j] << "\t"
	     << read_bins[k] << "\t"
	     << read_bins_a[k] << "\t"
	     << read_bins_b[k] << "\t"
	     << classes[k] << std::endl;
	++k;
      }

  if (!outf.good())
    {
      cerr << "Error: writing file " + file_name << endl;
      exit(-1);
    }

  outf.close();
  if (VERBOSE)
    cout << "Read count file: "
	 << file_name << std::endl;
}

string
strip_path_and_bed_suffix(const string &full_path)
{
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else ++start;
  size_t end = full_path.find_last_of('.');
  if (end == string::npos)
    end = full_path.length();
  return full_path.substr(start, end - start);
}

void
write_wigfile(const vector<vector<double> > &scores,
	      const vector<vector<SimpleGenomicRegion> > &bin_bounds,
	      const string &wigfile_name)
{
  std::ofstream wigout(wigfile_name.c_str());

  if (!wigout.good())
    {
      cerr << "write_bed_file: cannot open " << wigfile_name << endl;
      exit(-1);
    }

  for (size_t i = 0; i < bin_bounds.size(); ++i)
    for (size_t j = 0; j < bin_bounds[i].size(); ++j)
      wigout << bin_bounds[i][j] << "\t" << scores[i][j] << endl;

  if (!wigout.good())
    {
      cerr << "write_bed_file: error when writing " << wigfile_name << endl;
      exit(-1);
    }

  wigout.close();
}

void
write_wigfile(const vector<double> &fg_scores,
              const vector<double> &bg_scores,
              const vector<vector<SimpleGenomicRegion> > &bin_bounds,
              const string &wigfile_name)
{
  std::ofstream wigout(wigfile_name.c_str());

  if (!wigout.good())
    {
      cerr << "write_bed_file: cannot open " << wigfile_name << endl;
      exit(-1);
    }

  size_t k = 0;
  for (size_t i = 0; i < bin_bounds.size(); ++i)
    for (size_t j = 0; j < bin_bounds[i].size(); ++j){
      wigout << bin_bounds[i][j] << "\t" << fg_scores[k] << "\t"
             << bg_scores[k] << endl;
      ++k;
    }

  if (!wigout.good())
    {
      cerr << "write_bed_file: error when writing " << wigfile_name << endl;
      exit(-1);
    }

  wigout.close();
}

void
write_bed_file(const vector<vector<GenomicRegion> > &regions,
	       const string &bed_file)
{
  std::ofstream out(bed_file.c_str());
  if (!out.good())
    {
      cerr << "write_bed_file: cannot open " << bed_file << endl;
      exit(-1);
    }

  for (vector<vector<GenomicRegion> >::const_iterator i =
	 regions.begin(); i != regions.end(); ++i)
    copy(i->begin(), i->end(), std::ostream_iterator<GenomicRegion>(out, "\n"));

  if (!out.good())
    {
      cerr << "write_bed_file: error when writing " << bed_file << endl;
      exit(-1);
    }
  out.close();
}


// for two-state segmentation
void
build_domains(const vector<vector<SimpleGenomicRegion> > &bins,
              const vector<vector<bool> > &classes,
              const vector<vector<double> > &scores, // posterior score of classes[i]
              const double score_cutoff,
              vector<vector<GenomicRegion> > &domains,
              const size_t undef_domain_cutoff = std::numeric_limits<size_t>::max())
{
  static const size_t FG_LABEL = 1;
  static const size_t BG_LABEL = 0;
  static const size_t UN_LABEL = 2;

  vector<string> LABEL_NAMES(3);
  LABEL_NAMES[FG_LABEL] = "ENRICHED";
  LABEL_NAMES[BG_LABEL] = "BACKGROUND";
  LABEL_NAMES[UN_LABEL] = "UNCONFIDENT";

  vector<vector<size_t> > labels(classes.size());
  vector<vector<double> > local_scores(scores);

  // STEP I: Relabeling bins
  for (size_t i = 0; i < classes.size(); ++i)
    {
      labels[i].resize(classes[i].size(), UN_LABEL);

      for (size_t j = 0; j < classes[i].size(); ++j)
	if (scores[i][j] >= score_cutoff)
	  labels[i][j] = static_cast<size_t>(classes[i][j]);
    }

  for (size_t i = 0; i < bins.size(); ++i)
    {
      const int lim = static_cast<int>(bins[i].size());

      for (int j = 1; j < lim; ++j)
	if (labels[i][j] == UN_LABEL &&
	    labels[i][j-1] != UN_LABEL &&
	    static_cast<size_t>(classes[i][j]) ==  labels[i][j-1])
	  labels[i][j] = labels[i][j-1];

      for (int j = lim - 2; j >= 0; --j)
	if (labels[i][j] == UN_LABEL &&
	    labels[i][j+1] != UN_LABEL &&
	    static_cast<size_t>(classes[i][j]) == labels[i][j+1])
	  labels[i][j] = labels[i][j+1];

      // deal with undefined regions
      int start = 0;
      while (start < lim)
        {
	  while (start < lim && labels[i][start] != UN_LABEL)
	    ++start;

	  int end = start + 1;
	  while (end < lim && labels[i][end] == UN_LABEL)
	    ++end;

	  if (start >= lim)
	    break; // no undefined bins in this big region

	  if (bins[i][end - 1].get_end() - bins[i][start].get_start()
	      > undef_domain_cutoff) // size of undefined region is big
            {
	      start = end;
	      continue;
            }

	  //// of undefined region is small
	  //  determine the likely state of these undefined bins
	  const int fg_n = std::accumulate(classes[i].begin() + start,
					   classes[i].begin() + end,
					   static_cast<int>(0));
	  size_t label = (fg_n * 2 >= end - start) ? FG_LABEL : BG_LABEL;

	  const size_t pre_label =
	    (start >  0) ? labels[i][start-1] : UN_LABEL;
	  const size_t next_label =
	    (end < lim) ? labels[i][end] : UN_LABEL;
	  if (pre_label == next_label && pre_label != UN_LABEL)
	    label = pre_label;

	  for (int j = start; j < end; ++j)
            {
	      labels[i][j] = label;
	      if (static_cast<size_t>(classes[i][j]) != label)
		local_scores[i][j] = 1 - scores[i][j];
            }

	  start = end;
        }
    }


  // STEP II: Build domains
  domains.resize(bins.size(), vector<GenomicRegion>());
  for (size_t i = 0; i < bins.size(); ++i)
    {
      domains[i].push_back(
			   GenomicRegion(bins[i].front().get_chrom(),
					 bins[i].front().get_start(),
					 bins[i].front().get_end(),
					 "", 0, '+'));

      double current_score = local_scores[i].front();
      for (size_t j = 1; j < labels[i].size(); ++j)
	if (labels[i][j] == labels[i][j - 1])
	  current_score += local_scores[i][j];
	else
	  {
	    domains[i].back().set_end(bins[i][j - 1].get_end());
	    domains[i].back().set_name(LABEL_NAMES[labels[i][j - 1]]);
	    domains[i].back().set_score(current_score);

	    domains[i].push_back(
				 GenomicRegion(bins[i][j].get_chrom(),
					       bins[i][j].get_start(),
					       bins[i][j].get_end(),
					       "", 0, '+'));
	    current_score = local_scores[i][j];
	  }
      domains[i].back().set_end(bins[i].back().get_end());
      domains[i].back().set_name(LABEL_NAMES[labels[i].back()]);
      domains[i].back().set_score(current_score);
    }
}

// for three-state segmentation
void
build_domains(const vector<vector<SimpleGenomicRegion> > &bins,
              const vector<vector<size_t> > &classes,
              const vector<vector<double> > &scores, // posterior score of classes[i]
              const double score_cutoff,
              vector<vector<GenomicRegion> > &domains,
              const size_t undef_domain_cutoff = std::numeric_limits<size_t>::max())
{
  static const size_t FG_LABEL = 0;
  static const size_t MG_LABEL = 1;
  static const size_t BG_LABEL = 2;
  static const size_t UN_LABEL = 3;

  vector<string> LABEL_NAMES(4);
  LABEL_NAMES[FG_LABEL] = "SAMPLE-I-ENRICHED";
  LABEL_NAMES[MG_LABEL] = "NO-DIFFERENCE";
  LABEL_NAMES[BG_LABEL] = "SAMPLE-II-ENRICHED";
  LABEL_NAMES[UN_LABEL] = "UNCONFIDENT";

  vector<vector<size_t> > labels(classes.size());
  vector<vector<double> > local_scores(scores);

  // STEP I: Relabeling bins
  for (size_t i = 0; i < classes.size(); ++i)
    {
      labels[i].resize(classes[i].size(), UN_LABEL);

      for (size_t j = 0; j < classes[i].size(); ++j)
	if (scores[i][j] >= score_cutoff)
	  labels[i][j] = classes[i][j];
    }

  for (size_t i = 0; i < bins.size(); ++i)
    {
      const int lim = static_cast<int>(bins[i].size());

      for (int j = 1; j < lim; ++j)
	if (labels[i][j] == UN_LABEL &&
	    labels[i][j-1] != UN_LABEL &&
	    classes[i][j] ==  labels[i][j-1])
	  labels[i][j] = classes[i][j];

      for (int j = lim - 2; j >= 0; --j)
	if (labels[i][j] == UN_LABEL &&
	    labels[i][j+1] != UN_LABEL &&
	    classes[i][j] == labels[i][j+1])
	  labels[i][j] = classes[i][j];

      // deal with undefined regions
      int start = 0;
      while (start < lim)
        {
	  while (start < lim && labels[i][start] != UN_LABEL)
	    ++start;

	  int end = start + 1;
	  while (end < lim && labels[i][end] == UN_LABEL)
	    ++end;

	  if (start >= lim)
	    break; // no undefined bins in this big region

	  if (bins[i][end - 1].get_end() - bins[i][start].get_start()
	      > undef_domain_cutoff) // size of undefined region is big
            {
	      start = end;
	      continue;
            }

	  //// of undefined region is small
	  //  determine the likely state of these undefined bins
	  vector<size_t> label_bin_nums(3, 0);
	  for (int j = start; j < end; ++j)
	    ++label_bin_nums[classes[i][j]];

	  size_t label =
	    std::max_element(label_bin_nums.begin(),
			     label_bin_nums.end())
	    - label_bin_nums.begin();

	  const size_t pre_label =
	    (start >  0) ? labels[i][start-1] : UN_LABEL;
	  const size_t next_label =
	    (end < lim) ? labels[i][end] : UN_LABEL;

	  if (pre_label == next_label && pre_label != UN_LABEL)
	    label = pre_label;

	  for (int j = start; j < end; ++j)
            {
	      labels[i][j] = label;
	      if (classes[i][j] != label)
		local_scores[i][j] = (1 - scores[i][j]) / 2;
            }
	  start = end;
        }
    }


  // STEP II: Build domains
  domains.resize(bins.size(), vector<GenomicRegion>());
  for (size_t i = 0; i < bins.size(); ++i)
    {
      domains[i].push_back(
			   GenomicRegion(bins[i].front().get_chrom(),
					 bins[i].front().get_start(),
					 bins[i].front().get_end(),
					 "", 0, '+'));

      double current_score = local_scores[i].front();
      for (size_t j = 1; j < labels[i].size(); ++j)
	if (labels[i][j] == labels[i][j - 1])
	  current_score += local_scores[i][j];
	else
	  {
	    domains[i].back().set_end(bins[i][j - 1].get_end());
	    domains[i].back().set_name(LABEL_NAMES[labels[i][j - 1]]);
	    domains[i].back().set_score(current_score);

	    domains[i].push_back(
				 GenomicRegion(bins[i][j].get_chrom(),
					       bins[i][j].get_start(),
					       bins[i][j].get_end(),
					       "", 0, '+'));
	    current_score = local_scores[i][j];
	  }
      domains[i].back().set_end(bins[i].back().get_end());
      domains[i].back().set_name(LABEL_NAMES[labels[i].back()]);
      domains[i].back().set_score(current_score);
    }
}

struct IsEnrichedDomain : public std::unary_function<GenomicRegion, bool>
{
  bool operator()(const GenomicRegion &dom) const
  {return dom.get_name().find("ENRICHED") != string::npos;}
};

void
pick_domains(const vector<vector<SimpleGenomicRegion> > &bins,
             const vector<vector<double> > &read_counts,
             const vector<vector<double> > &scales,
             const vector<Distro> &distros,
             vector<vector<GenomicRegion> > &domains,
             const double cdf_cutoff = 0.4)
{
  static const size_t FG_LABEL = 1;
  static const size_t BG_LABEL = 0;
  static const size_t UN_LABEL = 2;

  vector<string> LABEL_NAMES(3);
  LABEL_NAMES[FG_LABEL] = "ENRICHED";
  LABEL_NAMES[BG_LABEL] = "BACKGROUND";
  LABEL_NAMES[UN_LABEL] = "UNCONFIDENT";

  vector<vector<double> > domain_means;

  double max_count(0);

  domain_means.resize(domains.size());
  for (size_t i = 0; i < read_counts.size(); ++i)
    {
      domain_means[i].resize(domains[i].size());

      size_t k = 0;
      double domain_read_count = 0;
      double domain_size = 0;
      for (size_t j = 0; j < read_counts[i].size(); ++j)
        {
	  if (!domains[i][k].overlaps(bins[i][j]))
            {
	      domain_means[i][k] = domain_read_count / domain_size;
	      domain_read_count = 0;
	      domain_size = 0;
	      ++k;
            }
	  domain_read_count += read_counts[i][j];
	  domain_size += scales[i][j];

	  if (read_counts[i][j] > max_count)
	    max_count = read_counts[i][j];
        }
      domain_means[i][k] = domain_read_count / domain_size;
    }

  vector<double> fg_cdfs(static_cast<size_t>(max_count + 1)),
    bg_cdfs(static_cast<size_t>(max_count + 1));

  fg_cdfs[0] = exp(distros.front().log_likelihood(0));
  bg_cdfs[0] = exp(distros.back().log_likelihood(0));

  for (size_t i = 1; i < fg_cdfs.size(); ++i)
    {
      fg_cdfs[i] = fg_cdfs[i - 1] + exp(distros.front().log_likelihood(i));
      bg_cdfs[i] = bg_cdfs[i - 1] + exp(distros.back().log_likelihood(i));
    }

  for (size_t i = 0; i < domains.size(); ++i)
    for (size_t j = 0; j < domains[i].size(); ++j)
      {
	string name = domains[i][j].get_name();
	const size_t c = static_cast<size_t>(floor(domain_means[i][j]));
	if (name == LABEL_NAMES[FG_LABEL] && fg_cdfs[c] < cdf_cutoff)
	  name = LABEL_NAMES[UN_LABEL];
	if (name == LABEL_NAMES[BG_LABEL] && 1 - bg_cdfs[c] < cdf_cutoff)
	  name = LABEL_NAMES[UN_LABEL];
	domains[i][j].set_name(name + "\t"
			       + smithlab::toa(domain_means[i][j]));
      }
  for (size_t i = 0; i < domains.size(); ++i)
    domains[i].erase(std::stable_partition(domains[i].begin(),
                                           domains[i].end(),
                                           IsEnrichedDomain()),
                     domains[i].end());
}

void
pick_domains(const vector<vector<SimpleGenomicRegion> > &bins,
             const vector<vector<double> > &read_counts,
             const vector<vector<double> > &scales,
             const vector<SplitDistro> &distros,
             vector<vector<GenomicRegion> > &domains,
             const double cdf_cutoff = 0.2)
{
  static const size_t FG_LABEL = 1;
  static const size_t BG_LABEL = 0;
  static const size_t UN_LABEL = 2;

  vector<string> LABEL_NAMES(3);
  LABEL_NAMES[FG_LABEL] = "ENRICHED";
  LABEL_NAMES[BG_LABEL] = "BACKGROUND";
  LABEL_NAMES[UN_LABEL] = "UNCONFIDENT";

  vector<vector<double> > domain_means;

  double max_count(0), min_count(std::numeric_limits<double>::max());

  domain_means.resize(domains.size());
  for (size_t i = 0; i < read_counts.size(); ++i)
    {
      domain_means[i].resize(domains[i].size());

      size_t k = 0;
      double domain_read_count = 0;
      double domain_size = 0;
      for (size_t j = 0; j < read_counts[i].size(); ++j)
        {
	  if (!domains[i][k].overlaps(bins[i][j]))
            {
	      domain_means[i][k] = domain_read_count / domain_size;
	      domain_read_count = 0;
	      domain_size = 0;
	      ++k;
            }
	  domain_read_count += read_counts[i][j];
	  domain_size += scales[i][j];

	  if (read_counts[i][j] > max_count)
	    max_count = read_counts[i][j];

	  if (read_counts[i][j] < min_count)
	    min_count = read_counts[i][j];
        }
      domain_means[i][k] = domain_read_count / domain_size;
    }

  const double offset = min_count;

  vector<double> fg_cdfs(static_cast<size_t>(max_count - min_count + 1)),
    bg_cdfs(static_cast<size_t>(max_count - min_count + 1));
  fg_cdfs[0] = exp(distros.front().log_likelihood(0 + offset));
  bg_cdfs[0] = exp(distros.back().log_likelihood(0 + offset));
  for (size_t i = 1; i < fg_cdfs.size(); ++i)
    {
      fg_cdfs[i] = fg_cdfs[i - 1] + exp(distros.front().log_likelihood(i + offset));
      bg_cdfs[i] = bg_cdfs[i - 1] + exp(distros.back().log_likelihood(i + offset));
    }

  for (size_t i = 0; i < domains.size(); ++i)
    for (size_t j = 0; j < domains[i].size(); ++j)
      {
	string name = domains[i][j].get_name();
	const size_t c = static_cast<size_t>(floor(domain_means[i][j]));

	if (name == LABEL_NAMES[FG_LABEL] &&
	    fg_cdfs[static_cast<size_t>(c - offset)] < cdf_cutoff)
	  name = LABEL_NAMES[UN_LABEL];

	if (name == LABEL_NAMES[BG_LABEL] &&
	    1 - bg_cdfs[static_cast<size_t>(c - offset)] < cdf_cutoff)
	  name = LABEL_NAMES[UN_LABEL];

	domains[i][j].set_name(name + "\t"
			       + smithlab::toa(domain_means[i][j]));
      }
  for (size_t i = 0; i < domains.size(); ++i)
    domains[i].erase(std::stable_partition(domains[i].begin(),
                                           domains[i].end(),
                                           IsEnrichedDomain()),
                     domains[i].end());
}

void
pick_domains_3s(const vector<vector<SimpleGenomicRegion> > &bins,
                const vector<vector<double> > &read_counts,
                const vector<vector<double> > &scales,
                const vector<SplitDistro> &distros,
                vector<vector<GenomicRegion> > &domains,
                const double cdf_cutoff = 0.4)
{
  static const size_t FG_LABEL = 0;
  static const size_t MG_LABEL = 1;
  static const size_t BG_LABEL = 2;
  static const size_t UN_LABEL = 3;

  vector<string> LABEL_NAMES(4);
  LABEL_NAMES[FG_LABEL] = "SAMPLE-I-ENRICHED";
  LABEL_NAMES[MG_LABEL] = "NO-DIFFERENCE";
  LABEL_NAMES[BG_LABEL] = "SAMPLE-II-ENRICHED";
  LABEL_NAMES[UN_LABEL] = "UNCONFIDENT";

  vector<vector<double> > domain_means;

  double max_count(0), min_count(std::numeric_limits<double>::max());

  domain_means.resize(domains.size());
  for (size_t i = 0; i < read_counts.size(); ++i)
    {
      domain_means[i].resize(domains[i].size());

      size_t k = 0;
      double domain_read_count = 0;
      double domain_size = 0;
      for (size_t j = 0; j < read_counts[i].size(); ++j)
        {
	  if (!domains[i][k].overlaps(bins[i][j]))
            {
	      domain_means[i][k] = domain_read_count / domain_size;
	      domain_read_count = 0;
	      domain_size = 0;
	      ++k;
            }
	  domain_read_count += read_counts[i][j];
	  domain_size += scales[i][j];

	  if (read_counts[i][j] > max_count)
	    max_count = read_counts[i][j];

	  if (read_counts[i][j] < min_count)
	    min_count = read_counts[i][j];
        }
      domain_means[i][k] = domain_read_count / domain_size;
    }

  const double offset = min_count;

  vector<double> fg_cdfs(static_cast<size_t>(max_count - min_count + 1)),
    bg_cdfs(static_cast<size_t>(max_count - min_count + 1));
  fg_cdfs[0] = exp(distros.front().log_likelihood(0 + offset));
  bg_cdfs[0] = exp(distros.back().log_likelihood(0 + offset));
  for (size_t i = 1; i < fg_cdfs.size(); ++i)
    {
      fg_cdfs[i] = fg_cdfs[i - 1] + exp(distros.front().log_likelihood(i + offset));
      bg_cdfs[i] = bg_cdfs[i - 1] + exp(distros.back().log_likelihood(i + offset));
    }

  for (size_t i = 0; i < domains.size(); ++i)
    for (size_t j = 0; j < domains[i].size(); ++j)
      {
	string name = domains[i][j].get_name();
	const double cf = floor(domain_means[i][j]);
	const double cc = ceil(domain_means[i][j]);

	if (name == LABEL_NAMES[FG_LABEL] &&
	    cc < distros.back().get_mean())
	  name = LABEL_NAMES[BG_LABEL];

	if (name == LABEL_NAMES[BG_LABEL] &&
	    cf > distros.front().get_mean())
	  name = LABEL_NAMES[FG_LABEL];

	if (name == LABEL_NAMES[FG_LABEL] &&
	    fg_cdfs[static_cast<size_t>(cf - offset)] < cdf_cutoff)
	  name = LABEL_NAMES[UN_LABEL];

	if (name == LABEL_NAMES[BG_LABEL] &&
	    1 - bg_cdfs[static_cast<size_t>(cc - offset)] < cdf_cutoff)
	  name = LABEL_NAMES[UN_LABEL];

	domains[i][j].set_name(name + "\t"
			       + smithlab::toa(domain_means[i][j]));
      }
  for (size_t i = 0; i < domains.size(); ++i)
    domains[i].erase(std::stable_partition(domains[i].begin(),
                                           domains[i].end(),
                                           IsEnrichedDomain()),
                     domains[i].end());
}

