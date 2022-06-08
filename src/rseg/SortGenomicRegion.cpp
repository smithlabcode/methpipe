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

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "GenomicRegion.hpp"
#include "SortGenomicRegion.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::sort;

typedef GenomicRegion* GenomicRegionPointer;

struct region_pointer_less {
  bool operator()(const GenomicRegionPointer a,
		  const GenomicRegionPointer b) const {
    return (*a) < (*b);
  }
};

void
SortGenomicRegion::sort_regions(vector<GenomicRegion> &regions) {
  vector<GenomicRegionPointer> sorter;
  for (vector<GenomicRegion>::iterator i = regions.begin();
       i != regions.end(); ++i) sorter.push_back(&(*i));
  sort(sorter.begin(), sorter.end(), region_pointer_less());

  vector<GenomicRegion> r;
  r.reserve(regions.size());
  for (vector<GenomicRegionPointer>::const_iterator i(sorter.begin());
       i != sorter.end(); ++i)
    r.push_back(*(*i));
  r.swap(regions);
}

void
SortGenomicRegion::sort_regions_collapse_chrom(vector<GenomicRegion> &regions) {
  static const string FAKE_NAME("X");
  const string chrom(regions.front().get_chrom());

  vector<pair<size_t, bool> > boundaries;
  for (size_t i = 0; i < regions.size(); ++i) {
    boundaries.push_back(make_pair(regions[i].get_start(), false));
    boundaries.push_back(make_pair(regions[i].get_end(), true));
  }
  regions.clear();
  sort(boundaries.begin(), boundaries.end());

  GenomicRegion holder(chrom, 0, 0, FAKE_NAME, 0, '+');
  size_t count = 0;
  for (size_t i = 0; i < boundaries.size(); ++i)
    if (boundaries[i].second) {
      --count;
      if (count == 0) {
	holder.set_end(boundaries[i].first);
	regions.push_back(holder);
      }
    }
    else {
      if (count == 0)
	holder.set_start(boundaries[i].first);
      ++count;
    }
}

void
SortGenomicRegion::sort_regions_collapse(vector<GenomicRegion> &regions) {
  vector<vector<GenomicRegion> > separated_by_chrom;
  separate_chromosomes(regions, separated_by_chrom);
  regions.clear();
  for (size_t i = 0; i < separated_by_chrom.size(); ++i) {
    sort_regions_collapse_chrom(separated_by_chrom[i]);
    regions.insert(regions.end(),
		   separated_by_chrom[i].begin(),
		   separated_by_chrom[i].end());
    separated_by_chrom[i].clear();
  }
}

typedef SimpleGenomicRegion* SimpleGenomicRegionPointer;

struct simple_region_pointer_less {
  bool operator()(const SimpleGenomicRegionPointer a,
		  const SimpleGenomicRegionPointer b) const {
    return (*a) < (*b);
  }
};

void
SortGenomicRegion::sort_regions(vector<SimpleGenomicRegion> &regions) {
  vector<SimpleGenomicRegionPointer> sorter;
  for (vector<SimpleGenomicRegion>::iterator i = regions.begin();
       i != regions.end(); ++i) sorter.push_back(&(*i));
  sort(sorter.begin(), sorter.end(), simple_region_pointer_less());

  vector<SimpleGenomicRegion> r;
  r.reserve(regions.size());
  for (vector<SimpleGenomicRegionPointer>::const_iterator i(sorter.begin());
       i != sorter.end(); ++i)
    r.push_back(*(*i));
  r.swap(regions);
}

void
SortGenomicRegion::sort_regions_collapse_chrom(vector<SimpleGenomicRegion> &regions) {
  static const string FAKE_NAME("X");
  const string chrom(regions.front().get_chrom());

  vector<pair<size_t, bool> > boundaries;
  for (size_t i = 0; i < regions.size(); ++i) {
    boundaries.push_back(make_pair(regions[i].get_start(), false));
    boundaries.push_back(make_pair(regions[i].get_end(), true));
  }
  regions.clear();
  sort(boundaries.begin(), boundaries.end());

  SimpleGenomicRegion holder(chrom, 0, 0);
  size_t count = 0;
  for (size_t i = 0; i < boundaries.size(); ++i)
    if (boundaries[i].second) {
      --count;
      if (count == 0) {
	holder.set_end(boundaries[i].first);
	regions.push_back(holder);
      }
    }
    else {
      if (count == 0)
	holder.set_start(boundaries[i].first);
      ++count;
    }
}

void
SortGenomicRegion::sort_regions_collapse(vector<SimpleGenomicRegion> &regions) {
  vector<vector<SimpleGenomicRegion> > separated_by_chrom;
  separate_chromosomes(regions, separated_by_chrom);
  regions.clear();
  for (size_t i = 0; i < separated_by_chrom.size(); ++i) {
    sort_regions_collapse_chrom(separated_by_chrom[i]);
    regions.insert(regions.end(),
		   separated_by_chrom[i].begin(),
		   separated_by_chrom[i].end());
    separated_by_chrom[i].clear();
  }
}

