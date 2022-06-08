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

#ifndef SORT_GENOMIC_REGION_HPP
#define SORT_GENOMIC_REGION_HPP

#include "GenomicRegion.hpp"

namespace SortGenomicRegion {
void
sort_regions(std::vector<GenomicRegion> &regions);

void
sort_regions_collapse_chrom(std::vector<GenomicRegion> &regions);

void
sort_regions_collapse(std::vector<GenomicRegion> &regions);

void
sort_regions(std::vector<SimpleGenomicRegion> &regions);

void
sort_regions_collapse_chrom(std::vector<SimpleGenomicRegion> &regions);

void
sort_regions_collapse(std::vector<SimpleGenomicRegion> &regions);
}
#endif
