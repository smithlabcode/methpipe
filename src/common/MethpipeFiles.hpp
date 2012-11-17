/*
  Copyright (C) 2012 University of Southern California
  Authors: Andrew D. Smith, Song Qiang

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with This program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef METHPIPE_FILES_HPP
#define METHPIPE_FILES_HPP

#include <vector>
#include <string>
#include <utility>
#include "GenomicRegion.hpp"

namespace methpipe
{
    enum FILETYPE {OLD, NEW};

    void
    load_cpgs(const std::string &cpgs_file,
              std::vector<SimpleGenomicRegion> &cpgs,
              std::vector<std::pair<double, double> > &meths,
              std::vector<size_t> &reads);

    void
    load_cpgs(const std::string &cpgs_file,
              std::vector<GenomicRegion> &cpgs,
              std::vector<std::pair<double, double> > &meths,
              std::vector<size_t> &reads);

    bool
    is_methpipe_file_single(const std::string &file);
}
#endif

