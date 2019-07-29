/*
  Copyright (C) 2012 University of Southern California
  Authors: Andrew D. Smith, Song Qiang, Benjamin Decato

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

namespace methpipe {
  enum FILETYPE {OLD, NEW};

  std::string
  skip_header(std::istream &in);

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

  void
  load_cpgs_old(const std::string &cpgs_file,
                std::vector<SimpleGenomicRegion> &cpgs,
                std::vector<std::pair<double, double> > &meths,
                std::vector<size_t> &reads);

  void
  load_cpgs_old(const std::string &cpgs_file,
                std::vector<GenomicRegion> &cpgs,
                std::vector<std::pair<double, double> > &meths,
                std::vector<size_t> &reads);

  std::istream&
  read_site(std::istream &in, std::string &chrom, size_t &pos,
            std::string &strand, std::string &seq,
            double &meth, size_t &coverage);

  std::istream&
  read_site(std::istream &in, std::string &chrom, size_t &pos,
            char &strand, std::string &seq,
            double &meth, size_t &coverage);

  std::istream&
  read_site(std::istream &in, std::string &chrom, size_t &pos,
            std::string &strand, std::string &seq,
            double &meth, size_t &coverage, bool &is_array_data);

  std::ostream&
  write_site(std::ostream &out, const std::string &chrom, const size_t &pos,
             const std::string &strand, const std::string &seq,
             const double &meth, const size_t &coverage);

  // re-locate the file handler point to the first line
  // that are at or behind location chrom, pos
  void
  seek_site(std::istream &in, const std::string &chrom,
            const size_t pos);

  bool
  is_methpipe_file_single(const std::string &file);

  bool
  is_methpipe_file_array(const std::string &file);

  // files to support old format
  std::istream&
  read_site_old(std::istream &in, std::string &chrom, size_t &pos,
                std::string &strand, std::string &seq,
                double &meth, size_t &coverage);

  std::ostream &
  write_site_old(std::ostream &out, const std::string &chrom,
                 const size_t &pos, const std::string &strand,
                 const std::string &seq, const double &meth,
                 const size_t &coverage);

  //    functions for methdiff results I/O
  std::ostream &
  write_methdiff_site(std::ostream &out, const std::string &chrom,
                      const size_t pos, const std::string &strand,
                      const std::string &seq, const double diffscore,
                      const size_t meth_a, const size_t unmeth_a,
                      const size_t meth_b, const size_t unmeth_b);

  std::istream &
  read_methdiff_site(std::istream &in, std::string &chrom,
                     size_t &pos, std::string &strand,
                     std::string &seq, double &diffscore,
                     size_t &meth_a, size_t &unmeth_a,
                     size_t &meth_b, size_t &unmeth_b);
}
#endif
