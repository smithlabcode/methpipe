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

#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "GenomicRegion.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeFiles.hpp"

using std::vector;
using std::string;
using std::pair;
using std::ios_base;

void
methpipe::load_cpgs(const string &cpgs_file,
                    vector<SimpleGenomicRegion> &cpgs,
                    vector<pair<double, double> > &meths,
                    vector<size_t> &reads)
{
  string chrom, prev_chrom;
  size_t pos, prev_pos = 0;
  string strand, seq;
  double meth;
  size_t coverage;

  std::ifstream in(cpgs_file.c_str());
  while (in >> chrom >> pos >> strand >> seq >> meth >> coverage) {
    // sanity check
    if (chrom.empty() || strand.empty() || seq.empty()
        || meth < 0 || meth > 1 || coverage < 0) {
      std::ostringstream oss;
      oss << chrom << "\t" << pos << "\t" << strand << "\t"
          << seq << "\t" << meth << "\t" << coverage << "\n";
      throw SMITHLABException("Invalid input line:" + oss.str());
    }
    
    // order check
    if (prev_chrom > chrom || (prev_chrom == chrom && prev_pos > pos)) {
      throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");
    }
    prev_chrom = chrom;
    prev_pos = pos;
    
    // append site
    cpgs.push_back(SimpleGenomicRegion(chrom, pos, pos+1));
    reads.push_back(coverage);
    meths.push_back(std::make_pair(0.0, 0.0));
    // plus 0.5 to make sure the value is rounded correctly
    meths.back().first = static_cast<size_t>(meth * coverage + 0.5);
    meths.back().second = static_cast<size_t>(coverage  - meths.back().first);
  }
}

void
methpipe::load_cpgs(const string &cpgs_file,
                    vector<GenomicRegion> &cpgs,
                    vector<pair<double, double> > &meths,
                    vector<size_t> &reads)
{
  string chrom, prev_chrom;
  size_t pos, prev_pos = 0;
  string strand, seq;
  double meth;
  size_t coverage;

  std::ifstream in(cpgs_file.c_str());
  while (in >> chrom >> pos >> strand >> seq >> meth >> coverage) {
    // sanity check
    if (chrom.empty() || strand.empty() || seq.empty()
        || meth < 0 || meth > 1 || coverage < 0) {
      std::ostringstream oss;
      oss << chrom << "\t" << pos << "\t" << strand << "\t"
          << seq << "\t" << meth << "\t" << coverage << "\n";
      throw SMITHLABException("Invalid input line:" + oss.str());
    }
    
    // order check
    if (prev_chrom > chrom || (prev_chrom == chrom && prev_pos > pos)) {
      throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");
    }
    prev_chrom = chrom;
    prev_pos = pos;
    
    // append site
    cpgs.push_back(GenomicRegion(chrom, pos, pos+1, seq, 0, strand[0]));
    reads.push_back(coverage);
    meths.push_back(std::make_pair(0.0, 0.0));
    // plus 0.5 to make sure the value is rounded correctly
    meths.back().first = static_cast<size_t>(meth * coverage + 0.5);
    meths.back().second = static_cast<size_t>(coverage  - meths.back().first);
  }
}

bool
methpipe::is_methpipe_file_single(const string &file)
{

  std::ifstream in(file.c_str());
  string line;
  if (std::getline(in, line)) {
    in.close();
  
    vector<string> fields;
    smithlab::split_whitespace(line, fields);
    
    if (fields.size() == 6
        && fields[1].find_first_not_of("0123456789") == string::npos
        && fields[2].find_first_not_of("+-") == string::npos
        && fields[3].find_first_not_of("ACGTacgtpH") == string::npos
        && fields[4].find_first_not_of("0123456789.") == string::npos
        && atof(fields[4].c_str()) >= 0.0 &&  atof(fields[4].c_str()) <= 1.0
        && fields[5].find_first_not_of("0123456789") == string::npos)
      return true;
  }
  return false;
}


static void
move_to_start_of_line(std::istream &in)
{
    char next;
    while (in.good() && in.get(next) && next != '\n')
    {
        in.unget();
        in.unget();
    }
    if (in.bad()) 
        // hope this only happens when hitting the start of the file
        in.clear();
}

void
methpipe::seek_site(std::istream &in, const std::string &chr,
                    const size_t idx) 
{
    in.seekg(0, ios_base::beg);
    const size_t begin_pos = in.tellg();
    in.seekg(0, ios_base::end);
    size_t step_size = (static_cast<size_t>(in.tellg()) - begin_pos)/2;
  
    in.seekg(0, ios_base::beg);
    string low_chr;
    size_t low_idx = 0;
    in >> low_chr >> low_idx;

    // MAGIC: need the -2 here to get past the EOF and possibly a '\n'
    in.seekg(-2, ios_base::end);
    move_to_start_of_line(in);
    string high_chr;
    size_t high_idx;
    in >> high_chr >> high_idx;
  
    size_t pos = step_size;
    in.seekg(pos, ios_base::beg);
    move_to_start_of_line(in);
  
    while (step_size > 0)
    {
        string mid_chr;
        size_t mid_idx = 0;
        in >> mid_chr >> mid_idx;
        step_size /= 2;
        if (chr < mid_chr || (chr == mid_chr && idx <= mid_idx)) {
            std::swap(mid_chr, high_chr);
            std::swap(mid_idx, high_idx);
            pos -= step_size;
        }
        else {
            std::swap(mid_chr, low_chr);
            std::swap(mid_idx, low_idx);
            pos += step_size;
        }
        in.seekg(pos, ios_base::beg);
        move_to_start_of_line(in);
    }
}

bool
methpipe::read_site(std::istream &in, string &chrom, size_t &pos,
                    string &strand, string &seq,
                    double &meth, size_t &coverage) {
  string line;
  getline(in, line);

  std::istringstream is(line);

  string pos_str, meth_str, cov_str;
  if (!(is >> chrom >> pos_str >> strand >>
        seq >> meth_str >> cov_str)) {
    return false;
  }

  is.clear();
  is.str(pos_str);
  if (!(is >> pos))
    return false;

  is.clear();
  is.str(meth_str);
  meth= strtod(meth_str.c_str(), NULL);

  is.clear();
  is.str(cov_str);
  if (!(is >> coverage))
    return false;

  if (isnan(meth) && coverage == 0)
    meth = 0;
  else
    return false;

  return in.good();
}

bool
methpipe::write_site(std::ostream &out,
                     const string &chrom, const size_t &pos,
                     const string &strand, const string &seq,
                     const double &meth, const size_t &coverage) {
  return (out << chrom << "\t" << pos << "\t" << strand
          << "\t" << seq << "\t" << meth << "\t"
          << coverage << '\n');
}
