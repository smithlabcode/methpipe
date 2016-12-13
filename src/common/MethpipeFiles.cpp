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

#include <cmath>
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
using std::round;
using std::istream;
using std::ostream;
using std::endl;

string
methpipe::skip_header(std::istream &in){
  string line;
  getline(in, line);
  while((line.substr(0,1)) == "#"){
    getline(in, line);
  }
  return line;
}

void
methpipe::load_cpgs(const string &cpgs_file,
                    vector<SimpleGenomicRegion> &cpgs,
                    vector<pair<double, double> > &meths,
                    vector<size_t> &reads) {
  string chrom, prev_chrom;
  size_t pos, prev_pos = 0;
  string strand, seq;
  double meth;
  size_t coverage;

  std::ifstream in(cpgs_file.c_str());
  string line = skip_header(in); // added
  std::istringstream iss(line); //added

  iss >> chrom >> pos >> strand >> seq >> meth >> coverage;
  prev_chrom = chrom;
  prev_pos = pos;

  // append site
  cpgs.push_back(SimpleGenomicRegion(chrom, pos, pos+1));
  reads.push_back(coverage);
  meths.push_back(std::make_pair(0.0, 0.0));
  meths.back().first = static_cast<size_t>(round(meth * coverage));
  meths.back().second = static_cast<size_t>(coverage  - meths.back().first);


  while (in >> chrom >> pos >> strand >> seq >> meth >> coverage) { //changed
    // sanity check
    if (chrom.empty() || strand.empty() || seq.empty() ||
        meth < 0.0 || meth > 1.0) {
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
    meths.back().first = static_cast<size_t>(round(meth * coverage));
    meths.back().second = static_cast<size_t>(coverage  - meths.back().first);
  }
}

void
methpipe::load_cpgs(const string &cpgs_file,
                    vector<GenomicRegion> &cpgs,
                    vector<pair<double, double> > &meths,
                    vector<size_t> &reads) {
  string chrom, prev_chrom;
  size_t pos, prev_pos = 0;
  string strand, seq;
  double meth;
  size_t coverage;

  std::ifstream in(cpgs_file.c_str());
  string line = skip_header(in); //added
  std::istringstream iss(line); //added

  iss >> chrom >> pos >> strand >> seq >> meth >> coverage;
  prev_chrom = chrom;
  prev_pos = pos;

  // append site
  cpgs.push_back(GenomicRegion(chrom, pos, pos+1, seq, 0, strand[0]));
  reads.push_back(coverage);
  meths.push_back(std::make_pair(0.0, 0.0));
  meths.back().first = static_cast<size_t>(round(meth * coverage));
  meths.back().second = static_cast<size_t>(coverage  - meths.back().first);

  while (in >> chrom >> pos >> strand >> seq >> meth >> coverage) {//changed
    // sanity check
    if (chrom.empty() || strand.empty() || seq.empty()
        || meth < 0.0 || meth > 1.0) {
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
    meths.back().first = static_cast<size_t>(round(meth * coverage));
    meths.back().second = static_cast<size_t>(coverage  - meths.back().first);
  }
}

void
methpipe::load_cpgs_old(const string &cpgs_file,
                        vector<SimpleGenomicRegion> &cpgs,
                        vector<pair<double, double> > &meths,
                        vector<size_t> &reads)
{
  std::ifstream in(cpgs_file.c_str());
  string line = skip_header(in); // added
  std::istringstream iss(line); //added
  GenomicRegion r, prev_site;

  iss >> r;

  // append site
  const string site_name(r.get_name());
  const size_t coverage = atoi(
                               site_name.substr(site_name.find_first_of(":") + 1).c_str());
  const double meth = r.get_score();
  cpgs.push_back(SimpleGenomicRegion(r));
  reads.push_back(coverage);
  meths.push_back(std::make_pair(meth, 0.0));
  meths.back().first = static_cast<size_t>(round(meth * coverage));
  meths.back().second = static_cast<size_t>(coverage  - meths.back().first);

  prev_site = r;


  while (in >> r) {
    // order check
    if (r < prev_site) {
      throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");
    }
    // append site
    const string site_name(r.get_name());
    const size_t coverage = atoi(
                                 site_name.substr(site_name.find_first_of(":") + 1).c_str());
    const double meth = r.get_score();
    cpgs.push_back(SimpleGenomicRegion(r));
    reads.push_back(coverage);
    meths.push_back(std::make_pair(meth, 0.0));
    meths.back().first = static_cast<size_t>(round(meth * coverage));
    meths.back().second = static_cast<size_t>(coverage  - meths.back().first);

    prev_site = r;
  }
}

void
methpipe::load_cpgs_old(const string &cpgs_file,
                        vector<GenomicRegion> &cpgs,
                        vector<pair<double, double> > &meths,
                        vector<size_t> &reads) {
  std::ifstream in(cpgs_file.c_str());
  string line = skip_header(in); // added
  std::istringstream iss(line); //added
  GenomicRegion r, prev_site;

  iss >> r;

  // append site
  const string site_name(r.get_name());
  const size_t coverage = atoi(
                               site_name.substr(site_name.find_first_of(":") + 1).c_str());
  const double meth = r.get_score();
  cpgs.push_back(r);
  reads.push_back(coverage);
  meths.push_back(std::make_pair(meth, 0.0));
  meths.back().first = static_cast<size_t>(round(meth * coverage));
  meths.back().second = static_cast<size_t>(coverage  - meths.back().first);

  prev_site = r;


  while (in >> r) {
    // order check
    if (r < prev_site) {
      throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");
    }
    // append site
    const string site_name(r.get_name());
    const size_t coverage = atoi(
                                 site_name.substr(site_name.find_first_of(":") + 1).c_str());
    const double meth = r.get_score();
    cpgs.push_back(r);
    reads.push_back(coverage);
    meths.push_back(std::make_pair(meth, 0.0));
    meths.back().first = static_cast<size_t>(round(meth * coverage));
    meths.back().second = static_cast<size_t>(coverage  - meths.back().first);

    prev_site = r;
  }
}


bool
methpipe::is_methpipe_file_single(const string &file) {

  std::ifstream in(file.c_str());
  if (!in)
    throw SMITHLABException("could not open file: " + file);

  string line;
  if (!std::getline(in, line))
    throw SMITHLABException("could not read file: " + file);
  line = skip_header(in);
  std::istringstream iss(line);

  string chrom, strand, name;
  size_t pos = 0, coverage = 0;
  double meth = 0.0;
  iss >> chrom >> pos >> strand >> name >> meth >> coverage;

  if (strand != "+" && strand != "-") return false;

  if (meth < 0.0 || meth > 1.0) return false;

  if (name.find_first_not_of("ACGTacgtpHXx") != string::npos)
    return false;

  return true;
}


static void
move_to_start_of_line(std::istream &in) {
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
                    const size_t idx) {
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

  while (step_size > 0) {
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


std::istream&
methpipe::read_site(std::istream &in, string &chrom, size_t &pos,
                    string &strand, string &seq,
                    double &meth, size_t &coverage) {
  string line = skip_header(in);
  std::istringstream iss(line);
  iss >> chrom >> pos >> strand >> seq >> meth >> coverage;
  return in;
}


std::istream&
methpipe::read_site(std::istream &in, string &chrom, size_t &pos,
                    char &strand, string &seq,
                    double &meth, size_t &coverage) {
  string line = skip_header(in);
  std::istringstream iss(line);
  iss >> chrom >> pos >> strand >> seq >> meth >> coverage;
  return in;
}


ostream &
methpipe::write_site(ostream &out,
                     const string &chrom, const size_t &pos,
                     const string &strand, const string &seq,
                     const double &meth, const size_t &coverage) {
  return (out << chrom << "\t" << pos << "\t" << strand
          << "\t" << seq << "\t" << (coverage == 0 ? 0.0 : meth) << "\t"
          << coverage << '\n');
}


istream&
methpipe::read_site_old(istream &in, string &chrom, size_t &pos,
                        string &strand, std::string &seq,
                        double &meth, size_t &coverage)
{
  string line = skip_header(in); //added
  std::istringstream iss(line);//added
  size_t dummy;
  string name;
  if (iss >> chrom >> pos >> dummy >> name >> meth >> strand) {
    seq = name.substr(0, name.find(":"));
    coverage = atoi(name.substr(name.find(":") + 1).c_str());
  }
  else if (!line.empty())
    throw SMITHLABException("Invalid input line:" + line);
  return in;
}


ostream &
methpipe::write_site_old(ostream &out, const string &chrom, const size_t &pos,
                         const string &strand, const string &seq,
                         const double &meth, const size_t &coverage) {
  return (out << chrom << "\t" << pos << "\t" << pos + 1 << "\t"
          << (seq + ":" + smithlab::toa(coverage)) << "\t"
          << (coverage == 0 ? 0.0 : meth) << "\t" << strand << endl);
}


// IO functions for methdiff output
ostream &
methpipe::write_methdiff_site(ostream &out, const std::string &chrom,
                              const size_t pos, const std::string &strand,
                              const std::string &seq, const double diffscore,
                              const size_t meth_a, const size_t unmeth_a,
                              const size_t meth_b, const size_t unmeth_b) {
  return (out << chrom << "\t" << pos << "\t" << strand << "\t" << seq << "\t"
          << diffscore << "\t" << meth_a << "\t" << unmeth_a << "\t"
          << meth_b << "\t" << unmeth_b << endl);
}


istream &
methpipe::read_methdiff_site(istream &in, std::string &chrom,
                             size_t &pos, std::string &strand,
                             std::string &seq, double &diffscore,
                             size_t &meth_a, size_t &unmeth_a,
                             size_t &meth_b, size_t &unmeth_b) {
  return (in >> chrom >> pos >> strand >> seq >>
          diffscore >> meth_a >> unmeth_a >> meth_b >> unmeth_b);
}
