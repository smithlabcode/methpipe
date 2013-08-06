/*
 *    Part of RMAP software
 *
 *    Copyright (C) 2009 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "load_paired_end_reads.hpp"
#include "SeedMaker.hpp"
#include "smithlab_os.hpp"
#include "QualityScore.hpp"

#include <cstring>
#include <fstream>

using std::vector;
using std::string;
using std::ptr_fun;
using std::not1;
using std::min;

using std::cerr;
using std::endl;

static const int INPUT_BUFFER_SIZE = 10000;

static char
to_base_symbol(char c) {return (isvalid(c)) ? toupper(c) : 'N';}

static size_t
get_read_word(const string &read) {
  const size_t trunc_to = min(read.length(), SeedMaker::max_seed_part);
  // Need to replace the Ns because otherwise they will destroy the
  // conversion from DNA to integers. Could replace with random
  // bases, but everyone hates non-deterministic programs.
  string s;
  replace_copy(read.end() - trunc_to, read.end(), std::back_inserter(s), 'N', 'A');
  return SeedMaker::make_read_word(s);
}

static void
check_and_add(string &read, const int max_diffs,
	      size_t &read_width, 
	      vector<FastRead> &fast_reads_left,
	      vector<FastRead> &fast_reads_right,
	      vector<size_t> &read_words, vector<size_t> &read_index, 
	      size_t &read_count) {
  const size_t half_width = read.length()/2;
  if (read_width == 0) {
    if (half_width*2 != read.length())
      throw SMITHLABException("PE reads not even length (must specify length)");
    read_width = read.length()/2;
  }
  else if (read.length() < 2*read_width)
    throw SMITHLABException("Incorrect read width");
  
  if (read_count == 0)
    FastRead::set_read_width(read_width);
  
  // clean the read
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));
  
  // check for quality
  const bool good_read = (count(read.begin(), read.end(), 'N') <= 2*max_diffs);
  if (good_read) {
    fast_reads_left.push_back(FastRead(read.begin(), read.begin() + read_width));
    string right_read(read.begin() + half_width,
		      read.begin() + half_width + read_width);
    revcomp_inplace(right_read);
    fast_reads_right.push_back(FastRead(right_read));
    read_words.push_back(get_read_word(right_read));
    read_index.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_fasta_file(const string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   vector<FastRead> &fast_reads_left,
			   vector<FastRead> &fast_reads_right,
			   vector<size_t> &read_words, vector<size_t> &read_index) {
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  
  char buffer[INPUT_BUFFER_SIZE + 1];
  
  size_t read_count = 0;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      if (buffer[0] != '>') {
	string read(buffer);
	check_and_add(read, max_diffs, read_width, fast_reads_left, fast_reads_right,
		      read_words, read_index, read_count);
      }
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

inline bool
is_fastq_name_line(size_t line_count) {
  return ((line_count % 4) == 0);
}

inline bool
is_fastq_sequence_line(size_t line_count) {
  return ((line_count % 4) == 1);
}

inline bool
is_fastq_score_name_line(size_t line_count) {
  return ((line_count % 4) == 2);
}

inline bool
is_fastq_score_line(size_t line_count) {
  return ((line_count % 4) == 3);
}

void
load_reads_from_fastq_file(const string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   vector<FastRead> &fast_reads_left,
			   vector<FastRead> &fast_reads_right,
			   vector<size_t> &read_words, vector<size_t> &read_index) {
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];

  size_t read_count = 0, line_count = 0;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      
      // if (is_fastq_name_line(line_count))
      //   if (buffer[0] != '@')
      //     throw SMITHLABException("invalid FASTQ name line: " + string(buffer));
      if (is_fastq_sequence_line(line_count)) {
	string read(buffer);
	check_and_add(read, max_diffs, read_width, fast_reads_left, fast_reads_right, 
		      read_words, read_index, read_count);
      }
      //  if (is_fastq_score_name_line(line_count))
      //    if (buffer[0] != '+')
      //      throw SMITHLABException("invalid FASTQ score name line: " + string(buffer));
      //  if (is_fastq_score_line(line_count))
      //    ; //!!!!!!!!!!!!!
      ++line_count;
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}

static void
matrix_revcomp(const vector<vector<double> >::iterator a, 
	       const vector<vector<double> >::iterator b) {
  for (vector<vector<double> >::iterator i(a); i != b; ++i)
    reverse(i->begin(), i->end());
  reverse(a, b);
}

static void
check_and_add(const FASTQScoreType score_format, const size_t max_diffs,
	      const string &score_line, string &read, size_t &read_width, 
	      vector<FastReadWC> &fast_reads_left, 
	      vector<FastReadWC> &fast_reads_right, 
	      vector<size_t> &read_words, 
	      vector<size_t> &read_index, size_t &read_count) {

  const size_t half_width = read.length()/2;
  if (read_width == 0) {
    if (half_width*2 != read.length())
      throw SMITHLABException("PE reads not even length (must specify length)");
    read_width = read.length()/2;
  }
  else if (read.length() < 2*read_width)
    throw SMITHLABException("Incorrect read width");
  
  if (read_count == 0)
    FastReadWC::set_read_width(read_width);
  
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));
  size_t bad_count = 0;
  vector<vector<double> > scores;
  for (size_t i = 0; i < score_line.length(); ++i) {
    const double error_prob = 
      quality_char_to_error_probability(score_format, score_line[i]);
    const double other_probs = 1.0 - error_prob/(smithlab::alphabet_size - 1);
    scores.push_back(vector<double>(smithlab::alphabet_size, other_probs));
    scores[i][base2int(read[i])] = error_prob;
    bad_count += (error_prob > 0.5);
  }

  const bool good_read = (bad_count <= 2*max_diffs);
  if (good_read) {
    fast_reads_left.push_back(FastReadWC(scores.begin(), scores.begin() + read_width));
    matrix_revcomp(scores.begin() + half_width,
		   scores.begin() + half_width + read_width);
    fast_reads_right.push_back(FastReadWC(scores.begin() + half_width,
					  scores.begin() + half_width + read_width));
    string right_read(read.begin() + half_width,
		      read.begin() + half_width + read_width);
    revcomp_inplace(right_read);
    read_words.push_back(get_read_word(right_read));
    read_index.push_back(read_count);
  }
  ++read_count;
}

void
load_reads_from_fastq_file(const string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   vector<FastReadWC> &fast_reads_left,
			   vector<FastReadWC> &fast_reads_right,
			   vector<size_t> &read_words, vector<size_t> &read_index) {

  FASTQScoreType score_format = fastq_score_type(filename);
  
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];

  size_t read_count = 0, line_count = 0;
  string sequence;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      
      //       if (is_fastq_name_line(line_count))
      // 	if (buffer[0] != '@')
      // 	  throw SMITHLABException("invalid FASTQ name line: " + string(buffer));
      if (is_fastq_sequence_line(line_count)) {
	sequence = string(buffer);
      }
      //       if (is_fastq_score_name_line(line_count))
      // 	if (buffer[0] != '+')
      // 	  throw SMITHLABException("invalid FASTQ score name line: " + string(buffer));
      if (is_fastq_score_line(line_count)) {
	const string score_line(buffer);
	check_and_add(score_format, max_diffs, score_line, sequence, read_width, 
		      fast_reads_left, fast_reads_right,
		      read_words, read_index, read_count);
      }
      ++line_count;
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}


static void
check_and_add(const FASTQScoreType score_format, const size_t max_diffs,
	      const string &score_line, string &read, size_t &read_width, 
	      vector<FastReadQuality> &fast_reads_left, 
	      vector<FastReadQuality> &fast_reads_right, 
	      vector<size_t> &read_words, 
	      vector<size_t> &read_index, size_t &read_count) {
  
  const size_t half_width = read.length()/2;
  if (read_width == 0) {
    if (half_width*2 != read.length())
      throw SMITHLABException("PE reads not even length (must specify length)");
    read_width = read.length()/2;
  }
  else if (read.length() < 2*read_width)
    throw SMITHLABException("Incorrect read width");
  if (read_count == 0)
    FastReadQuality::set_read_width(read_width);
  
  // clean the read
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));
  size_t bad_count = 0;
  vector<vector<double> > scores;
  for (size_t i = 0; i < score_line.length(); ++i) {
    const double error_prob = 
      quality_char_to_error_probability(score_format, score_line[i]);
    const double other_probs = 1.0 - error_prob/(smithlab::alphabet_size - 1);
    scores.push_back(vector<double>(smithlab::alphabet_size, other_probs));
    scores[i][base2int(read[i])] = error_prob;
    bad_count += (error_prob > 0.5);
  }
  
  const bool good_read = (bad_count <= 2*max_diffs);
  
  if (good_read) {
    fast_reads_left.push_back(FastReadQuality(scores.begin(), scores.begin() + read_width));
    matrix_revcomp(scores.begin() + half_width,
		   scores.begin() + half_width + read_width);
    fast_reads_right.push_back(FastReadQuality(scores.begin() + half_width,
					       scores.begin() + half_width + read_width));
    string right_read(read.begin() + half_width,
		      read.begin() + half_width + read_width);
    revcomp_inplace(right_read);
    read_words.push_back(get_read_word(right_read));
    read_index.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_fastq_file(const string &filename, const size_t max_diffs,
			   size_t &read_width, 
			   vector<FastReadQuality> &fast_reads_left,
			   vector<FastReadQuality> &fast_reads_right,
			   vector<size_t> &read_words, vector<size_t> &read_index) {
  FASTQScoreType score_format = fastq_score_type(filename);

  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];
  
  size_t read_count = 0, line_count = 0;
  string sequence;

  
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos/mac carriage returns before newlines
      const size_t last_pos = in.gcount() - 2; //strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      
      //       if (is_fastq_name_line(line_count))
      // 	;
      if (is_fastq_sequence_line(line_count))
	sequence = string(buffer);
      //       if (is_fastq_score_name_line(line_count))
      // 	;
      if (is_fastq_score_line(line_count)) {
	const string score_line(buffer);
	check_and_add(score_format, max_diffs, score_line, sequence, 
		      read_width, fast_reads_left, fast_reads_right,
		      read_words, read_index, read_count);
      }
      ++line_count;
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static void
check_and_add(const FASTQScoreType score_format, const size_t max_diffs,
	      const string &score_line, size_t &read_width, 
	      vector<FastReadWC> &fast_reads_left, 
	      vector<FastReadWC> &fast_reads_right, 
	      vector<size_t> &read_words, 
	      vector<size_t> &read_index, size_t &read_count) {
  
  // parse the score line
  vector<string> parts;
  smithlab::split_whitespace(score_line, parts);
  if (parts.size() % smithlab::alphabet_size != 0)
    throw SMITHLABException("bad format:\n" + score_line);

  // check the read width
  const size_t half_width = (parts.size()/2)/smithlab::alphabet_size;
  if (read_width == 0) {
    if (half_width*2 != parts.size()/smithlab::alphabet_size)
      throw SMITHLABException("PE reads not even length (must specify length)");
    read_width = half_width;
  }
  else if (parts.size()/smithlab::alphabet_size < 2*read_width)
    throw SMITHLABException("Incorrect read width");
  if (read_count == 0)
    FastReadWC::set_read_width(read_width);
  
  // convert to numerical values
  const size_t lim = parts.size()/smithlab::alphabet_size;
  vector<vector<double> > error_probs(lim, vector<double>(smithlab::alphabet_size));
  for (size_t i = 0; i < lim; ++i)
    for (size_t j = 0; j < smithlab::alphabet_size; ++j)
      error_probs[i][j] = atof(parts[i*smithlab::alphabet_size + j].c_str());
  
  // convert to probability
  size_t bad_count = 0;
  for (size_t i = 0; i < lim; ++i) {
    for (size_t j = 0; j < smithlab::alphabet_size; ++j)
      error_probs[i][j] = 
	quality_score_to_error_probability(score_format, error_probs[i][j]);
    bad_count += (*min_element(error_probs[i].begin(), error_probs[i].end()) > 0.995);
  }
  
  const bool good_read = (bad_count <= 2*max_diffs);
  if (good_read) {
    fast_reads_left.push_back(FastReadWC(error_probs.begin(), error_probs.begin() + read_width));
    matrix_revcomp(error_probs.begin() + half_width,
		   error_probs.begin() + half_width + read_width);
    fast_reads_right.push_back(FastReadWC(error_probs.begin() + half_width,
					  error_probs.begin() + half_width + read_width));
    string read;
    for (size_t i = 0; i < error_probs.size(); ++i)
      read += int2base(min_element(error_probs[i].begin(), 
				   error_probs[i].end()) - error_probs[i].begin());
    string right_read_already_rc(read.begin() + half_width,
				 read.begin() + half_width + read_width);
    //  NO NEED TO REVCOMP AGAIN!! revcomp_inplace(right_read);
    read_words.push_back(get_read_word(right_read_already_rc));
    read_index.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_prb_file(const string &filename, const size_t max_diffs,
			 size_t &read_width, 
			 vector<FastReadWC> &fast_reads_left,
			 vector<FastReadWC> &fast_reads_right,
			 vector<size_t> &read_words, vector<size_t> &read_index) {

  FASTQScoreType score_format = FASTQ_Solexa;

  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];

  size_t read_count = 0;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      const string score_line(buffer);
      check_and_add(score_format, max_diffs, score_line, read_width, 
		    fast_reads_left, fast_reads_right, 
		    read_words, read_index, read_count);
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}
 
 
static void
check_and_add(const FASTQScoreType score_format, const size_t max_diffs,
	      const string &score_line, size_t &read_width, 
	      vector<FastReadQuality> &fast_reads_left, 
	      vector<FastReadQuality> &fast_reads_right, 
	      vector<size_t> &read_words, 
	      vector<size_t> &read_index, size_t &read_count) {
  
  
  // parse the score line
  vector<string> parts;
  smithlab::split_whitespace(score_line, parts);
  if (parts.size() % smithlab::alphabet_size != 0)
    throw SMITHLABException("bad format:\n" + score_line);

  // check the read width
  const size_t half_width = (parts.size()/2)/smithlab::alphabet_size;
  if (read_width == 0) {
    if (half_width*2 != parts.size()/smithlab::alphabet_size)
      throw SMITHLABException("PE reads not even length (must specify length)");
    read_width = half_width;
  }
  else if (parts.size()/smithlab::alphabet_size < 2*read_width)
    throw SMITHLABException("Incorrect read width");
  if (read_count == 0)
    FastReadQuality::set_read_width(read_width);
  
  // convert to numerical values
  const size_t lim = parts.size()/smithlab::alphabet_size;
  vector<vector<double> > error_probs(lim, vector<double>(smithlab::alphabet_size));
  for (size_t i = 0; i < lim; ++i)
    for (size_t j = 0; j < smithlab::alphabet_size; ++j)
      error_probs[i][j] = atof(parts[i*smithlab::alphabet_size + j].c_str());
  
  // convert to probability
  size_t bad_count = 0;
  for (size_t i = 0; i < lim; ++i) {
    for (size_t j = 0; j < smithlab::alphabet_size; ++j)
      error_probs[i][j] = 
	quality_score_to_error_probability(score_format, error_probs[i][j]);
    bad_count += (*min_element(error_probs[i].begin(), error_probs[i].end()) > 0.5);
  }
  
  const bool good_read = (bad_count <= 2*max_diffs);
  if (good_read) {
    fast_reads_left.push_back(FastReadQuality(error_probs.begin(), error_probs.begin() + read_width));
    matrix_revcomp(error_probs.begin() + half_width,
		   error_probs.begin() + half_width + read_width);
    fast_reads_right.push_back(FastReadQuality(error_probs.begin() + half_width,
					       error_probs.begin() + half_width + read_width));
    string read;
    for (size_t i = 0; i < error_probs.size(); ++i)
      read += int2base(min_element(error_probs[i].begin(), 
				   error_probs[i].end()) - error_probs[i].begin());
    string right_read_already_rc(read.begin() + half_width,
				 read.begin() + half_width + read_width);
    //  NO NEED TO REVCOMP AGAIN!! revcomp_inplace(right_read);
    read_words.push_back(get_read_word(right_read_already_rc));
    read_index.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_prb_file(const string &filename, const size_t max_diffs,
			 size_t &read_width, 
			 vector<FastReadQuality> &fast_reads_left,
			 vector<FastReadQuality> &fast_reads_right,
			 vector<size_t> &read_words, vector<size_t> &read_index) {
  FASTQScoreType score_format = FASTQ_Solexa;
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];

  size_t read_count = 0;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      const string score_line(buffer);
      check_and_add(score_format, max_diffs, score_line, read_width, 
		    fast_reads_left, fast_reads_right,
		    read_words, read_index, read_count);
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static size_t
get_read_word_cheat(const string &read) {
  const size_t trunc_to = min(read.length(), SeedMaker::max_seed_part);
  // Need to replace the Ns because otherwise they will destroy the
  // conversion from DNA to integers. Could replace with random
  // bases, but everyone hates non-deterministic programs.
  string s;
  replace_copy(read.begin(), read.begin() + trunc_to, std::back_inserter(s), 'N', 'A');
  return SeedMaker::make_read_word(s);
}

static void
check_and_add_cheat(string &read, const int max_diffs,
		    size_t &read_width, 
		    vector<FastRead> &fast_reads_left,
		    vector<FastRead> &fast_reads_right,
		    vector<size_t> &read_words_l, 
		    vector<size_t> &read_words_r, 
		    vector<size_t> &read_index_l, 
		    vector<size_t> &read_index_r, 
		    size_t &read_count) {
  const size_t half_width = read.length()/2;
  if (read_width == 0) {
    if (half_width*2 != read.length())
      throw SMITHLABException("PE reads not even length (must specify length)");
    read_width = read.length()/2;
  }
  else if (read.length() < 2*read_width)
    throw SMITHLABException("Incorrect read width");
  
  if (read_count == 0)
    FastRead::set_read_width(read_width);
  
  // clean the read
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));
  
  // check for quality
  const bool good_read = (count(read.begin(), read.end(), 'N') <= 2*max_diffs);
  if (good_read) {
    const string right_read(read.begin() + half_width,
			    read.begin() + half_width + read_width);
    fast_reads_right.push_back(FastRead(right_read));
    read_words_r.push_back(get_read_word_cheat(right_read));
    read_index_r.push_back(read_count);
    const string left_read(read.begin(), read.begin() + read_width);
    fast_reads_left.push_back(FastRead(left_read));
    read_words_l.push_back(get_read_word_cheat(left_read));
    read_index_l.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_fasta_file_cheat(const string &filename, const size_t max_diffs,
				 size_t &read_width, 
				 vector<FastRead> &fast_reads_left,
				 vector<FastRead> &fast_reads_right,
				 vector<size_t> &read_words_l, 
				 vector<size_t> &read_words_r, 
				 vector<size_t> &read_index_l,
				 vector<size_t> &read_index_r) {
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  
  char buffer[INPUT_BUFFER_SIZE + 1];
  
  size_t read_count = 0;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      if (buffer[0] != '>') {
	string read(buffer);
	check_and_add_cheat(read, max_diffs, read_width, fast_reads_left, fast_reads_right,
			    read_words_l, read_words_r, read_index_l, read_index_r, read_count);
      }
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// inline bool
// is_fastq_name_line(size_t line_count) {
//   return ((line_count % 4) == 0);
// }

// inline bool
// is_fastq_sequence_line(size_t line_count) {
//   return ((line_count % 4) == 1);
// }

// inline bool
// is_fastq_score_name_line(size_t line_count) {
//   return ((line_count % 4) == 2);
// }

// inline bool
// is_fastq_score_line(size_t line_count) {
//   return ((line_count % 4) == 3);
// }

void
load_reads_from_fastq_file_cheat(const string &filename, const size_t max_diffs,
				 size_t &read_width, 
				 vector<FastRead> &fast_reads_left,
				 vector<FastRead> &fast_reads_right,
				 vector<size_t> &read_words_l, 
				 vector<size_t> &read_words_r, 
				 vector<size_t> &read_index_l,
				 vector<size_t> &read_index_r
				 ) {
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];

  size_t read_count = 0, line_count = 0;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      
      // if (is_fastq_name_line(line_count))
      //   if (buffer[0] != '@')
      //     throw SMITHLABException("invalid FASTQ name line: " + string(buffer));
      if (is_fastq_sequence_line(line_count)) {
	string read(buffer);
	check_and_add_cheat(read, max_diffs, read_width, fast_reads_left, fast_reads_right, 
			    read_words_l, read_words_r, read_index_l, read_index_r, read_count);
      }
      //  if (is_fastq_score_name_line(line_count))
      //    if (buffer[0] != '+')
      //      throw SMITHLABException("invalid FASTQ score name line: " + string(buffer));
      //  if (is_fastq_score_line(line_count))
      //    ; //!!!!!!!!!!!!!
      ++line_count;
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}

// static void
// matrix_revcomp(const vector<vector<double> >::iterator a, 
// 	       const vector<vector<double> >::iterator b) {
//   for (vector<vector<double> >::iterator i(a); i != b; ++i)
//     reverse(i->begin(), i->end());
//   reverse(a, b);
// }

static void
check_and_add_cheat(const FASTQScoreType score_format, const size_t max_diffs,
		    const string &score_line, string &read, size_t &read_width, 
		    vector<FastReadWC> &fast_reads_left, 
		    vector<FastReadWC> &fast_reads_right, 
		    vector<size_t> &read_words_l, 
		    vector<size_t> &read_words_r, 
		    vector<size_t> &read_index_l, 
		    vector<size_t> &read_index_r, 
		    size_t &read_count) {
  
  const size_t half_width = read.length()/2;
  if (read_width == 0) {
    if (half_width*2 != read.length())
      throw SMITHLABException("PE reads not even length (must specify length)");
    read_width = read.length()/2;
  }
  else if (read.length() < 2*read_width)
    throw SMITHLABException("Incorrect read width");
  
  if (read_count == 0)
    FastReadWC::set_read_width(read_width);
  
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));
  size_t bad_count = 0;
  vector<vector<double> > scores;
  for (size_t i = 0; i < score_line.length(); ++i) {
    const double error_prob = 
      quality_char_to_error_probability(score_format, score_line[i]);
    const double other_probs = 1.0 - error_prob/(smithlab::alphabet_size - 1);
    scores.push_back(vector<double>(smithlab::alphabet_size, other_probs));
    scores[i][base2int(read[i])] = error_prob;
    bad_count += (error_prob > 0.5);
  }

  const bool good_read = (bad_count <= 2*max_diffs);
  if (good_read) {
    fast_reads_right.push_back(FastReadWC(scores.begin() + half_width,
					  scores.begin() + half_width + read_width));
    const string right_read(read.begin() + half_width,
			    read.begin() + half_width + read_width);
    read_words_r.push_back(get_read_word_cheat(right_read));
    read_index_r.push_back(read_count);

    fast_reads_left.push_back(FastReadWC(scores.begin(), scores.begin() + read_width));
    const string left_read(read.begin(), read.begin() + read_width);
    read_words_l.push_back(get_read_word_cheat(left_read));
    read_index_l.push_back(read_count);
  }
  ++read_count;
}

void
load_reads_from_fastq_file_cheat(const string &filename, const size_t max_diffs,
				 size_t &read_width, 
				 vector<FastReadWC> &fast_reads_left,
				 vector<FastReadWC> &fast_reads_right,
				 vector<size_t> &read_words_l, 
				 vector<size_t> &read_words_r, 
				 vector<size_t> &read_index_l,
				 vector<size_t> &read_index_r) {
  
  FASTQScoreType score_format = fastq_score_type(filename);
  
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];

  size_t read_count = 0, line_count = 0;
  string sequence;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      
      //       if (is_fastq_name_line(line_count))
      // 	if (buffer[0] != '@')
      // 	  throw SMITHLABException("invalid FASTQ name line: " + string(buffer));
      if (is_fastq_sequence_line(line_count)) {
	sequence = string(buffer);
      }
      //       if (is_fastq_score_name_line(line_count))
      // 	if (buffer[0] != '+')
      // 	  throw SMITHLABException("invalid FASTQ score name line: " + string(buffer));
      if (is_fastq_score_line(line_count)) {
	const string score_line(buffer);
	check_and_add_cheat(score_format, max_diffs, score_line, sequence, read_width, 
			    fast_reads_left, fast_reads_right, read_words_l, read_words_r, 
			    read_index_l, read_index_r, read_count);
      }
      ++line_count;
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}


static void
check_and_add_cheat(const FASTQScoreType score_format, const size_t max_diffs,
		    const string &score_line, string &read, size_t &read_width, 
		    vector<FastReadQuality> &fast_reads_left, 
		    vector<FastReadQuality> &fast_reads_right, 
		    vector<size_t> &read_words_l, 
		    vector<size_t> &read_words_r, 
		    vector<size_t> &read_index_l, 
		    vector<size_t> &read_index_r, 
		    size_t &read_count) {
  
  const size_t half_width = read.length()/2;
  if (read_width == 0) {
    if (half_width*2 != read.length())
      throw SMITHLABException("PE reads not even length (must specify length)");
    read_width = read.length()/2;
  }
  else if (read.length() < 2*read_width)
    throw SMITHLABException("Incorrect read width");
  if (read_count == 0)
    FastReadQuality::set_read_width(read_width);
  
  // clean the read
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));
  size_t bad_count = 0;
  vector<vector<double> > scores;
  for (size_t i = 0; i < score_line.length(); ++i) {
    const double error_prob = 
      quality_char_to_error_probability(score_format, score_line[i]);
    const double other_probs = 1.0 - error_prob/(smithlab::alphabet_size - 1);
    scores.push_back(vector<double>(smithlab::alphabet_size, other_probs));
    scores[i][base2int(read[i])] = error_prob;
    bad_count += (error_prob > 0.5);
  }
  
  const bool good_read = (bad_count <= 2*max_diffs);
  
  if (good_read) {
    fast_reads_right.push_back(FastReadQuality(scores.begin() + half_width,
					       scores.begin() + half_width + read_width));
    const string right_read(read.begin() + half_width,
			    read.begin() + half_width + read_width);
    read_words_r.push_back(get_read_word_cheat(right_read));
    read_index_r.push_back(read_count);

    fast_reads_left.push_back(FastReadQuality(scores.begin(), scores.begin() + read_width));
    const string left_read(read.begin(), read.begin() + read_width);
    read_words_l.push_back(get_read_word_cheat(left_read));
    read_index_l.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_fastq_file_cheat(const string &filename, const size_t max_diffs,
				 size_t &read_width, 
				 vector<FastReadQuality> &fast_reads_left,
				 vector<FastReadQuality> &fast_reads_right,
				 vector<size_t> &read_words_l, 
				 vector<size_t> &read_words_r, 
				 vector<size_t> &read_index_l,
				 vector<size_t> &read_index_r) {
  FASTQScoreType score_format = fastq_score_type(filename);

  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];
  
  size_t read_count = 0, line_count = 0;
  string sequence;
  
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos/mac carriage returns before newlines
      const size_t last_pos = in.gcount() - 2; //strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      
      //       if (is_fastq_name_line(line_count))
      // 	;
      if (is_fastq_sequence_line(line_count))
	sequence = string(buffer);
      //       if (is_fastq_score_name_line(line_count))
      // 	;
      if (is_fastq_score_line(line_count)) {
	const string score_line(buffer);
	check_and_add_cheat(score_format, max_diffs, score_line, sequence, 
			    read_width, fast_reads_left, fast_reads_right,
			    read_words_l, read_words_r, 
			    read_index_l, read_index_r, read_count);
      }
      ++line_count;
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static void
check_and_add_cheat(const FASTQScoreType score_format, const size_t max_diffs,
		    const string &score_line, size_t &read_width, 
		    vector<FastReadWC> &fast_reads_left, 
		    vector<FastReadWC> &fast_reads_right, 
		    vector<size_t> &read_words_l, 
		    vector<size_t> &read_words_r, 
		    vector<size_t> &read_index_l, 
		    vector<size_t> &read_index_r, 
		    size_t &read_count) {
  
  // parse the score line
  vector<string> parts;
  smithlab::split_whitespace(score_line, parts);
  if (parts.size() % smithlab::alphabet_size != 0)
    throw SMITHLABException("bad format:\n" + score_line);

  // check the read width
  const size_t half_width = (parts.size()/2)/smithlab::alphabet_size;
  if (read_width == 0) {
    if (half_width*2 != parts.size()/smithlab::alphabet_size)
      throw SMITHLABException("PE reads not even length (must specify length)");
    read_width = half_width;
  }
  else if (parts.size()/smithlab::alphabet_size < 2*read_width)
    throw SMITHLABException("Incorrect read width");
  if (read_count == 0)
    FastReadWC::set_read_width(read_width);
  
  // convert to numerical values
  const size_t lim = parts.size()/smithlab::alphabet_size;
  vector<vector<double> > error_probs(lim, vector<double>(smithlab::alphabet_size));
  for (size_t i = 0; i < lim; ++i)
    for (size_t j = 0; j < smithlab::alphabet_size; ++j)
      error_probs[i][j] = atof(parts[i*smithlab::alphabet_size + j].c_str());
  
  // convert to probability
  size_t bad_count = 0;
  for (size_t i = 0; i < lim; ++i) {
    for (size_t j = 0; j < smithlab::alphabet_size; ++j)
      error_probs[i][j] = 
	quality_score_to_error_probability(score_format, error_probs[i][j]);
    bad_count += (*min_element(error_probs[i].begin(), error_probs[i].end()) > 0.995);
  }
  
  const bool good_read = (bad_count <= 2*max_diffs);
  if (good_read) {
    string read;
    for (size_t i = 0; i < error_probs.size(); ++i)
      read += int2base(min_element(error_probs[i].begin(), 
				   error_probs[i].end()) - error_probs[i].begin());
    fast_reads_right.push_back(FastReadWC(error_probs.begin() + half_width,
					  error_probs.begin() + half_width + read_width));
    const string right_read(read.begin() + half_width,
			    read.begin() + half_width + read_width);
    read_words_r.push_back(get_read_word_cheat(right_read));
    read_index_r.push_back(read_count);
    
    fast_reads_left.push_back(FastReadWC(error_probs.begin(), error_probs.begin() + read_width));
    const string left_read(read.begin(), read.begin() + read_width);
    read_words_l.push_back(get_read_word_cheat(left_read));
    read_index_l.push_back(read_count);
    
    //     fast_reads_left.push_back(FastReadWC(error_probs.begin(), error_probs.begin() + read_width));
    //     matrix_revcomp(error_probs.begin() + half_width,
    // 		   error_probs.begin() + half_width + read_width);
    //     fast_reads_right.push_back(FastReadWC(error_probs.begin() + half_width,
    // 					  error_probs.begin() + half_width + read_width));
    //     string read;
    //     for (size_t i = 0; i < error_probs.size(); ++i)
    //       read += int2base(min_element(error_probs[i].begin(), 
    // 				   error_probs[i].end()) - error_probs[i].begin());
    //     string right_read_already_rc(read.begin() + half_width,
    // 				 read.begin() + half_width + read_width);
    //     //  NO NEED TO REVCOMP AGAIN!! revcomp_inplace(right_read);
    //     read_words.push_back(get_read_word_cheat(right_read_already_rc));
    //     read_index.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_prb_file_cheat(const string &filename, const size_t max_diffs,
			       size_t &read_width, 
			       vector<FastReadWC> &fast_reads_left,
			       vector<FastReadWC> &fast_reads_right,
			       vector<size_t> &read_words_l, 
			       vector<size_t> &read_words_r, 
			       vector<size_t> &read_index_l,
			       vector<size_t> &read_index_r) {
  
  FASTQScoreType score_format = FASTQ_Solexa;

  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];

  size_t read_count = 0;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      const string score_line(buffer);
      check_and_add_cheat(score_format, max_diffs, score_line, read_width, 
			  fast_reads_left, fast_reads_right, 
			  read_words_l, read_words_r, read_index_l, read_index_r, read_count);
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}
 
 
static void
check_and_add_cheat(const FASTQScoreType score_format, const size_t max_diffs,
		    const string &score_line, size_t &read_width, 
		    vector<FastReadQuality> &fast_reads_left, 
		    vector<FastReadQuality> &fast_reads_right, 
		    vector<size_t> &read_words_l, 
		    vector<size_t> &read_words_r, 
		    vector<size_t> &read_index_l, 
		    vector<size_t> &read_index_r, 
		    size_t &read_count) {
  
  // parse the score line
  vector<string> parts;
  smithlab::split_whitespace(score_line, parts);
  if (parts.size() % smithlab::alphabet_size != 0)
    throw SMITHLABException("bad format:\n" + score_line);

  // check the read width
  const size_t half_width = (parts.size()/2)/smithlab::alphabet_size;
  if (read_width == 0) {
    if (half_width*2 != parts.size()/smithlab::alphabet_size)
      throw SMITHLABException("PE reads not even length (must specify length)");
    read_width = half_width;
  }
  else if (parts.size()/smithlab::alphabet_size < 2*read_width)
    throw SMITHLABException("Incorrect read width");
  if (read_count == 0)
    FastReadQuality::set_read_width(read_width);
  
  // convert to numerical values
  const size_t lim = parts.size()/smithlab::alphabet_size;
  vector<vector<double> > error_probs(lim, vector<double>(smithlab::alphabet_size));
  for (size_t i = 0; i < lim; ++i)
    for (size_t j = 0; j < smithlab::alphabet_size; ++j)
      error_probs[i][j] = atof(parts[i*smithlab::alphabet_size + j].c_str());
  
  // convert to probability
  size_t bad_count = 0;
  for (size_t i = 0; i < lim; ++i) {
    for (size_t j = 0; j < smithlab::alphabet_size; ++j)
      error_probs[i][j] = 
	quality_score_to_error_probability(score_format, error_probs[i][j]);
    bad_count += (*min_element(error_probs[i].begin(), error_probs[i].end()) > 0.5);
  }
  
  const bool good_read = (bad_count <= 2*max_diffs);
  if (good_read) {
    fast_reads_left.push_back(FastReadQuality(error_probs.begin(), error_probs.begin() + read_width));
    fast_reads_right.push_back(FastReadQuality(error_probs.begin() + half_width,
					       error_probs.begin() + half_width + read_width));
    string read;
    for (size_t i = 0; i < error_probs.size(); ++i)
      read += int2base(min_element(error_probs[i].begin(), 
				   error_probs[i].end()) - error_probs[i].begin());
    string right_read(read.begin() + half_width,
		      read.begin() + half_width + read_width);
    string right_left(read.begin(), read.begin() + read_width);
    //  NO NEED TO REVCOMP AGAIN!! revcomp_inplace(right_read);
    read_words_r.push_back(get_read_word_cheat(right_read));
    read_words_l.push_back(get_read_word_cheat(right_left));
    read_index_r.push_back(read_count);
    read_index_l.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_prb_file_cheat(const string &filename, const size_t max_diffs,
			       size_t &read_width, 
			       vector<FastReadQuality> &fast_reads_left,
			       vector<FastReadQuality> &fast_reads_right,
			       vector<size_t> &read_words_l, 
			       vector<size_t> &read_words_r, 
			       vector<size_t> &read_index_l,
			       vector<size_t> &read_index_r) {
  FASTQScoreType score_format = FASTQ_Solexa;
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw SMITHLABException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];

  size_t read_count = 0;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw SMITHLABException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      const string score_line(buffer);
      check_and_add_cheat(score_format, max_diffs, score_line, read_width, 
			  fast_reads_left, fast_reads_right,
			  read_words_l, read_words_r, read_index_l, read_index_r, read_count);
    }
    in.peek();
  }
  if (fast_reads_left.empty())
    throw SMITHLABException("no high-quality reads in file:\"" + filename + "\"");
}
