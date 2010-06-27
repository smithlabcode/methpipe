/*    sortreads: a program for sorting read sequences relative to an order
 *    given in a BED file (i.e. mapping locations)
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

#include <string>
#include <vector>
#include <unistd.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"

#include <sys/time.h>
#include <sys/resource.h>

using std::tr1::unordered_map;

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::numeric_limits;
using std::ofstream;

static const int INPUT_BUFFER_SIZE = 10000;

static void
skip_non_whitespace(const char *buffer, size_t &pos) {
  while (buffer[pos] != '\0' && buffer[pos] != '\t' && buffer[pos] != ' ') ++pos;
  if (buffer[pos] == '\0')
    throw RMAPException("bad format:\n" + string(buffer));
  while (buffer[pos] != '\0' && (buffer[pos] == '\t' || buffer[pos] == ' ')) ++pos;
  if (buffer[pos] == '\0')
    throw RMAPException("bad format:\n" + string(buffer));
}

static void
extract_read_name(const char *buffer, string &read_name) {
  size_t pos = 0;
  skip_non_whitespace(buffer, pos);
  skip_non_whitespace(buffer, pos);
  skip_non_whitespace(buffer, pos);
  const size_t start_pos = pos;
  while (buffer[pos] != '\0' && buffer[pos] != '\t' && buffer[pos] != ' ') ++pos;
  if (buffer[pos] == '\0')
    throw RMAPException("bad format:\n" + string(buffer));
  read_name = string(buffer + start_pos, buffer + pos);
}

static void
get_read_names(string filename, unordered_map<string, size_t> &read_name_index) {
  vector<char> buffer(INPUT_BUFFER_SIZE);
  
  std::ifstream in(filename.c_str());
  if (!in) throw RMAPException("cannot open input file " + filename);
  
  size_t read_count = 0;
  while (!in.eof()) {
    in.getline(&buffer[0], INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE - 1)
	throw RMAPException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      string read_name;
      extract_read_name(&buffer[0], read_name);
      read_name_index[read_name] = read_count++;
    }
    in.peek();
  }
  if (read_name_index.empty())
    throw RMAPException("no mapping locations in file:\"" + filename + "\"");
}

static size_t
get_file_id(const size_t tmp_file_size, const size_t idx) {
  return static_cast<size_t>(std::floor(idx*1.0/tmp_file_size));
}

inline bool
is_fasta_name_line(size_t line_count) {
  return ((line_count & 1ul) == 0);
}

static void
partition_fasta_file(const bool VERBOSE, string filename, 
		     const unordered_map<string, size_t> &read_name_index,
		     vector<ofstream*> &outfiles, const size_t tmp_file_size) {
  vector<char> buffer(INPUT_BUFFER_SIZE);
  
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw RMAPException("cannot open input file " + filename);
  
  const size_t filesize = get_filesize(filename);
  const size_t one_percent = filesize/100;
  
  size_t line_count = 0, file_id = std::numeric_limits<size_t>::max();
  size_t prev_total_read = 0;  
  string read_name;
  
  if (VERBOSE) cerr << "\r0%";
  while (!in.eof()) {
    in.getline(&buffer[0], INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE - 1)
	throw RMAPException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      if (is_fasta_name_line(line_count)) {
	read_name = string(&buffer[0] + 1);
	const unordered_map<string, size_t>::const_iterator 
	  idx(read_name_index.find(read_name));
	if (idx != read_name_index.end())
	  file_id = get_file_id(tmp_file_size, idx->second);
	else file_id = std::numeric_limits<size_t>::max();
      }
      else if (file_id != std::numeric_limits<size_t>::max())
	(*outfiles[file_id]) << ">" << read_name << "\n" 
			     << &buffer[0] << '\n';
      ++line_count;
    }
    const size_t total_read = in.tellg();
    if (VERBOSE && (total_read - prev_total_read > one_percent)) {
      cerr << '\r' << percent(total_read, filesize) << '%';
      prev_total_read = total_read;
    }
    in.peek();
  }
  if (VERBOSE) cerr << "\r100%" << endl;
  
}

inline bool
is_fastq_name_line(size_t line_count) {
  return ((line_count & 3ul) == 0ul);
}
inline bool
is_fastq_seq_line(size_t line_count) {
  return ((line_count & 3ul) == 1ul);
}
inline bool
is_fastq_score_line(size_t line_count) {
  return ((line_count & 3ul) == 3ul);
}

static void
partition_fastq_file(const bool VERBOSE, string filename,
		     const unordered_map<string, size_t> &read_name_index,
		     vector<ofstream*> &outfiles, const size_t tmp_file_size) {
  vector<char> buffer(INPUT_BUFFER_SIZE);
  
  if (VERBOSE)
    cerr << "PROCESSING: " << filename << endl;
  
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw RMAPException("cannot open input file " + filename);
  
  const size_t filesize = get_filesize(filename);
  const size_t one_percent = filesize/100;
  
  size_t line_count = 0, file_id = std::numeric_limits<size_t>::max();
  size_t prev_total_read = 0;  
  string read_name, read_seq;
  
  if (VERBOSE) cerr << "\r0%";
  while (!in.eof()) {
    in.getline(&buffer[0], INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE - 1)
	throw RMAPException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      if (is_fastq_name_line(line_count)) {
	read_name = string(&buffer[0] + 1);
	const unordered_map<string, size_t>::const_iterator 
	  idx(read_name_index.find(read_name));
	if (idx != read_name_index.end())
	  file_id = get_file_id(tmp_file_size, idx->second);
	else file_id = std::numeric_limits<size_t>::max();
      }
      else if (is_fastq_seq_line(line_count))
	read_seq = string(&buffer[0]);
      else if (is_fastq_score_line(line_count) &&
	       file_id != std::numeric_limits<size_t>::max()) {
	(*outfiles[file_id]) << '@' << read_name << '\n'
			     << read_seq << "\n+\n"
			     << &buffer[0] << '\n';
      }
      ++line_count;
    }
    const size_t total_read = in.tellg();
    if (VERBOSE && (total_read - prev_total_read > one_percent)) {
      cerr << '\r' << percent(total_read, filesize) << '%';
      prev_total_read = total_read;
    }
    in.peek();
  }
  if (VERBOSE) cerr << "\r100%" << endl;
}


static void
relative_sort_reads_fasta(const bool VERBOSE,
			  const bool KEEP_TEMP_FILES,
			  const size_t read_len,
			  const string &filename, 
			  const unordered_map<string, size_t> &read_name_index,
			  ofstream &out) {
  
  if (VERBOSE)
    cerr << "[READING]";
  vector<string> names, sequences;
  read_fasta_file(filename.c_str(), names, sequences);
  if (!KEEP_TEMP_FILES)
    remove(filename.c_str());
  
  if (VERBOSE)
    cerr << "[SORTING]";
  vector<pair<size_t, size_t> > sorter(names.size());
  for (size_t i = 0; i < names.size(); ++i) {
    const unordered_map<string, size_t>::const_iterator 
      idx(read_name_index.find(names[i]));
    if (idx == read_name_index.end())
      throw RMAPException("could not find read in temp file");
    sorter[i] = make_pair(idx->second, i);
  }
  sort(sorter.begin(), sorter.end());
  
  if (VERBOSE)
    cerr << "[WRITING]";
  for (size_t i = 0; i < names.size(); ++i) {
    const size_t idx = sorter[i].second;
    out << '>' << names[idx] << '\n' << sequences[idx].substr(0, read_len) << '\n';
  }
}


static void
relative_sort_reads_fastq(const bool VERBOSE,
			  const bool KEEP_TEMP_FILES,
			  const size_t read_len,
			  const string &filename, 
			  const unordered_map<string, size_t> &read_name_index,
			  ofstream &out) {
  if (VERBOSE)
    cerr << "[READING]";
  vector<string> names, sequences, scores;
  read_fastq_file(filename.c_str(), names, sequences, scores);
  if (!KEEP_TEMP_FILES)
    remove(filename.c_str());
  
  if (VERBOSE)
    cerr << "[SORTING]";
  vector<pair<size_t, size_t> > sorter(names.size());
  for (size_t i = 0; i < names.size(); ++i) {
    const unordered_map<string, size_t>::const_iterator 
      idx(read_name_index.find(names[i]));
    if (idx == read_name_index.end())
      throw RMAPException("could not find read in temp file");
    sorter[i] = make_pair(idx->second, i);
  }
  sort(sorter.begin(), sorter.end());
  
  if (VERBOSE)
    cerr << "[WRITING]";
  for (size_t i = 0; i < names.size(); ++i) {
    const size_t idx = sorter[i].second;
    const size_t curr_len = sequences[idx].length();
    if (sequences[idx].length() > read_len)
      out << '@' << names[idx] << '\n' 
	  << sequences[idx].substr(0, read_len) << '\n'
	  << "+\n" 
	  << scores[idx].substr(0, read_len) << '\n';
    else
      out << '@' << names[idx] << '\n' 
	  << sequences[idx] + string(read_len - curr_len, 'N') << '\n'
	  << "+\n" 
	  << scores[idx] + string(read_len - curr_len, 'B') << '\n';
  }
}


int 
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string mapped_file;
    string outfile;
    size_t tmp_file_size = 100000ul;
    string tmp_dir = ".";

    bool KEEP_TEMP_FILES = false;
  
    size_t read_len = std::numeric_limits<size_t>::max();
  
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], "a program for sorting read sequences "
			   "relative to an order given in a BED file (i.e. "
			   "mapping locations)cpgcaller2",
			   "<fasta-reads-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      true, outfile);
    opt_parse.add_opt("mapped", 'm', "file of mapped locations", 
		      true , mapped_file);
    opt_parse.add_opt("size", 's', "size of temporary files", 
		      false , tmp_file_size);
    opt_parse.add_opt("width", 'w', "read width", 
		      false , read_len);
    opt_parse.add_opt("dir", 'd', "directory for temporary files", 
		      false , tmp_dir);
    opt_parse.add_opt("keep", 'K', "keep temporary files", false, KEEP_TEMP_FILES);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    vector<string> reads_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[LOADING MAPPED LOCATIONS]" << endl;
    unordered_map<string, size_t> read_name_index;
    get_read_names(mapped_file, read_name_index);
    if (VERBOSE)
      cerr << "[LOADED " << read_name_index.size() << " MAPPED READS]" << endl;
    const size_t n_files = get_file_id(tmp_file_size, read_name_index.size() - 1) + 1;
    
    const string pid = rmap::toa(getpid());
    vector<string> filenames;
    for (size_t i = 0; i < n_files; ++i)
      filenames.push_back(path_join(tmp_dir, "." + rmap::toa(i) + 
				    "." + pid + ".tmp"));
    

    if (VERBOSE)
      cerr << "CREATING: " << filenames.front() << " -- "
	   << filenames.back() << endl;
    vector<ofstream*> outfiles;
    for (size_t i = 0; i < filenames.size(); ++i) {
      outfiles.push_back(new ofstream(filenames[i].c_str()));
      if (!(*outfiles.back()))
	throw RMAPException("cannot open input file " + filenames[i]);
    }
    
    if (VERBOSE)
      cerr << "PARTITIONING READS" << endl;

    for (size_t i = 0; i < reads_files.size(); ++i) {
      const bool FASTQ = is_fastq(reads_files[i]);
      if (VERBOSE)
	cerr << "READS FILE FORMAT: " << reads_files[i] 
	     << "=" << ((FASTQ) ? "FASTQ" : "FASTA") << endl;
      if (FASTQ)
	partition_fastq_file(VERBOSE, reads_files[i], read_name_index, 
			     outfiles, tmp_file_size);
      else
	partition_fasta_file(VERBOSE, reads_files[i], read_name_index, 
			     outfiles, tmp_file_size);
      if (VERBOSE)
	cerr << "DONE: " << reads_files[i] << endl;
    }
    for (size_t i = 0; i < outfiles.size(); ++i)
      delete outfiles[i];
    
    ofstream out(outfile.c_str());
    for (size_t i = 0; i < n_files; ++i) {
      const bool FASTQ = is_fastq(filenames[i]);
      if (VERBOSE) cerr << "[SORTING=" << filenames[i] << "]";
      if (FASTQ)
	relative_sort_reads_fastq(VERBOSE, KEEP_TEMP_FILES, read_len,
				  filenames[i], read_name_index, out);
      else relative_sort_reads_fasta(VERBOSE, KEEP_TEMP_FILES, read_len,
				     filenames[i], read_name_index, out);
      if (VERBOSE)
	cerr << "[DONE]" << endl;
    }
    if (VERBOSE) cerr << endl;
  }
  catch (const RMAPException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
