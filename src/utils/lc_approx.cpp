/*    apxlc: a program for approximately and very quickly counting lines
 *
 *    Copyright (C) 2012 University of Southern California and
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
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"

using std::string;
using std::ios_base;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

static size_t
get_approx_line_count(const bool VERBOSE, const string &filename,
		      const size_t n_samples, size_t sample_size) {
  
  static const size_t megabyte = (1ul << 20);
  static const size_t kilobyte = (1ul << 10);
  
  const size_t filesize = get_filesize(filename);
  
  if (sample_size == 0)
    sample_size = std::min(megabyte/10, filesize/n_samples);
  
  const size_t increment = 
    std::floor((filesize - sample_size*n_samples)/
	       (n_samples - 1.0)) + sample_size;
  
  assert(filesize > n_samples && filesize > sample_size && 
	 filesize > n_samples*sample_size);
  
  if (VERBOSE) {
    cerr << "[PROCESSING FILE: " << filename << "]" << endl
	 << "[FILESIZE: " 
	 << static_cast<double>(filesize)/megabyte << "MB]" << endl
	 << "[CHUNK SIZE: " 
	 << static_cast<size_t>(1.0*sample_size/kilobyte) << "KB]" << endl
	 << "[NUM CHUNKS: " << n_samples << "]" << endl
	 << "[TOTAL SAMPLE: " 
	 << (1.0*n_samples*sample_size)/megabyte << "MB]" << endl;
  }
  std::ifstream in(filename.c_str(), ios_base::binary);
  if (!in)
    throw SMITHLABException("cannot open input file " + string(filename));
  
  vector<char> buffer(sample_size);
  double total_lines = 0.0;
  for (size_t i = 0; i < filesize && in.good(); i += increment) {
    in.seekg(i, ios_base::beg);
    in.read(&buffer.front(), sample_size);
    if (in.good())
      total_lines += (0.5 + count(buffer.begin(), buffer.end(), '\n'));
  }
  return (filesize*total_lines)/(n_samples*sample_size);
}



int 
main(int argc, const char **argv) {
  try {

    size_t n_samples = 100;
    size_t sample_size = 0;
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
			   "approximate line counting in large files",
			   "<file1> <file2> ..." );
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("samples", 'n', "number of samples", false, n_samples);
    opt_parse.add_opt("size", 'z', "sample size (bytes)", false, sample_size);

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
    if (leftover_args.size() < 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_FAILURE;
    }
    vector<string> filenames(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/
    //////////////////////////////////////////////////////////////
    
    for (size_t i = 0; i < filenames.size(); ++i)
      cout << filenames[i] << "\t" 
	   << get_approx_line_count(VERBOSE, filenames[i],
				    n_samples, sample_size) << endl;
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
