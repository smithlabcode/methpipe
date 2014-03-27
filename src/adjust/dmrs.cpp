
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "OptionParser.hpp"

#include "locus.hpp"
#include "merge.hpp"

using std::cout;  using std::endl;
using std::cerr;  using std::vector;
using std::string;  using std::pair;

int main(int argc, const char **argv) {
  try {
    /* FILES */
    string outfile;
    string bin_spec = "1:200:25";
    double cutoff = 0.01;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("dmrs", "a program to merge significantly "
                           "differentially methylated CpGs into DMRs",
                           "<bed-file-in-wand-format>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("cutoff", 'p', "P-value cutoff (default: 0.01)", 
		      false , cutoff);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string bed_filename = leftover_args.front();
    /*************************************************************************/
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    std::ifstream bed_file(bed_filename.c_str());
    
    if (!bed_file)
      throw "could not open file: " + bed_filename;
    
    merge(bed_file, out, cutoff);
        
  } catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
