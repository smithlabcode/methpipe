
#include "cpg.hpp"
#include "merge.hpp"

using std::vector; using std::istream;
using std::ostream; using std::string;

void
merge_methylomes(vector<string> names, 
                  vector<istream*> methylomes, 
                  ostream &count_table) {
  
  vector<string>::const_iterator it = names.begin();
  
  count_table << *it++;
  
  while (it != names.end())
    count_table << "\t" << *it++;
  
  while (true) {
    vector<istream*>::iterator meth_it = methylomes.begin();
    string encoding;
    getline(*(*meth_it), encoding);
    
    if (encoding.empty())
      break;
      
    count_table << "\n";
    
    Cpg cpg(encoding);
    
    count_table << cpg.chrom() << ":" << cpg.locus() << ":" << cpg.locus() + 1
                << "\t" << cpg.total() << "\t" << cpg.meth();
    
    meth_it++;
    while(meth_it != methylomes.end()) {
      getline(*(*meth_it), encoding);
      Cpg cpg(encoding);
      
      count_table << "\t" << cpg.total() << "\t" << cpg.meth();
      meth_it++;
    }
  }
}