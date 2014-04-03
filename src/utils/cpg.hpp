
#ifndef CPG_HPP_
#define CPG_HPP_

#include <string>

// This file contains that definition of the Cpg class designed to 
// store information about a single CpG site.

class Cpg {
public:
  Cpg() {};
  Cpg(std::string encoding);
  void set_chrom(std::string chrom) { chrom_ = chrom; }
  std::string chrom() const { return chrom_; }
  void set_locus(size_t locus) { locus_ = locus; }
  size_t locus() const { return locus_; }
  void set_total(size_t total) { total_ = total; }
  size_t total() const { return total_; }
  void set_meth(size_t meth) { meth_ = meth; }
  size_t meth() const { return meth_; }
private:
  std::string chrom_;
  size_t locus_;
  size_t total_;
  size_t meth_;
};

#endif //CPG_HPP_