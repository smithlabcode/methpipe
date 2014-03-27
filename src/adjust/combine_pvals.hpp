
#ifndef COMBINE_PVALS_HPP_
#define COMBINE_PVALS_HPP_

#include <vector>

#include "locus.hpp"
#include "bin_for_distance.hpp"

void combine_pvals(std::vector<LocusIterator> &loci_iterators, 
                    const BinForDistance &bin_for_distance);

#endif //COMBINE_PVALS_HPP_
