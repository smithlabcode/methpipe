
#ifndef STOUFFERLIPTAK_HPP_
#define STOUFFERLIPTAK_HPP_

#include<vector>

double 
stouffer_liptak(std::vector<double> &p_vals, 
                const std::vector< std::vector<double> > &cor_matrix = 
                                        std::vector< std::vector<double> >());

double 
stouffer_liptak_zaykin(std::vector<double> &pvals, 
                       const std::vector< std::vector<double> > &cor_matrix = 
                                        std::vector< std::vector<double> >());
  
#endif //STOUFFERLIPTAK_HPP_

