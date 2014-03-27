
#include <vector>

#include "gmock/gmock.h"

#include "stoufferliptak.hpp"

using ::testing::Eq; using ::testing::DoubleNear;

using std::vector;

const double max_abs_error = 0.001;

TEST(StoufferLiptak, WorksWithoutCorrelationMatrix) {
  vector<double> p_vals = {0.1, 0.2, 0.8, 0.12, 0.011};

  ASSERT_THAT(stouffer_liptak_zaykin(p_vals), DoubleNear(0.0168, max_abs_error));
  
  p_vals = {0.5, 0.5, 0.5, 0.5, 0.5};
  ASSERT_THAT(stouffer_liptak_zaykin(p_vals), DoubleNear(0.5, max_abs_error));

  p_vals = {0.5, 0.1, 0.5, 0.5, 0.5};
  ASSERT_THAT(stouffer_liptak_zaykin(p_vals), DoubleNear(0.283, max_abs_error));

  p_vals = {0.5, 0.1, 0.1, 0.1, 0.5};
  ASSERT_THAT(stouffer_liptak_zaykin(p_vals), DoubleNear(0.042, max_abs_error));

  p_vals = {1, 0.1, 0.1};
  ASSERT_THAT(stouffer_liptak_zaykin(p_vals), DoubleNear(0.748, max_abs_error));
  
  p_vals = {1, 1, 1};
  ASSERT_THAT(stouffer_liptak_zaykin(p_vals), DoubleNear(1.0, max_abs_error));
  
  p_vals = {0, 0};
  ASSERT_THAT(stouffer_liptak_zaykin(p_vals), DoubleNear(0, max_abs_error));
}

TEST(StoufferLiptak, WorksWithPositiveDefiniteCorrelationMatrix) {
  vector<double> pvals = {0.1, 0.2};
  vector< vector<double> > corr_matrix = {{0.1, 0.2}, {0.1, 0.4}};
  ASSERT_THAT(stouffer_liptak_zaykin(pvals, corr_matrix), 
              DoubleNear(0.085, max_abs_error));
  
  pvals = {0.1, 0.2, 0.3};
  corr_matrix = {{0.1, 0.2, 0.1}, {0.1, 0.4, 0.5}, {0.1, 0.2, 0.3}};
  ASSERT_THAT(stouffer_liptak_zaykin(pvals, corr_matrix), 
              DoubleNear(0.109, max_abs_error));
              
  pvals = {1, 1, 1};
  corr_matrix = {{0.1, 0.2, 0.1}, {0.1, 0.4, 0.5}, {0.1, 0.2, 0.3}};
  ASSERT_THAT(stouffer_liptak_zaykin(pvals, corr_matrix), 
              DoubleNear(1.0, max_abs_error));

  pvals = {0, 0, 0.9};
  corr_matrix = {{0.1, 0.2, 0.1}, {0.1, 0.4, 0.5}, {0.1, 0.2, 0.3}};
  ASSERT_THAT(stouffer_liptak_zaykin(pvals, corr_matrix), 
              DoubleNear( 0.002, max_abs_error));
}

TEST(StoufferLiptak, WorksWithNotPositiveDefiniteCorrelationMatrix) {
  vector<double> pvals = {0.1, 0.2};
  vector< vector<double> > corr_matrix = {{0.1, 0.2}, {0.3, 0.4}};
  ASSERT_THAT(stouffer_liptak_zaykin(pvals, corr_matrix), 
              DoubleNear(0.085, max_abs_error));
}

TEST(StoufferLiptak, CaseFromMouseFrontalCortex) {
  vector<double> pvals = {0.00309547, 0.999, 0.0211285};
  vector< vector<double> > corr_matrix = {  {1,        0.486335, 0}, 
                                            {0.486335, 1,        0.520099}, 
                                            {0,        0.520099, 1} };

  //std::cerr << "With adjustment: " 
  //          << stouffer_liptak_zaykin(pvals, corr_matrix) 
  //          << std::endl;
  //std::cerr << "Without adjustment: " 
  //          << stouffer_liptak_zaykin(pvals) 
  //          << std::endl;
  
  ASSERT_THAT(stouffer_liptak_zaykin(pvals, corr_matrix), 
              DoubleNear(0.2267, max_abs_error));
}
