/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
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
 */

// GSL headers.
#include <sstream>
#include <vector>

// GSL headers.
#include <gsl/gsl_matrix_double.h>

// Google Mock headers.
#include "gmock/gmock.h"

// Local headers.
#include "design.hpp"
#include "gsl_fitter.hpp"
#include "regression.hpp"

using ::testing::Eq; using ::testing::ElementsAre;
using ::testing::DoubleNear; using ::testing::Test;

using std::istringstream; using std::vector;

vector<double> 
gslvec_to_stlvec(gsl_vector *gsl_vec) {
  vector<double> stl_vec;
  for (size_t ind = 0; ind < gsl_vec->size; ++ind)
    stl_vec.push_back(gsl_vector_get(gsl_vec, ind));
  return stl_vec;
}

const double max_abs_error = 0.001;

MATCHER_P(DoubleVectorNear, value, "") {
  if (arg.size() != value.size()) 
    return false;
  
  for (size_t ind = 0; ind < arg.size(); ++ind) {
    if (fabs(arg[ind] - value[ind]) >= max_abs_error)
      return false;
  }
  
  return true;
}

class Unfitted2x2RegressionWithGivenParameters: public Test {
public:
  istringstream os;
  Design *design;
  Regression *regression;
  gsl_vector *parameters;
  
  void SetUp() override {
    os.str("f1\tf2\ns1\t1\t1\ns2\t1\t0");
    design = new Design(os);
    regression = new Regression(*design);
    
    vector<size_t> response_total(4, 10);
    vector<size_t> response_meth = {2,	3,	7,	8};

    regression->set_response(response_total, response_meth);
    
    parameters = gsl_vector_calloc(3);
    gsl_vector_set(parameters, 0, 2.0);
    gsl_vector_set(parameters, 1, 1.5);
    gsl_vector_set(parameters, 2, 0.9);

  }
  void TearDown() override {
    delete design;
    design = nullptr;
    delete regression;
    regression = nullptr;
    gsl_vector_free(parameters);
    parameters = nullptr;
  } 
};

class Unfitted10x1Regression: public Test {
public:
  istringstream design_encoding;
  Design *design;
  Regression *regression;
  
  void SetUp() override {
    design_encoding.str(
      "f\ns0\t1\ns1\t1\ns2\t1\ns3\t1\ns4\t1\ns5\t1\ns6\t1\ns7\t1\ns8\t1\ns9\t1"
      );
    design = new Design(design_encoding);
    regression = new Regression(*design);
  
    vector<size_t> response_total(10, 20);
    vector<size_t> response_meth = { 0,  0,  6,  9,  7,  6,  2,  8, 10,  5 };
    
    regression->set_response(response_total, response_meth);
  }

  void TearDown() override {
    delete design;
    design = nullptr;
    delete regression;
    regression = nullptr;
  }
};

class Fitted8x2Regression: public Test {
public:
  istringstream design_encoding;
  Design *design;
  Regression *regression;

  void SetUp() override { 
    design_encoding.str("f0\tf1\ns0\t1\t0\ns1\t1\t0\ns2\t1\t0\ns3\t1\t0\n"
                                "s4\t1\t1\ns5\t1\t1\ns6\t1\t1\ns7\t1\t1");
    design = new Design(design_encoding);
    regression = new Regression(*design);
  
    vector<size_t> response_total(8, 15);
    vector<size_t> response_meth = { 8, 5, 2, 2, 13, 14, 8, 5};
  
    regression->set_response(response_total, response_meth);
    
    vector<double> initial_parameters = {0.0, 0.0, 0.2};
    gsl_fitter(*regression, initial_parameters);
  }
  
  void TearDown() override {
    delete design;
    design = nullptr;
    delete regression;
    regression = nullptr;
  }

};

TEST(a_regression, initializes_design_matrix) {
  istringstream os("f1\tf2\ns1\t1\t1\ns2\t1\t0");
  Design design(os);
  Regression regression(design);
  ASSERT_THAT(regression.design(), Eq(design));
}

TEST(a_regression, accepts_response) {
  istringstream os("f1\tf2\ns1\t1\t1\ns2\t1\t0");
  Design design(os);
  Regression regression(design);
  vector<size_t> response_total(4, 100);
  vector<size_t> response_meth = {64,	54,	72,	32};

  regression.set_response(response_total, response_meth);
  ASSERT_THAT(regression.response_total(), ElementsAre(100, 100, 100, 100));
  ASSERT_THAT(regression.response_meth(), ElementsAre(64, 54, 72, 32));
}

TEST_F(Unfitted2x2RegressionWithGivenParameters, computes_p) {
  ASSERT_THAT(regression->p(0, parameters), DoubleNear(0.9707, max_abs_error));
}

TEST_F(Unfitted2x2RegressionWithGivenParameters, computes_loglikelihood) {
  ASSERT_THAT(regression->loglik(parameters), 
              DoubleNear(-18.556, max_abs_error));
}

TEST_F(Unfitted2x2RegressionWithGivenParameters, computes_gradient) {  
  gsl_vector *grad = gsl_vector_calloc(3);
  regression->gradient(parameters, grad);
  vector<double> stl_grad = gslvec_to_stlvec(grad);
  gsl_vector_free(grad);
  
  vector<double> expected = {-1.77654, -0.962869, -0.935081};
  ASSERT_THAT(stl_grad, DoubleVectorNear(expected));  
}

/*
TEST_F(Unfitted2x2RegressionWithGivenParameters, 
                                          computes_negative_loglik_and_grad) {
  
  double loglik_val;
  gsl_vector *d_loglik_val = gsl_vector_calloc(3);
  regression->neg_loglik_and_gradient(parameters, &loglik_val, d_loglik_val);

  vector<double> stl_grad = gslvec_to_stlvec(d_loglik_val);
  vector<double> expected = {1.77654, 0.962869, 0.935081};
  gsl_vector_free(d_loglik_val);

  ASSERT_THAT(loglik_val, DoubleNear(18.556, max_abs_error));  

  ASSERT_THAT(stl_grad, DoubleVectorNear(expected));
}*/

TEST_F(Unfitted10x1Regression, finishes_fitting_in_allotted_iterations) {

  vector<double> initial_parameters = {0.0, 0.1};
  
  ASSERT_THAT(gsl_fitter(*regression, initial_parameters), Eq(true));
}

TEST_F(Unfitted10x1Regression, checks_number_of_initial_parameters) {
  
  vector<double> initial_parameters = {0.0};
  
  ASSERT_ANY_THROW(gsl_fitter(*regression, initial_parameters));
}

TEST_F(Unfitted10x1Regression, creates_initial_paramters_if_none_are_given) {
  ASSERT_THAT(gsl_fitter(*regression), Eq(true));
}

TEST(regression_10x1_design, estimates_regression_parameters_) {
  istringstream design_matrix_encoding(
    "f\ns0\t1\ns1\t1\ns2\t1\ns3\t1\ns4\t1\ns5\t1\ns6\t1\ns7\t1\ns8\t1\ns9\t1");
  Design design(design_matrix_encoding);
  Regression regression(design);
  
  vector<size_t> response_total(10, 15);
  vector<size_t> response_meth = { 0,  0,  6,  9,  7,  6,  2,  8, 10,  5 };
    
  regression.set_response(response_total, response_meth);
  gsl_fitter(regression);
    
  vector<double> regression_parameters = regression.fitted_parameters();
  vector<double> actual = {-0.679533, -1.28847};
  
  ASSERT_THAT(regression_parameters, DoubleVectorNear(actual));
}


TEST(regression_8x1_design, estimates_regression_parameters) {
  istringstream design_matrix_encoding(
    "f\ns0\t1\ns1\t1\ns2\t1\ns3\t1\ns4\t1\ns5\t1\ns6\t1\ns7\t1");
  Design design(design_matrix_encoding);
  Regression regression(design);
  
  vector<size_t> response_total(8, 15);
  vector<size_t> response_meth = { 2, 0, 14, 8, 15, 11, 10, 1 };
  
  regression.set_response(response_total, response_meth);
  
  vector<double> initial_parameters = {0.0, 0.2};
  
  regression.set_response(response_total, response_meth);
  gsl_fitter(regression, initial_parameters);
    
  vector<double> regression_parameters = regression.fitted_parameters();
  
  vector<double> expected = {-0.00160943, -0.02214};
  
  ASSERT_THAT(regression_parameters, DoubleVectorNear(expected));
}

TEST_F(Fitted8x2Regression, estimates_regression_parameters) {
  vector<double> regression_parameters = regression->fitted_parameters();
  vector<double> expected = {-0.868078, 1.61086, -1.8563};
  ASSERT_THAT(regression_parameters, DoubleVectorNear(expected));
}

TEST_F(Fitted8x2Regression, estimates_distribution_parameters) {      
  vector<double> distribution_parameters = 
                                  regression->fitted_distribution_parameters();
  
  ASSERT_THAT(distribution_parameters.size(), Eq(9));
  
  vector<double> expected = { 0.2960, 0.2960, 0.2960, 0.2960, 
                              0.6776, 0.6776, 0.6776, 0.6776, 
                              0.1351 };
 
  ASSERT_THAT(distribution_parameters, DoubleVectorNear(expected));
}

TEST_F(Fitted8x2Regression, maximum_loglikelihood) {  
  ASSERT_THAT(regression->maximum_likelihood(), 
              DoubleNear(-69.8045, max_abs_error));
}

TEST_F(Fitted8x2Regression, estimates_min_methdiff) {
  ASSERT_THAT(regression->min_methdiff(1), 
              DoubleNear(0.381948, max_abs_error));
}

TEST_F(Fitted8x2Regression, estimates_log_fold_change) {
  ASSERT_THAT(regression->log_fold_change(1),
              DoubleNear(1.61086, max_abs_error));
}
