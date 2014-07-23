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
 
#include <sstream>
#include <string>

#include "gmock/gmock.h"

#include "design.hpp"

using ::testing::Eq; using ::testing::ElementsAre;
using ::testing::Ne;

using std::string;  using std::istringstream;
using std::ostringstream;

TEST(a_design, sets_number_of_samples_and_factors_to_zero_by_default) {
  Design design;
  ASSERT_THAT(design.num_factors(), Eq(0));
  ASSERT_THAT(design.num_samples(), Eq(0));
}

TEST(a_design, parses_factor_names_during_initialization) {
  istringstream iss( "f1 f2\n"
                    "s1 1 1\n"
                    "s2 1 0" );
  
  Design design(iss);
  
  ASSERT_THAT(design.factor_names(), ElementsAre("f1", "f2"));
  ASSERT_THAT(design.num_factors(), Eq(2));
}

TEST(a_design, parses_sample_names_during_initialization) {
  istringstream iss( "f1 f2\n"
                    "s1 1 1\n"
                    "s2 1 0"  );
  
  Design design(iss);
  ASSERT_THAT(design.sample_names(), ElementsAre("s1", "s2"));
}

TEST(a_design, parses_matrix_from_file) {
  istringstream iss(  "f1 f2\n"
                    "s1 1  1\n"
                    "s2 1  0" );
  
  Design design(iss);
  ASSERT_THAT(design.matrix().size(), Eq(2));
  ASSERT_THAT(design.matrix()[0], ElementsAre(1, 1));
  ASSERT_THAT(design.matrix()[1], ElementsAre(1, 0));
}

TEST(a_design, parses_matrix_from_file_disregarding_extra_lines) {
  istringstream iss(  "f1 f2\n"
                    "s1 1  1\n"
                    "s2 1  0\n"
                    "\n" );
  
  Design design(iss);
  ASSERT_THAT(design.matrix().size(), Eq(2));
  ASSERT_THAT(design.matrix()[0], ElementsAre(1, 1));
  ASSERT_THAT(design.matrix()[1], ElementsAre(1, 0));
}

TEST(a_design, verifies_that_factor_levels_are_binary) {
  istringstream iss(  "f1 f2\n"
                    "s1 1  5\n"
                    "s2 1  0");
  
  ASSERT_ANY_THROW(Design design(iss));
}

TEST(a_design, verifies_factor_name_exists_for_every_matrix_columns) {
  istringstream iss("f1\n"
                    "s1 1 1\n"
                    "s2 1 0");

  ASSERT_ANY_THROW(Design design(iss));
}

TEST(a_design, can_access_elements_through_function_call_operator) {
  istringstream iss(  "f1 f2\n"
                    "s1 1  1\n"
                    "s2 1  0");
                    
  Design design(iss);
  ASSERT_THAT(design(0, 0), Eq(1.0));
  ASSERT_THAT(design(1, 1), Eq(0.0));
}

TEST(a_design, outputs_its_string_representation) {
  string encoding = "f1\tf2\ns1\t1\t1\ns2\t1\t0";
  istringstream iss(encoding);
  
  Design design(iss);
  ostringstream oss;
  oss << design;
  ASSERT_THAT(oss.str(), Eq(encoding + '\n'));
}

TEST(a_design, removes_factor) {
  istringstream iss(  "f1 f2\n"
                    "s1 1  1\n"
                    "s2 1  0");
                    
  Design design(iss);
  design.remove_factor(0);
  ostringstream oss;
  oss << design;
  ASSERT_THAT(oss.str(), Eq("f2\ns1\t1\ns2\t0\n"));
}

TEST(a_design, removes_factor_by_factor_name) {
  istringstream iss( "f1 f2\n"
                     "s1 1 1\n"
                     "s2 1 0" );
  Design design(iss);
  
  design.remove_factor_name("f1");
  ostringstream oss;
  oss << design;
  ASSERT_THAT(oss.str(), Eq("f2\ns1\t1\ns2\t0\n"));
}

TEST(a_design, throws_exception_when_removing_nonexistent_factor_name) {
  istringstream iss(   "f1 f2\n"
                     "s1 1  1\n"
                     "s2 1  0");
  Design design(iss);
  ASSERT_ANY_THROW(design.remove_factor_name("I don't exist"));
}

TEST(a_design, equals_to_another_design_initialized_from_equal_string) {
  istringstream iss_a (   "f1 f2\n"
                        "s1 1  1\n"
                        "s2 1  0" );
  Design design_a(iss_a);
  
  istringstream iss_b( "f1 f2\n"
                      "s1 1 1\n"
                      "s2 1 0" );
  Design design_b(iss_b);
  
  ASSERT_THAT(design_a, Eq(design_b));
}

TEST(a_design, after_dropping_factor_design_is_not_equal_to_original) {
  string encoding =   "f1 f2\n"
                    "s1 1  1\n"
                    "s2 1  0";
  istringstream iss(encoding);
  
  Design design(iss);
  
  Design reduced_design = design;
  reduced_design.remove_factor_name("f2");
  
  ASSERT_THAT(reduced_design, Ne(design));
}
