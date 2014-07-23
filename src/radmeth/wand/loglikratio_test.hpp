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

// This file declares the loglikratio_test function which performs the 
// log-likelihood ratio test.

#ifndef LOGLIKRATIO_TEST_
#define LOGLIKRATIO_TEST_

// Performs a log-likelihood ratio test given likelihoods of the full and
// the null models. It is assumed that the null model has one factor fewer
// than the full model and so the number of degrees of freedom is 1.
double loglikratio_test(double null_loglik, double full_loglik);

#endif // LOGLIKRATIO_TEST_
