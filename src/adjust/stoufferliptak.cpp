
#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include "stoufferliptak.hpp"

using std::for_each;  using std::vector;
using std::cout;  using std::endl;
using std::back_inserter;

static double 
pval_to_zscores(double x) {
  if (x == 1)
    x = 0.9999;
  if (x == 0)
    x = 1 - 0.9999;
  return gsl_cdf_ugaussian_Pinv(1 - x);
}

gsl_matrix* 
matrix_to_gslmatrix(const vector< vector<double> > &matrix) {
  //assuming the matrix is not ragged
  size_t nrows = matrix.size();
  size_t ncols = matrix[0].size();
  gsl_matrix* gmatrix = gsl_matrix_alloc(nrows, ncols);
  for (size_t row = 0; row < nrows; ++row)
    for (size_t col = 0; col < ncols; ++col)
      gsl_matrix_set(gmatrix, row, col, matrix[row][col]);
  return gmatrix;
}

gsl_vector* 
vector_to_gsl_vector(const vector<double> &vector) {
  //assuming the matrix is not ragged
  size_t nelems = vector.size();
  gsl_vector* gvector = gsl_vector_alloc(nelems);
  for (size_t ind = 0; ind < nelems; ++ind)
    gsl_vector_set(gvector, ind, vector[ind]);
  return gvector;
}

double
stouffer_liptak(vector<double> &pvals,
                    const vector< vector<double> > &corr_matrix) {
  
  const size_t num_pvals = pvals.size();

  vector<double> zscores;
  
  transform(pvals.begin(), pvals.end(), back_inserter(zscores), 
            pval_to_zscores);
  
  gsl_vector* gsl_zsocres = vector_to_gsl_vector(zscores); //freed
  
  if (!corr_matrix.empty()) {
    gsl_set_error_handler_off();
    gsl_matrix *gsl_corr_matrix = matrix_to_gslmatrix(corr_matrix); //freed
    
    //Compute Cholesky decomposition.
    int err_status = gsl_linalg_cholesky_decomp(gsl_corr_matrix);
   
    //If successful, get Cholesky matrix by zeroing out 
    //above the main diagonal.
    if (err_status != GSL_EDOM) {
      for (size_t row = 0; row < num_pvals; ++row)
        for (size_t col = row + 1; col < num_pvals; ++col)
          gsl_matrix_set(gsl_corr_matrix, row, col, 0);
    
      //get LU-decomposition in order to invert Cholesky matrix 
      gsl_matrix *chol_inv = gsl_matrix_alloc(num_pvals, num_pvals); //freed
      gsl_permutation *p = gsl_permutation_alloc (num_pvals); //freed
      int signum;
      gsl_linalg_LU_decomp (gsl_corr_matrix, p, &signum);
      gsl_linalg_LU_invert (gsl_corr_matrix, p, chol_inv);

      //Matrix to store Cholesky decomposition
      //gsl_matrix *chol_inv = gsl_matrix_alloc(num_pvals, num_pvals);
      
      //freed
      gsl_vector* gsl_uncorrelated_zscores = gsl_vector_alloc(num_pvals); 
      gsl_blas_dgemv( CblasNoTrans, 1.0, chol_inv, gsl_zsocres, 0.0, 
                      gsl_uncorrelated_zscores );

      gsl_vector_memcpy(gsl_zsocres, gsl_uncorrelated_zscores);
      
      gsl_matrix_free(gsl_corr_matrix);
      gsl_matrix_free(chol_inv);
      gsl_permutation_free(p);
      gsl_vector_free(gsl_uncorrelated_zscores);
    }
  }
  
  double sum = 0;
  for (size_t ind = 0; ind < num_pvals; ++ind)
    sum += gsl_vector_get(gsl_zsocres, ind);
    
  gsl_vector_free(gsl_zsocres);
  
  return 1 - gsl_cdf_gaussian_P(sum/(std::sqrt(num_pvals)), 1);
}

double 
stouffer_liptak_zaykin(std::vector<double> &pvals, 
                       const std::vector< std::vector<double> > &cor_matrix) {
  double correction = 0;
  size_t num_pvals = pvals.size();
  for (size_t row_ind = 0; row_ind < cor_matrix.size(); ++row_ind)
    for (size_t col_ind = row_ind + 1; col_ind < cor_matrix.size(); ++col_ind)
      correction += cor_matrix[row_ind][col_ind];
  
  vector<double> zscores;
  
  transform(pvals.begin(), pvals.end(), back_inserter(zscores), 
            pval_to_zscores);
  
  double sum = 0;
  for (size_t ind = 0; ind < num_pvals; ++ind)
    sum += zscores[ind];
  
  double test_statistic = sum/sqrt(double(num_pvals) + 2*correction);
  
  return 1 - gsl_cdf_gaussian_P(test_statistic, 1);
}

