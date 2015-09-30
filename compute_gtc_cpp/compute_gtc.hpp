// Header for GTC algorithm and support functions
#include<iostream>
#include<iomanip> 
#include<fstream>
#include<sstream>
#include<string>
#include<math.h>
#include<vector>
#include<omp.h>

#ifndef __GLOG_H_INCLUDED
#include <glog/logging.h>
#include <gflags/gflags.h>
#define __GLOG_H_INCLUDED
#endif

#ifndef __BOOST_TIMER_H_INCLUDED
#include<boost/timer/timer.hpp>
#define __BOOST_TIMER_H_INCLUDED
#endif

#ifndef __COMPLEX_H_INCLUDED
#define __COMPLEX_H_INCLUDED
#include <complex.h>
#endif

#ifndef __CS_H_INCLUDED
#define __CS_INCLUDED
#define CS_LONG
#define CS_COMPLEX
#include <cs.h>
#endif

#ifndef __KLU_H_INCLUDED
#define __KLU_H_INCLUDED
#include <klu.h>
#endif


CS_INT thevenin_voltage_source( cs * y_nc, cs * y_link_up, cs * y_link_lo, cs * y_vc, std::vector<double>& zth_vc_real, std::vector<double>& zth_vc_imag, std::vector<double>& gtc_vc_real, std::vector<double>& gtc_vc_imag); 
   // Obtain Thevenin Equivalents seen from voltage sources in a network.
   // IN: y_nc is an (n-k by n-k)  admittance matrix with all voltage sources short-circuited.
   // IN: y_vc is an (k by k)  admittance matrix with all the nodes of constant voltage.
   // IN: y_link_up is a (n-k by k) upper right admittance matrix linking y_nc with y_vc.
   // IN: y_link_lo is a (k by n-k) lower left admittance matrix linking y_nc with y_vc.
   // IN/OUT: zth_vc_real is a vector of length k holding the real part of Thevenin Impedances in order according to y_vc.
   // IN/OUT: zth_vc_imag is a vector of length k holding the imaginary part of Thevenin Impedances in order according to y_vc.
   // IN/OUT: gtc_vc_real is a vector of length k*k which at return holds the real parts of gtc coefficients ordered suchs that consecutive elements corresponds to the rows of the GTC matrix.
   // IN/OUT: gtc_vc_imag is a vector of length k*k which at return holds the imaginary parts of gtc coefficients ordered suchs that consecutive elements corresponds to the rows of the GTC matrix.
   // If either gtc_vc_real or gtc_vc_imag have size smaller than k*k thevenin_voltage_sources will return thevenin impedances only
  


cs_cl * cs_load_matrix(const std::string s_MatFileName); // read data file to sparse matrix

CS_INT cs_element(cs *A, CS_INT i, CS_INT j, cs_complex_t *c);// Extract single element from compressed column matrix

CS_INT cs_spcolumn(cs *A, CS_INT k, cs *b); // Extract sparse column vector from compressed column matrix:

cs *cs_vcat(cs *A, cs_complex_t *b, const CS_INT n); // vertical concatenation, sparse mat -dense row

cs * cs_extract_diags(cs *A); // extract diagonal of matrix A and store result in diagonal csc.

cs * cs_submatrix(cs *A, CS_INT row_a, CS_INT col_a, CS_INT row_b, CS_INT col_b);


