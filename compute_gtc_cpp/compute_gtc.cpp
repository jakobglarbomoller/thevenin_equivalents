// Sources for GTC algorithm and support functions
#include "compute_gtc.hpp"


   cs_complex_t cs_one = cs_complex_t(1,0);
   cs_complex_t cs_zero = cs_complex_t(0,0);


CS_INT thevenin_voltage_source( cs * y_nc, cs * y_link_up, cs * y_link_lo, cs * y_vc, std::vector<double>& zth_vc_real, std::vector<double>& zth_vc_imag, std::vector<double>& gtc_vc_real, std::vector<double>& gtc_vc_imag){
   // Obtain Thevenin Equivalents seen from voltage sources in a network.
   // IN: y_nc is an (n-k by n-k)  admittance matrix with all voltage Sources short-circuited.
   // IN: y_vc is an (k by k)  admittance matrix with all the nodes of constant voltage.
   // IN: y_link_up is a (n-k by k) upper right admittance matrix linking y_nc with y_vc.
   // IN: y_link_lo is a (k by n-k) lower left admittance matrix linking y_nc with y_vc.
   // OUT: zth_vc_real is a vector of length k holding the real part of Thevenin Impedances in order according to y_vc.
   // OUT: zth_vc_imag is a vector of length k holding the imaginary part of Thevenin Impedances in order according to y_vc.
   // OUT: gtc_vc_real is a vector of length k*k which at return holds the real parts of gtc coefficients ordered suchs that consecutive elements corresponds to the rows of the GTC matrix.
   // OUT: gtc_vc_imag is a vector of length k*k which at return holds the imaginary parts of gtc coefficients ordered suchs that consecutive elements corresponds to the rows of the GTC matrix.
   // If either gtc_vc_real or gtc_vc_imag have size smaller than k*k thevenin_voltage_sources will return thevenin impedances only


   boost::timer::cpu_timer timer;
   boost::timer::cpu_times times;
   double d_wall_time;
   double d_cpu_time;
   timer.start();

   // check inputs
   if ((!CS_CSC(y_nc))||(!CS_CSC(y_link_up))||(!CS_CSC(y_link_lo))||(!CS_CSC(y_vc))){
      std::cout << "ERROR: admittance submatrices must be on cs compressed column format" << std::endl;
      return -1;
   }

   int nc_size = y_nc->n;
   int vc_size = y_vc->n;
   int numel_gtc = vc_size*vc_size;

   if((zth_vc_real.size() < vc_size) || (zth_vc_real.size() < vc_size)) {
      std::cout << "ERROR: container for real or imaginary Thevenin impedance unallocated. " << std::endl;
      return -1;
   }

   bool computeGTCs = true;
   if ((gtc_vc_real.size() < numel_gtc) || (gtc_vc_imag.size() < numel_gtc)) {
      computeGTCs = false;
      std::cout << "WARNING: GTC computations omitted. Output container not allocated." << std::endl;
   }

   int order = 2;     // 0: natural ordering (no column permutation, only row permutation: LU = PA)
   int qr = 0;        // 0: use LU factorization not QR factorization 
   double tol = 0.001;

      // sym. ordering return estimate on non-zeros in L and U S->unz, S->lnz and fill-in reducing order S->q
   css *y_nc_symbolic_ordering = cs_sqr(order, y_nc, qr);

      // LU decomposition
   csn *y_nc_lu_factorization = cs_lu(y_nc, y_nc_symbolic_ordering, tol);

   if ((y_nc_symbolic_ordering == NULL) || (y_nc_lu_factorization == NULL)) {
      std::cout << "ERROR: LU factorization failed" << std::endl;
      cs_sfree(y_nc_symbolic_ordering);
      cs_nfree(y_nc_lu_factorization);
      return -1;
   }

   cs_dropzeros(y_nc_lu_factorization->L);
   cs_dropzeros(y_nc_lu_factorization->U);
      // get row permutation
   CS_INT *row_permute, *col_permute;
      // get column permutation

  col_permute = y_nc_lu_factorization->pinv;  
//  row_permute = cs_pinv(y_nc_lu_factorization->pinv, nc_size);
  row_permute = cs_pinv(y_nc_lu_factorization->pinv, nc_size);  

      // permute y_link matrices
   cs *p_y_link_up = cs_permute(y_link_up, col_permute, NULL ,1);
   cs *p_y_link_lo = cs_permute(y_link_lo, NULL, row_permute, 1);


   if ((p_y_link_up == NULL) || (p_y_link_lo == NULL)){
      std::cout << "ERROR: permutation of y_link failed" << std::endl;
      cs_spfree(p_y_link_up);
      cs_spfree(p_y_link_lo);
      cs_sfree(y_nc_symbolic_ordering);
      cs_nfree(y_nc_lu_factorization);
      return -1;
   }

   // store p_y_link_lo in transposed form to exploit lower triangular solver
   cs * p_y_link_lo_t = p_y_link_lo;
   p_y_link_lo_t = cs_transpose(p_y_link_lo, -1);

   cs * L_nc_t = cs_transpose(y_nc_lu_factorization->L,-1);



   if ((p_y_link_lo_t == NULL)){
      std::cout << "ERROR: Transposition of U or y_link_lo failed" << std::endl;
      cs_spfree(p_y_link_up);
      cs_spfree(p_y_link_lo);
      cs_spfree(p_y_link_lo_t);
      cs_sfree(y_nc_symbolic_ordering);
      cs_nfree(y_nc_lu_factorization);
      return -1;
   }


   timer.stop();
   times = timer.elapsed();
   d_wall_time = times.wall;
   d_cpu_time = times.system+times.user;
   LOG(INFO) << "Time spent on LU-decomposition:              " << std::setw(8) << d_wall_time * 1e-9 << " s, cpu time: " << d_cpu_time * 1e-9 << std::endl;
   timer.start();


// forward and backward susbstitutions:

   cs *Y_eq_tripl, *Y_eq;
   Y_eq_tripl = cs_spalloc(vc_size,vc_size,2*sizeof(cs_complex_t),1,1);
   cs_entry(Y_eq_tripl,0,0,cs_zero);
   Y_eq = cs_compress(Y_eq_tripl);
   cs_spfree(Y_eq_tripl);

   double tollerance = 1e-13;
   int num_threads;



#pragma omp parallel 
{
   num_threads = omp_get_num_threads();
}

   cs * p_BL_sp_local[num_threads];
   cs * p_BU_sp_local[num_threads];
   cs * p_y_eq_local[num_threads];
   cs * p_y_eq_local_c[num_threads];

//   cs * y_eq_triplet;
//   y_eq_triplet = cs_spalloc(vc_size,vc_size,10*sizeof(cs_complex_t),1,1);

#pragma omp parallel
{
   int this_thread = omp_get_thread_num();

   // allocate variables used to solve solve BL^T = U^T \ Ylink_lo^T
//   cs * BL_sp_local;
   p_BL_sp_local[this_thread] = cs_spalloc(vc_size, nc_size, nc_size*sizeof(cs_complex_t), 1, 1);
   cs *u_nc_t_local = cs_transpose(y_nc_lu_factorization->U, -1);
   CS_INT xi_l[2*nc_size];

   // allocate variables used to solve BU = L \ Ylink_up
//   cs * BU_sp_local;
   p_BU_sp_local[this_thread] = cs_spalloc(nc_size, vc_size, nc_size*sizeof(cs_complex_t), 1, 1);
   cs * L_nc_local = cs_transpose(L_nc_t,-1);
   CS_INT xi_u[2*nc_size];

// allocate variables for subsequent multiplication
   cs * BL_local_c;
   cs * BU_local_c;

   // solve BL^T = U^T \ Ylink_lo^T
#pragma omp for schedule(static)
   for (int k = 0; k < vc_size; k++) {
      cs_complex_t bL[nc_size];

      // sparse lower triangular solve for next column of BL^T
      if (k >= 0 && k < p_y_link_lo_t->n) {
         if (-1 != cs_spsolve(u_nc_t_local,p_y_link_lo_t,k,xi_l,bL,NULL,1))
         { 
           // update local output matrix if result is greater than tolerance
           for (int i = 0; i < nc_size; i++) {
              if(abs(bL[i]) > tollerance) cs_entry(p_BL_sp_local[this_thread],k,i,bL[i]);
            }
         } // end if spsolve success
      } // end if k < ylink->n

      cs_complex_t bU[nc_size];
      // sparse lower triangular solve for next row of BU^T
      if (k >= 0 && k < p_y_link_up->n) {
         if (-1 != cs_spsolve(L_nc_local,p_y_link_up,k,xi_u,bU,NULL,1));
      } // end if j < y->n
      for (int i = 0; i < nc_size; i++) {
         if(abs(bU[i]) > tollerance) cs_entry(p_BU_sp_local[this_thread],i,k,bU[i]);
      }
   } // end for k

#pragma omp barrier
// mat mul


   BL_local_c = cs_compress(p_BL_sp_local[this_thread]);
   cs_spfree(p_BL_sp_local[this_thread]);
   p_BL_sp_local[this_thread] = BL_local_c;


   BU_local_c = cs_compress(p_BU_sp_local[this_thread]);
   cs_spfree(p_BU_sp_local[this_thread]);
   p_BU_sp_local[this_thread] = BU_local_c;

   p_y_eq_local[this_thread] = cs_spalloc(vc_size,vc_size,2*sizeof(cs_complex_t),1,1);
   cs_entry(p_y_eq_local[this_thread],0,0,cs_zero);
   p_y_eq_local_c[this_thread] = cs_compress(p_y_eq_local[this_thread]);


#pragma omp barrier

   int nxt_idx;

   for (int k = 0; k < num_threads; k++)
   {
      nxt_idx = (this_thread+k) % (num_threads);
      cs *local_y_eq_update;
      cs *local_y_eq_sum;

      local_y_eq_update =  cs_multiply(p_BL_sp_local[nxt_idx],p_BU_sp_local[this_thread]);

      local_y_eq_sum = cs_add(p_y_eq_local_c[this_thread],local_y_eq_update,1.0,1.0);
      cs_spfree(p_y_eq_local_c[this_thread]);
      cs_spfree(local_y_eq_update);
      p_y_eq_local_c[this_thread] = local_y_eq_sum;
     //( cs_spfree(local_y_eq_sum);
   }


#pragma omp barrier
#pragma omp single
{
   for (int k = 0; k < num_threads; k++)
   {
      cs *tmp_Y_eq = Y_eq;
//      sum_of_Y_eq = cs_add(Y_eq,p_y_eq_local_c[k],1.0,1.0);
      Y_eq = cs_add(tmp_Y_eq,p_y_eq_local_c[k],1.0,1.0);
//      Y_eq = sum_of_Y_eq;
      cs_spfree(tmp_Y_eq);
   }
}

   cs_spfree(p_BU_sp_local[this_thread]);
   cs_spfree(p_BL_sp_local[this_thread]);

   cs_spfree(u_nc_t_local);
   cs_spfree(L_nc_local);

   cs_spfree(p_y_eq_local[this_thread]);
   cs_spfree(p_y_eq_local_c[this_thread]);

} // end parallel


   cs *cs_gtc, *cs_eye, *cs_eye_trip, *gtc, *cs_z_thevenin; //, *Y_eq

   cs* Y_eq_tmp = Y_eq;
   Y_eq = cs_add(y_vc,Y_eq_tmp,1.0*cs_one,-1.0*cs_one);
   cs_spfree(Y_eq_tmp);

   cs_z_thevenin = cs_extract_diags(Y_eq);

   cs_eye_trip = cs_spalloc(cs_z_thevenin->nzmax,cs_z_thevenin->nzmax,cs_z_thevenin->nzmax*sizeof(cs_complex_t),1,1);

   for (CS_INT i = 0; i < cs_z_thevenin->nzmax; i++) cs_entry(cs_eye_trip,i,i,cs_one);
   cs_eye = cs_compress(cs_eye_trip);
   cs_spfree(cs_eye_trip);

#pragma omp parallel for
   for (int i = 0; i < cs_z_thevenin->nzmax; i++) 
   {
      cs_z_thevenin->x[i] = cs_one / cs_z_thevenin->x[i];
      zth_vc_real[i] = cs_z_thevenin->x[i].real();
      zth_vc_imag[i] = cs_z_thevenin->x[i].imag();
   }

   gtc = cs_multiply(cs_z_thevenin,Y_eq);
   cs_gtc = cs_add(cs_eye,gtc,1.0*cs_one,-1.0*cs_one);

   cs_dropzeros(cs_gtc);

#pragma omp parallel for
   for (int j = 0; j < cs_gtc->m; j++ ) {
      for (CS_INT p = cs_gtc->p[j]; p < cs_gtc->p[j + 1]; p++) {
         gtc_vc_real[cs_gtc->i[p]*vc_size+j]= cs_gtc->x[p].real();
         gtc_vc_imag[cs_gtc->i[p]*vc_size+j]= cs_gtc->x[p].imag();
      }
   }


timer.stop();
times = timer.elapsed();
d_wall_time = times.wall;
d_cpu_time = (times.user + times.system);
LOG(INFO) << "Time spent on inner matrix product:          " << std::setw(8) << d_wall_time * 1e-9 << " s, cpu time: " << d_cpu_time * 1e-9 << std::endl;

LOG(INFO) << " # nonzeros in GTC: " << cs_gtc->nzmax << std::endl;
LOG(INFO) << "number of processors used: " << num_threads << std::endl;



   // free memory]
   cs_sfree(y_nc_symbolic_ordering);
   cs_nfree(y_nc_lu_factorization);
   cs_spfree(L_nc_t);
   cs_free(row_permute);
   cs_spfree(p_y_link_lo);
   cs_spfree(p_y_link_lo_t);

   cs_spfree(p_y_link_up);
   cs_spfree(cs_gtc);
   cs_spfree(cs_z_thevenin);
   cs_spfree(Y_eq);
   cs_spfree(cs_eye);
   cs_spfree(gtc);

   return 1;
}


cs_cl * cs_load_matrix(const std::string s_MatFileName) {
// -Reads formatted text file named s_MatFileName (list of indices and real and complex values)
// -Loads the file into a compressed column sparse matrix structure of real and complex values.
// -The matrix is returend by pointer.

   std::string buffer;
   int i, j;
   int i_num_rows = 0;
   int i_num_cols = 0;
   int i_num_elem = 0;
   double dummy1, dummy2;
   double f_reY, f_imY;

 // read indices used for matrix allocation
   std::ifstream in_init(s_MatFileName);
   assert(in_init.is_open());   
   while (std::getline(in_init, buffer)) {
      std::stringstream ss(buffer);
      if (ss >> i >> j >> dummy1 >> dummy2 ) {
         if (i > i_num_rows) i_num_rows = i;
         if (j > i_num_cols) i_num_cols = j;
         ++i_num_elem;
         }
      };    
   in_init.close();

// allocate and populate sparse complex matrix
   std::ifstream in_readData(s_MatFileName);
   assert(in_readData.is_open());   
   cs *A_triplet, *A;

   A_triplet = cs_spalloc(i_num_rows, i_num_cols, i_num_elem, 1, 1);
   if (A_triplet==NULL) {
      std::cout << "ERROR: reading "<< s_MatFileName<<" memory allocation failed." << std::endl;
      return NULL;
}
   cs_complex_t y_ij;
//   std::complex<double> y_ij;
   while (std::getline(in_readData, buffer)) {
      std::stringstream ss(buffer);
      if (ss >> i >> j >> f_reY >> f_imY ) {
         y_ij = cs_complex_t(f_reY,f_imY);
//         y_ij = f_reY+f_imY*1i;
         cs_entry(A_triplet, i, j, y_ij);
      };    
   };   
   in_readData.close();

   A = cs_compress(A_triplet);
   
   cs_dropzeros(A);
   cs_spfree(A_triplet);

   return A;
};


  // Extract single element from compressed column matrix:
  CS_INT cs_element(cs *A, CS_INT i, CS_INT j, cs_complex_t *c) {
    CS_INT  m, n, *Ap, *Ai;
    cs_complex_t *Ax;
    cs_complex_t val;
    if (!CS_CSC(A)) return 0;
    m = A->m; n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
    if ((i < 0) || (i >= m)) return 0;
    if ((j < 0) || (j >= n)) return 0;
    *c = cs_complex_t(0, 0);
    for (CS_INT p = Ap[j]; p < Ap[j + 1]; p++) {
      if (Ai[p] == i) {
	*c += Ax[p];
	break;
      }
    }
    return 1;
  }


  // Extract a submatrix of a compressed column matrix:
cs * cs_submatrix(cs *A, CS_INT row_a, CS_INT col_a, CS_INT row_b, CS_INT col_b) 
{
    CS_INT  m, n, *Ap, *Ai;
    cs_complex_t *Ax;
    m = A->m; n = A->n; Ap = A->p; Ai = A->i; Ax = A->x;
    cs *c_triplet, *C;
    c_triplet = cs_spalloc(m, n, 5*sizeof (cs_complex_t),1,1);

    if (!CS_CSC(A)) return NULL;
    if ((row_a < 0) || (row_a > row_b) || (row_b >= m)) return NULL;
    if ((col_a < 0) || (col_a > col_b) || (col_b >= n)) return NULL;

    for (CS_INT j = col_a; j < col_b + 1; j++)
    {
       for (CS_INT p = Ap[j]; p < Ap[j + 1]; p++) 
       {
          if ((Ai[p] >= row_a) && (Ai[p] <= row_b + 1)) 

          {
             cs_entry(c_triplet,Ai[p],j,Ax[p]);
           }
        }
     }

     C = cs_compress(c_triplet);
     cs_spfree(c_triplet);
     return C;
}

// extract diagonals of sparse complex matrix and return it in a diagonal matrix.
cs * cs_extract_diags(cs *A) {

   if (!CS_CSC(A)) return NULL;

   CS_INT n, m;  
   cs_complex_t D_entry;
   m = A->m; n = A->n;

   cs *D_triplet, *D;
   D_triplet = cs_spalloc(CS_MIN(n,m),CS_MIN(n,m),CS_MIN(n,m)*sizeof(cs_complex_t),1,1);
   for (int i = 0; i < CS_MIN(n,m); i++) {
      cs_element(A,i,i,&D_entry);
      cs_entry(D_triplet, i, i, D_entry);
   }
   D = cs_compress(D_triplet);

   cs_spfree(D_triplet);
   return D;
}


