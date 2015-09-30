// Test GTC algorithm
#include<iostream>
#include "compute_gtc.hpp"


// MAIN
int main( int argc, char* argv[] ){

// initialize logging    
   google::InitGoogleLogging(argv[0]);

   bool no_output = false;

   std::string load_case = "case14";
   std::cout << "Case: " << load_case << std::endl;

   cs *y_nc;
   cs * y_link_lo;
   cs * y_vc;
   cs * zth_ref;
   cs * GTC_ref;
   cs * y_link_up;
   LOG(INFO) << "Case: " << load_case << std::endl;

   y_nc = cs_load_matrix(         "test_cases/" + load_case + "_ync.txt");
   y_link_up = cs_load_matrix(  "test_cases/" + load_case + "_ylink.txt");
   y_link_lo = cs_load_matrix("test_cases/" + load_case + "_ylink_t.txt");
   y_vc = cs_load_matrix(         "test_cases/" + load_case + "_yvc.txt");

   LOG(INFO) << " # buses: " << y_nc->n+y_vc->n << std::endl;
   LOG(INFO) << " # vc nodes: " << y_vc->n << std::endl;

   int vc_size = y_vc->n;

   std::vector<double> zth_vc_real(vc_size,0.0);
   std::vector<double> zth_vc_imag(vc_size,0.0);

   std::vector<double> gtc_vc_real(1,0.0);
   std::vector<double> gtc_vc_imag(1,0.0);


   gtc_vc_real.resize(vc_size*vc_size,0.0);
   gtc_vc_imag.resize(vc_size*vc_size,0.0);

   boost::timer::cpu_timer timer;
   boost::timer::cpu_times times;
   double d_wall_time;
   double d_cpu_time;
   timer.start();


   if(1==thevenin_voltage_source(y_nc, y_link_up, y_link_lo,y_vc, zth_vc_real, zth_vc_imag, gtc_vc_real, gtc_vc_imag)){

      // show output
   if (!no_output){
      zth_ref = cs_load_matrix(     "test_cases/" + load_case + "_zths.txt");
      GTC_ref = cs_load_matrix(     "test_cases/" + load_case + "_gtcs.txt");


      if (!(gtc_vc_real.size() < vc_size*vc_size)) {
            cs_complex_t gtc_kj_ref, gtc_kj;
            for (int j=0;j<vc_size;j++){
               for (int k=0;k<vc_size;k++){
                  cs_element(GTC_ref, j, k, &gtc_kj_ref); 
                  gtc_kj = cs_complex_t( gtc_vc_real[j*vc_size+k], gtc_vc_imag[j*vc_size+k] );
                  double tol = abs(gtc_kj_ref- gtc_kj );
                  bool ok = (tol < 5e-11);
                  if ((abs(gtc_kj) > tol)){
                    std::cout << " GTC "<< std::setw(2) << j << "," << std::setw(2) << k << ": " << std::setw(30) << gtc_kj_ref << " | "
                    << std::setw(30) << gtc_kj  
                    << " | " << std::setw(25)<< tol << " "<< ok << " " << std::endl;
                 }
               }
            }
      }
      std::cout << std::endl;
   
      cs_complex_t  zth_ref_k, zth;
      for (int k = 0;k < vc_size; k++) {
         zth = cs_complex_t(zth_vc_real[k], zth_vc_imag[k]);
         cs_element(zth_ref, 0, k, &zth_ref_k);
         std::cout << "Zth " << std::setw(3) << k << " of " << vc_size-1 << ": " << std::setw(25) << zth_ref_k << " | "
         << std::setw(25) <<  zth  
         << " | " << std::setw(12)<<(abs(zth_ref_k- zth )/abs(zth_ref_k)) << " " << std::endl;
      }

      cs_spfree(zth_ref);
      cs_spfree(GTC_ref);

   }//end if not no_output

   } else { std::cout << "problem!!!!" << std::endl;}


   timer.stop();
   times = timer.elapsed();
   d_wall_time = times.wall;
   d_cpu_time = times.system+times.user;
   LOG(INFO) << "Time spent on main " << load_case << ": " << std::setw(8) << d_wall_time * 1e-9 << " s, cpu time: " << d_cpu_time * 1e-9 << std::endl;



   cs_spfree(y_nc);
   cs_spfree(y_link_up);
   cs_spfree(y_link_lo);
   cs_spfree(y_vc);

   google::ShutdownGoogleLogging();
   google::ShutDownCommandLineFlags();


   return 0;
}
