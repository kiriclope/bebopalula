#include <iostream> 
#include <iomanip> 
#include <sstream>
#include <fstream>
#include <cmath> 
#include <random>
#include <ctime> 
#include <sys/types.h> 
#include <sys/stat.h> 
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <gsl/gsl_permutation.h> 

#include "globals.h" 
#include "mean_field.h" 
#include "net_utils.h" 
#include "con_utils.h" 


int main(int argc , char** argv) { 
  
  get_args(argc , argv) ;
  
  IF_STRUCTURE = IF_RING || IF_SPEC || IF_LOW_RANK || IF_GAUSS ; 
  
  get_param() ; 
  init_globals() ; 

  create_con_dir() ; 
  
  cout << "coucou" << endl ; 
  init_ksi() ;
  
  write_to_file(ksi_path, "ksi", ksi , n_neurons) ; 
  // write_to_file(con_path, "ksi_scaled", ksi_scaled , n_per_pop[0]) ; 
    
  if(RANK==2) {
    init_ksi_1() ; 
    write_to_file(ksi_path, "ksi_1", ksi_1 , n_neurons) ; 
    // write_to_file(con_path, "ksi_1_scaled", ksi_1_scaled , n_per_pop[0]) ; 
  }
  
  delete_globals() ; 
  
}
