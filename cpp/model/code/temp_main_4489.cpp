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

#include "temp_globals_4489.h" 
#include "mean_field.h" 
#include "net_utils.h" 
#include "con_utils.h"
#include "rmvnorm.h"
#include "mat_utils.h" 
#include "stp_utils.h"
#include "tasks_utils.h"
#include "hyst_utils.h" 
#include "lif_utils.h" 
#include "binary_utils.h"


int main(int argc , char** argv) { 
  
  get_args(argc , argv) ;
  
  IF_STRUCTURE = IF_RING || IF_SPEC || IF_LOW_RANK || IF_GAUSS ;
  int IF_SIM= IF_BIN || IF_LIF || IF_RATE ;  
  IF_STIM = IF_DPA || IF_DUAL || IF_DRT || IF_STEP || IF_CHRISTOS ; 
  
  get_param() ; 
  init_globals() ; 
  
  if(IF_GEN_KSI)
    gen_ksi() ;
  
  if(IF_SINGLE_NEURON==0)
    if(IF_GEN_CON) 
      gen_con_sparse_vec() ; 
    else 
      get_con_sparse_vec() ; 
  
  if(IF_SIM) {
    
    create_dir() ; 
    mean_field_rates() ; 
    
    if(IF_LIF)
      if(IF_SINGLE_NEURON)
	run_single_neuron() ;
      else
	run_sim_lif() ;
    
    if(IF_BIN) 
      run_sim_bin() ; 
  }
}
