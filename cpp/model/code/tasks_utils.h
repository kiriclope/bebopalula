#ifndef __TASKSUTILS__ 
#define __TASKSUTILS__ 

void track_input() {
  double dum1=0, dum2=0, dum3=0, dum4=0 ; 
    
  if(t_time-TIME_STEADY > 2000.0 && dum1==0) {
    dum1=1 ;
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  KAPPA_EXT / sqrt(K) * cos( 2.0 * (double) i * M_PI / (double) n_per_pop[which_pop[i]] ) ) ; 
  }
  
  if(t_time-TIME_STEADY > 6000.0 && dum2==0) { 
    dum2=1 ; 
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  KAPPA_EXT / sqrt(K) * cos( 2.0 * (double) i * M_PI / (double) n_per_pop[which_pop[i]] + M_PI/2.0) ) ; 
  }

  if(t_time-TIME_STEADY > 10000.0 && dum3==0) { 
    dum3=1 ; 
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  KAPPA_EXT / sqrt(K) * cos( 2.0 * (double) i * M_PI / (double) n_per_pop[which_pop[i]] + M_PI) ) ; 
  }
    
  if(t_time-TIME_STEADY > 14000.0 && dum4==0) {
    dum4=1 ;
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  KAPPA_EXT / sqrt(K) * cos( 2.0 * (double) i * M_PI / (double) n_per_pop[which_pop[i]] + 3*M_PI/2.0 ) ) ; 
  }
  
}

void christos_task() { 
  // CUE 
  if(t_time-TIME_STEADY >= T_CUE_ON && t_time-TIME_STEADY < T_CUE_OFF  && !SWITCH_ON) {
    for(i=0;i<n_neurons; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]]
	+ sqrt_Ka[0] * A_CUE[which_pop[i]] 
	* ( 1.0 + EPS_CUE[which_pop[i]] * cos( theta[i] - 2.0 * PHI_CUE * M_PI) ) ; 
    SWITCH_ON = 1 ; 
  }
  if(t_time-TIME_STEADY >= T_CUE_OFF && SWITCH_ON) { 
    for(i=0;i<n_neurons; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
    SWITCH_ON = 0 ; 
  } 
  
  // ERASE BUMP 
  if(t_time-TIME_STEADY >= T_ERASE_ON && t_time-TIME_STEADY < T_ERASE_OFF  && !SWITCH_OFF) {
    if(IF_DIST)
      for(i=0;i<n_neurons; i++) 
	ff_inputs[i] = ext_inputs_scaled[which_pop[i]]
	  + A_ERASE[which_pop[i]] 
	  * (1.0 + EPS_ERASE[which_pop[i]] * cos( theta[i] - 2.0 * PHI_ERASE * M_PI)) ;
    else
      for(i=0;i<n_neurons; i++) 
	ff_inputs[i] = ext_inputs_scaled[which_pop[i]]
	  + sqrt_Ka[0] * A_ERASE[which_pop[i]] 
	  * (1.0 + EPS_ERASE[which_pop[i]] * cos( theta[i] - 2.0 * PHI_ERASE * M_PI)) ;
    
    SWITCH_OFF = 1 ; 
  }
  if(t_time-TIME_STEADY >= T_ERASE_OFF && SWITCH_OFF) { 
    for(i=0;i<n_neurons; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
    SWITCH_OFF = 0 ; 
  }   
}

void DPA_task() {
  double kappa_ext = KAPPA_EXT / sqrt_Ka[0] ; 
  double kappa_dist = KAPPA_DIST / sqrt_Ka[0] ;  
  double kappa_test = KAPPA_TEST / sqrt_Ka[0] ; 
 
  // Sample stimulus ON 
  if(t_time-TIME_STEADY >= T_SAMPLE_ON && t_time-TIME_STEADY < T_SAMPLE_OFF  && !SWITCH_ON) {     
    if(IF_SPEC) 
      for(i=0;i<n_per_pop[0]; i++) 
	ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  kappa_ext * cos( theta[i] -  2.0	* PHI_EXT  * M_PI ) ) ; 
    if(IF_LOW_RANK) {
      for(i=0;i<n_per_pop[0]; i++) 
	ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  kappa_ext * sample[i] ) ;
      for(i=n_per_pop[0];i<n_per_pop[1]; i++)
	ff_inputs[i] = ext_inputs_scaled[1] * ( 1.0 +  kappa_ext * unif(rand_gen) ) ; 
    }
    
    if(IF_RING)
      for(i=0;i<n_per_pop[0]; i++) 
	ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  KAPPA_EXT * cos( theta[i] - 2.0 * PHI_EXT  * M_PI) ) ;
    
    SWITCH_ON = 1 ; 
  }

  // Sample stimulus OFF
  if(t_time-TIME_STEADY >= T_SAMPLE_OFF && SWITCH_ON) { 
    for(i=0;i<n_neurons; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
    SWITCH_ON = 0 ; 
  }
  
  // Test stimulus ON 
  if(t_time-TIME_STEADY > T_TEST_ON && t_time-TIME_STEADY <= T_TEST_ON + 2.0*DT) { 
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  kappa_test * distractor[i] ) ;
    
    for(i=n_per_pop[0];i<n_per_pop[1]; i++)
      ff_inputs[i] = ext_inputs_scaled[1] * ( 1.0 +  kappa_ext * unif(rand_gen) ) ; 
  } 
  
  // Test stimulus OFF 
  if(t_time-TIME_STEADY > T_TEST_OFF && t_time-TIME_STEADY <= T_TEST_OFF + 2.0*DT) 
    for(i=0;i<n_neurons; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
}

void DRT_task() {
  double kappa_dist = KAPPA_DIST / sqrt_K ; 
  double kappa_cue = KAPPA_CUE / sqrt_K ; 
  
  // Distractor ON 
  if(t_time-TIME_STEADY > T_DIST_ON && t_time-TIME_STEADY <= T_DIST_ON + 2.0*DT) {
      if(IF_SPEC) { 
	for(i=0;i<n_per_pop[0]; i++) { 
	  if(RANK==1) 
	    ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + kappa_dist * unif(rand_gen) ) ; 
	  if(RANK==2) 
	    ff_inputs[idx_perm[i]] = ext_inputs_scaled[0] * ( 1.0 + kappa_dist * cos( theta[i] + PHI_DIST) ) ; 
	}
      }
      
      if(IF_LOW_RANK) 
	for(i=0;i<n_per_pop[0]; i++) 
	  ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  kappa_dist * distractor[i] ) ; 
  }
  
  // Distractor OFF
  if(t_time-TIME_STEADY > T_DIST_OFF && t_time-TIME_STEADY <= T_DIST_OFF + 2.0*DT) {
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[0] ;    
  }

  // Cue/Reward ON
  if(t_time-TIME_STEADY > T_RWD_ON && t_time-TIME_STEADY <= T_RWD_ON + 2.0*DT) {
      if(IF_SPEC) { 
	for(i=0;i<n_per_pop[0]; i++) { 
	  if(RANK==1) 
	    ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + kappa_cue * unif(rand_gen) ) ; 
	  if(RANK==2) 
	    ff_inputs[idx_perm[i]] = ext_inputs_scaled[0] * ( 1.0 + kappa_dist * cos( theta[i] + PHI_DIST) ) ; 
	}
      }
      
      if(IF_LOW_RANK) 
	for(i=0;i<n_per_pop[0]; i++) 
	  ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  kappa_cue * distractor[i] ) ; 
  }
  
  // Rwd OFF
  if(t_time-TIME_STEADY > T_RWD_OFF && t_time-TIME_STEADY <= T_RWD_OFF + 2.0*DT) {
    for(i=0;i<n_per_pop[0]; i++) 
      ff_inputs[i] = ext_inputs_scaled[0] ;    
  }
  
} 

void christos_dist() {
  
  if(t_time-TIME_STEADY >= T_STEP_ON && t_time-TIME_STEADY < T_STEP_OFF  && !SWITCH_ON) {
    for(i=0;i<n_neurons; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]]
	+ sqrt_Ka[0] * A_CUE[which_pop[i]] 
	* ( 1.0 + EPS_CUE[which_pop[i]] * cos( theta[i] - 2.0 * PHI_EXT * M_PI) ) ; 
    
    SWITCH_ON = 1 ; 
  }
  
  if(t_time-TIME_STEADY >= T_STEP_OFF && SWITCH_ON) { 
    for(i=0;i<n_neurons; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
    SWITCH_ON = 0 ; 
  } 

  if(t_time-TIME_STEADY >= T_DIST_ON && t_time-TIME_STEADY < T_DIST_OFF  && !SWITCH_ON) {
    for(i=0;i<n_neurons; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]]
	+ sqrt_Ka[0] * A_CUE[which_pop[i]] 
	* ( 1.0 + EPS_CUE[which_pop[i]] * cos( theta[i] - 2.0 * PHI_DIST * M_PI) ) ; 
    
    SWITCH_ON = 1 ; 
  }
  
  if(t_time-TIME_STEADY >= T_DIST_OFF && SWITCH_ON) { 
    for(i=0;i<n_neurons; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
    SWITCH_ON = 0 ; 
  } 
  
}

void step_input() {
  double kappa_ext = KAPPA_EXT / sqrt_K ; 
  
  // Sample stimulus ON 
  if(t_time-TIME_STEADY >= T_STEP_ON && t_time-TIME_STEADY < T_STEP_OFF  && !SWITCH_ON) { 
    for(i=0;i<n_neurons; i++) { 
      if(IF_STRUCTURE) { 
	if(IF_RING) 
	  ff_inputs[i] = ext_inputs_scaled[which_pop[i]]
	    + sqrt_Ka[0] * A_CUE[which_pop[i]] 
	    * ( 1.0 + EPS_CUE[which_pop[i]] * cos( theta[i] - 2.0 * PHI_EXT * M_PI) ) ; 
	  /* ff_inputs[i] = ext_inputs_scaled[which_pop[i]] */
	  /*   * ( 1.0 +  KAPPA_EXT * cos( theta[i] - PHI_EXT ) ) ;  */
	if(IF_SPEC && i<n_per_pop[0]) 
	  ff_inputs[i] = ext_inputs_scaled[which_pop[i]]
	    * ( 1.0 +  kappa_ext * cos( theta[i] - PHI_EXT ) ) ; 
	if(IF_LOW_RANK) 
	  ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  kappa_ext * ksi[i] ) ;
	if(IF_GAUSS) 
	  ff_inputs[i] = ext_inputs_scaled[which_pop[i]]
	    + sqrt_Ka[0] * A_CUE[which_pop[i]] 
	    * ( 1.0 + EPS_CUE[which_pop[i]] * cos( theta[i] - 2.0 * PHI_EXT * M_PI) ) ; 
      }
      else
	  ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * A_STEP[which_pop[i]] ; 
    } 
    SWITCH_ON = 1 ; 
  }
  // Sample stimulus OFF 
  if(t_time-TIME_STEADY >= T_STEP_OFF && SWITCH_ON) { 
    for(i=0;i<n_neurons; i++) 
	if(IF_TUNED_FF) 
	    ff_inputs[i] = ext_inputs_scaled[which_pop[i]]
	      * ( 1.0 +  KAPPA_EXT/sqrt_Ka[0] * cos(theta[i]- 2.0 * PHI_EXT * M_PI)) ; 
	else
	  ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
    SWITCH_ON = 0 ; 
  } 
} 
  
void tasks_inputs() { 
  if(IF_STEP)
    step_input() ;
  
  if(IF_CHRISTOS)
    christos_task() ;
    
  if(IF_DPA || IF_DUAL) 
    DPA_task() ; 
  
  if(IF_DUAL || IF_DRT) 
    DRT_task() ; 
}

#endif 
