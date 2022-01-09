#ifndef __TASKSUTILS__ 
#define __TASKSUTILS__ 

void track_input() {
  uniform_real_distribution<double> unif(0.0, 1.0) ; 
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
  
void tasks_inputs() { 
  uniform_real_distribution<double> unif(0.0, 1.0) ; 
  
  double kappa_ext = KAPPA_EXT / sqrt_K ; 
  double kappa_dist = KAPPA_DIST / sqrt_K ; 
  
  if(IF_DPA || IF_DUAL) { 
    // Sample stimulus ON 
    if(t_time-TIME_STEADY > T_SAMPLE_ON && t_time-TIME_STEADY <= T_SAMPLE_ON + 2.0 * DT) { 
      for(i=0;i<n_per_pop[0]; i++) {
	if(IF_SPEC) 
	  ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  kappa_ext * cos( theta[i] + PHI_EXT ) ) ; 
	if(IF_LOW_RANK) 
	  ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 +  kappa_ext * ksi[i] ) ;	
      } 
    } 
    // Sample stimulus OFF
    if(t_time-TIME_STEADY > T_SAMPLE_OFF && t_time-TIME_STEADY <= T_SAMPLE_OFF + 2.0*DT) 
      for(i=0;i<n_per_pop[0]; i++) 
	ff_inputs[i] = ext_inputs_scaled[0] ; 
  }
  
  // Distractor ON 
  if(IF_DUAL || IF_DRT) { 
    if(t_time-TIME_STEADY > T_DIST_ON && t_time-TIME_STEADY <= T_DIST_ON + 2.0*DT) {
      for(i=0;i<n_per_pop[0]; i++) { 
	/* ff_inputs[i] = ( 1.0 + KAPPA_EXT / sqrt(K) ) * ext_inputs_scaled[which_pop[i]] ; */ 

	if(IF_SPEC) {
	  if(RANK==1) 
	    ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + kappa_dist * unif(rand_gen) ) ; 
	
	  if(RANK==2) 
	    ff_inputs[idx_perm[i]] = ext_inputs_scaled[0] * ( 1.0 + kappa_dist * cos( theta[i] + PHI_DIST) ) ; 
	}
	
	if(IF_LOW_RANK) {
	  if(RANK==1)
	    ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + kappa_dist * unif(rand_gen) ) ; 
	  if(RANK==2)
	    ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + kappa_dist* ksi_1[i] ) ; 
	} 	
      }
    }
    // Distractor OFF
    if(t_time-TIME_STEADY > T_DIST_OFF && t_time-TIME_STEADY <= T_DIST_OFF + 2.0*DT) {
      for(i=0;i<n_per_pop[0]; i++) 
	ff_inputs[i] = ext_inputs_scaled[0] ; 
    } 
  } 
  
  if(IF_DPA || IF_DUAL) {
    // Test stimulus ON 
    if(t_time-TIME_STEADY > T_TEST_ON && t_time-TIME_STEADY <= T_TEST_ON + 2.0*DT) { 
      for(i=0;i<n_per_pop[0]; i++) 
	ff_inputs[i] = ( 1.0 + kappa_ext * unif(rand_gen) ) * ext_inputs_scaled[0] ;       
    }
    
    // Test stimulus OFF
    if(t_time-TIME_STEADY > T_TEST_OFF && t_time-TIME_STEADY <= T_TEST_OFF + 2.0*DT) 
      for(i=0;i<n_per_pop[0]; i++) 
	ff_inputs[i] = ext_inputs_scaled[0] ; 
  }
  
}

#endif 
