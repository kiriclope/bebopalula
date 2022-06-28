#ifndef __VOLTUTILS__
#define __VOLTUTILS__


double print_PSP(int i_pop, int j_pop) {
  double tp, vpsp ;
  tp = log(TAU_MEM[i_pop]/TAU_SYN[j_pop + i_pop * n_pop])/(1.0/TAU_MEM[i_pop]-1.0/TAU_SYN[j_pop + i_pop * n_pop]) ; 
  vpsp = exp(-tp/TAU_SYN[j_pop + i_pop * n_pop]) ; 
  vpsp = 20.0 * vpsp / Vth ;  
  return vpsp ;
}
  

void open_lif_files() {
  string str_volt = path + "/mem_volt.dat" ;
  file_volt.open(str_volt.c_str(), ios::out | ios::ate); 

  string str_spike_times = path + "/spike_times.dat" ; 
  file_spike_times.open(str_spike_times.c_str(), ios::out | ios::ate);   
}


void close_lif_files() {
  file_volt.close() ;
  file_spike_times.close() ;
}


void save_volt() {

  file_volt << t_time - TIME_STEADY ;
  
  for(i=0;i<n_neurons;i++) { 
    if(t_spike[i_neuron] == t_time) 
      file_volt << " " << Vpeak ; 
    else 
      file_volt << " " << volt[i] ; 
  } 
  file_volt << endl ; 
  
}


void init_ksi_init(double mean, double var) { 
  
  for(i=0; i<n_pop; i++) 
    for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) 
      ksi_init[j] = mean + sqrt(var) * white_noise(ksi_gen) ; 
} 


void scale_ksi() { // ONLY USE IF CONN PROB NOT CHANGED IN con_mat.h
  cout << "Scaling low_rank synapses... " << endl ;

  double ksi_dum ;
  int dum ;
  unsigned long counter=0 ;
  
  for(i=0; i<n_pop; i++) // post
    for(j=0; j<n_pop; j++) { // pre
      dum = j+i*n_pop ; 
      if(IS_STRUCT_SYN[dum]) 
	for(k=cum_n_per_pop[i]; k<cum_n_per_pop[i+1]; k++) 	
	  for(l=cum_n_per_pop[j]; l<cum_n_per_pop[j+1]; l++) {
	    
	    ksi_dum = J[dum] / sqrt_Ka[j] ; 
	    ksi_dum += KAPPA * ksi[k]*ksi[l] / K ; 
	    
	    if(RANK==2) 
	      ksi_dum += KAPPA_1 * ksi_1[k]*ksi_1[l] / K ;
	    
	    ksi_scaled[k+l*n_neurons] = GAIN * cut_LR( ksi_dum ) ; 
	    if(IF_SYN_DYN) ksi_scaled[k+l*n_neurons] /= TAU_SYN[dum] ; 
	    if(IF_RESCALE) ksi_scaled[k+l*n_neurons] *= (Vth-Vr) ; 
	    
	    if(ksi_scaled[k+l*n_neurons]==0) counter++ ; 	    
	    if(k<3 && l<3) cout << ksi_scaled[k+l*n_per_pop[0]] << " " ; 
	  }

      cout << endl ; 
      cout << "zeros: " << counter / n_per_pop[i]  << endl ; 
      
    } 
  cout << "Done" << endl ;       
}

void scale_J() {
  cout << "1/sqrt(K) " << 1.0/sqrt_K << endl ; 
  cout << "Scaling J ... " << endl ; 
  for(i=0;i<n_pop;i++) { // postsynaptic 
    for(j=0;j<n_pop;j++) { // presynaptic 
      J_scaled[j+i*n_pop] = GAIN * J[j+i*n_pop] / sqrt_Ka[j] ; 
      if(IF_SYN_DYN) { 
	J_scaled[j+i*n_pop] /= TAU_SYN[j+i*n_pop] ; 
	// might be the right thing given my implementation 
	J_scaled[j+i*n_pop] *= EXP_DT_TAU_SYN[j+i*n_pop] ;
      } 
      if(IF_MATO_K) J_scaled[j+i*n_pop] *= sqrt_Ka[0] / sqrt_Ka[j] ; 
      if(IF_RESCALE) J_scaled[j+i*n_pop] *= (Vth-Vr) ; 
      cout << "pre " << j << " post " << i << " Jij " << J[j+i*n_pop] ; 
      cout << " Jij_scaled " << J_scaled[j+i*n_pop] ; 
      cout << " PSP " << J_scaled[j+i*n_pop] * print_PSP(i,j) << endl ; 
    } 
  } 
  /* if(IF_NMDA) { */ 
  /*   /\* J_scaled[0] *= R_NMDA[0] / (1.0 + R_NMDA[0]) ;  *\/ */
  /*   /\* J_scaled[2] *= R_NMDA[1] / (1.0 + R_NMDA[1]) ;  *\/ */
  /* }   */
}

void scale_J_nmda() { 
  cout << "Scaling J_nmda ... " << endl ; 
  for(i=0;i<n_pop;i++) { // postsynaptic 
    for(j=0;j<n_pop;j++) { // presynaptic 
      J_nmda[j+i*n_pop] = GAIN * J[j+i*n_pop] / sqrt_Ka[j] ; 
      J_nmda[j+i*n_pop] /= TAU_NMDA[j+i*n_pop] ; 
      J_nmda[j+i*n_pop] *= EXP_DT_TAU_NMDA[j+i*n_pop] ;
      if(IF_MATO_K) J_nmda[j+i*n_pop] *= sqrt_Ka[0] / sqrt_Ka[j] ; 
      if(IF_RESCALE) J_nmda[j+i*n_pop] *= (Vth-Vr) ; 
      cout << "pre " << j << " post " << i << " Jij " << J[j+i*n_pop] ; 
      cout << " Jij_nmda " << J_nmda[j+i*n_pop] ; 
      cout << " PSP " << J_nmda[j+i*n_pop] * print_PSP(i,j) << endl ; 
    } 
  } 
  J_nmda[0] *= 0.92 ; // 1.0 / (1.0 + R_NMDA[0]) ; 
  J_nmda[2] *= 7.4 / J[2] ; //1.0 / (1.0 + R_NMDA[1]) ; 
  J_nmda[3] = 0.0 ; // 1.0 / (1.0 + R_NMDA[1]) ; 
  J_nmda[4] = 0.0 ; // 1.0 / (1.0 + R_NMDA[1]) ; 
} 

void scale_ext_inputs() {
  cout << "sqrt(K) " << sqrt_K << endl ; 
  cout << "Scaling ext_inputs ... " << endl ;  
  for(i=0;i<n_pop;i++) { // postsynaptic
    ext_inputs_scaled[i] = ext_inputs[i] * sqrt_Ka[0] * m0 ; 
    if(IF_RESCALE) ext_inputs_scaled[i] *= (Vth-Vr) ; 
    cout << "raw " << ext_inputs[i] << " scaled " << ext_inputs_scaled[i] << " " ; 
  } 
  cout << endl ; 
  
  for(i=0; i<n_neurons; i++) 
    ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ; 
  if(IF_TUNED_FF) 
    for(i=0; i<n_neurons; i++) 
      ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  KAPPA_EXT/sqrt_Ka[0] * cos(theta[i]-PHI_EXT)) ; 
}

void init_lif_globals() { 

  duration = DURATION + TIME_STEADY + TIME_WINDOW ; 
  time_rec = TIME_REC + TIME_STEADY + TIME_WINDOW ; 
  time_rec_spikes = TIME_REC_SPIKES + TIME_STEADY ; 
  time_steady = TIME_STEADY ; 
  
  volt = new double [n_neurons]() ; // instantaneous individual membrane voltage 
  t_spike = new double [n_neurons]() ; // time of spike emission 

  scale_ext_inputs() ; 
  scale_J() ; 
  if(IF_NMDA) scale_J_nmda() ; 
  
  /* if(IF_LOW_RANK) scale_ksi() ; */ 
  
}

void delete_lif_globals() { 
  delete [] volt ; 
  delete [] t_spike ; 
}

void integrate_mem_volt() {

  if(IF_RK2) { 
    /* RK1 = -(volt[i_neuron]-Vl) / TAU_MEM[pre_pop] + net_inputs[i_neuron] ; */
    /* RK2 = -(volt[i_neuron]-Vl + DT*RK1) / TAU_MEM[pre_pop] + net_inputs_RK2[i_neuron] ; */
    /* volt[i_neuron] = volt[i_neuron] + DT/2.0 * ( RK1 + RK2 ) ;  */
    volt[i_neuron] /= one_plus_dt_over_tau_mem[pre_pop] ;
    volt[i_neuron] += net_inputs[i_neuron] * dt_over_dt_tau_mem[pre_pop] ; 
  } 
  else { 
    /* volt[i_neuron] *= one_minus_dt_over_tau_mem[pre_pop] ; */ 
    volt[i_neuron] *= EXP_DT_TAU_MEM[pre_pop] ; 
    volt[i_neuron] += dt_over_tau_mem[pre_pop] * net_inputs[i_neuron] ; // V(t+dt)
    /* volt[i_neuron] += DT * ( net_inputs[i_neuron] + Vl/TAU_MEM[pre_pop] ) ; // V(t+dt) */ 
  } 
} 

void print_rates() {
  
  if( t_window<TIME_WINDOW ) { 
    cout << int(percentage*100.0) << "% " ; 
    cout << "\r" ; 
    cout.flush() ; 
  }
  
  if( t_window>=TIME_WINDOW ) { 
    
    cout << int(percentage*100.0) << "% " ; 
    
    // mean rates 
    cout << " t " << t_time-time_steady << " ms |";
    
    if(HYST_J_EE!=0)
      cout << " J_EE " << J[0] ;
    if(HYST_M0!=0 || IF_LOOP_M0) 
      cout << fixed << setprecision(3) << " m0 " << m0 ; 
    
    cout << " rates:" << fixed << setprecision(3) ;
    if(IF_SINGLE_NEURON)
      for(i=0;i<n_pop;i++) 
	cout << " " << mean_rates[i]*1000./TIME_WINDOW ; 
    else
      for(i=0;i<n_pop;i++) 
	cout << " " << mean_rates[i]*1000./TIME_WINDOW/(double)n_per_pop[i] ; 
    
    if(IF_SPEC || IF_RING || IF_GAUSS) { 
      get_m1_phase() ;
      cout << " m1 " ; 
      for(i=0;i<n_pop;i++)
	cout << m1[i] << " " ;
      
      cout << "phase " ; 
      for(i=0;i<n_pop;i++)
	cout << phase[i] << " " ;
    }
    
    if(IF_LOW_RANK) {
      cout << " | overlaps:" ; 
      for(i=0;i<n_pop;i++) 
	/* cout << " " << overlaps[i] * DT * IS_STRUCT_SYN[i] ;  */
	cout << " " << overlaps[i] * 1000. / TIME_WINDOW / (double) n_per_pop[i] ; 
    }
    
    cout << "\r" ;
    cout.flush() ;
    /* cout << endl ;  */
  }
  
}


void update_postsyn_currents_LR() { 
  
  j=idx_post[i_neuron] ;
  post_pop = which_pop[id_post[j]] ;
  
  if(pre_pop==0) { // only E to E/I 
    
    if(IF_STP) {
      while(post_pop==0) { // only EtoE
	inputs[pre_pop][id_post[j]] += A_u_x_stp[i_neuron] * ksi_scaled[i_neuron+id_post[j]*n_per_pop[0]] ; 
	j++ ; 
	post_pop = which_pop[id_post[j]] ; 
      }
      
      while(j<idx_post[i_neuron] + n_post[i_neuron]) { // EtoI 
	inputs[pre_pop][id_post[j]] += J_scaled[pre_pop + post_pop * n_pop] ; 
	j++ ; 
      } 
      
    }// end STP 
    else { // no STP       
      while(post_pop==0) { // only EtoE
	inputs[pre_pop][id_post[j]] += ksi_scaled[i_neuron+id_post[j]*n_per_pop[0]] ; 
	j++ ; 
	post_pop = which_pop[id_post[j]] ; 
      }
      
      while(j<idx_post[i_neuron] + n_post[i_neuron]) { // EtoI 
	inputs[pre_pop][id_post[j]] += J_scaled[pre_pop + post_pop * n_pop] ; 
	j++ ; 
      } 
      
    } // end no STP
  }
  else 
    for(j=idx_post[i_neuron]; j<idx_post[i_neuron] + n_post[i_neuron]; j++) { 
      post_pop = which_pop[id_post[j]] ; 
      inputs[pre_pop][id_post[j]] += J_scaled[pre_pop + post_pop * n_pop] ; 
    } 
} 

void update_postsyn_currents() {
  int counter_E=0, counter_I=0 ;
  /* cout << "pre_pop "<< pre_pop << " i_neuron " << i_neuron << endl ; */
  
  double J_dum ; 
  for(j=idx_post[i_neuron]; j<idx_post[i_neuron] + n_post[i_neuron]; j++) { 
    post_pop = which_pop[id_post[j]] ;
    
    /* cout << "pre_pop "<< pre_pop << " i_neuron " << i_neuron ; */
    /* cout << " post_pop "<< post_pop << " j_neuron " << id_post[j] ;  */
    /* cout << " J " << J_scaled[pre_pop + post_pop * n_pop] ; */
    /* cout << "\r" ; */
    /* cout.flush() ; */
    
    /* if(post_pop==0) */
    /*   counter_E++ ; */
    /* if(post_pop==1) */
    /*   counter_I++ ; */
        
    J_dum = J_scaled[pre_pop + post_pop * n_pop] ; 
    if(IF_STP && stp_synapse[pre_pop + post_pop * n_pop]) 
      J_dum *= A_u_x_stp[i_neuron] ; 
    if(IF_RK2) 
      J_dum *= exp(-(ISI[i_neuron])/TAU_SYN[pre_pop + post_pop * n_pop]) ; 
    inputs[pre_pop][id_post[j]] += J_dum ; 
  }
  
  /* cout << "n post " << counter_E << " " << counter_I << endl ;  */
  
} 

void update_postsyn_currents_nmda() {
  double J_dum ; 
  if(pre_pop==0)
    for(j=idx_post[i_neuron]; j<idx_post[i_neuron] + n_post[i_neuron]; j++) { 
      post_pop = which_pop[id_post[j]] ; 
      J_dum = J_nmda[pre_pop + post_pop * n_pop] ; 
      if(IF_STP && stp_synapse[pre_pop + post_pop * n_pop]) 
	J_dum *= A_u_x_stp[i_neuron] ; 
      if(IF_RK2) 
	J_dum *= exp(-(ISI[i_neuron])/TAU_NMDA[pre_pop+post_pop*n_pop]) ; 
      inputs_nmda[pre_pop][id_post[j]] += J_dum ; 
    } 
} 

void update_net_inputs() {   
  // cout << "reseting net inputs" << endl ;
  if(SIGMA_FF>0)
    for(i=0;i<n_neurons;i++) {
      net_inputs[i] = ff_inputs[i] + sqrt(SIGMA_FF) * white_noise(rand_gen) ; 
      /* release_stp() ;  */
    } 
  else 
    for(i=0;i<n_neurons;i++) 
      net_inputs[i] = ff_inputs[i] ; 
  
  if(IF_RK2) 
    for(i=0;i<n_neurons;i++) 
      net_inputs_RK2[i] = ff_inputs[i] ; 
  
  // updating net inputs 
  for(i=0;i<n_pop;i++)  // presynaptic pop 
    for(j=0;j<n_neurons;j++) { // postsynaptic neuron 
      post_pop = which_pop[j] ; 
      net_inputs[j] += inputs[i][j] ; // must be before decay 
      if(IF_SYN_DYN) inputs[i][j] *= EXP_DT_TAU_SYN[i+post_pop*n_pop] ;
      
      if(IF_RK2) net_inputs_RK2[j] += inputs[i][j] * EXP_DT_TAU_SYN[i+post_pop*n_pop] ; 
      filter_inputs[i][j] += inputs[i][j] ; 
      
      if(IF_NMDA) { 
	net_inputs[j] += inputs_nmda[i][j] ; // must be before decay 
	inputs_nmda[i][j] *= EXP_DT_TAU_NMDA[i+post_pop*n_pop] ; 
	
	if(IF_RK2) net_inputs_RK2[j] += inputs_nmda[i][j] * EXP_DT_TAU_NMDA[i+post_pop*n_pop] ; 
	filter_inputs[i][j] += inputs_nmda[i][j] ; 
      }      
      if(!IF_SYN_DYN) inputs[i][j] = 0.0 ; // delta synapses 
    } 
} 

void update_single_neuron_inputs() { 
  net_inputs[i_neuron] = ff_inputs[i_neuron] ; 
}

void initial_conditions() { 
  
  /* double mean_I[2] = { -unif(rand_gen) *  (Vth-Vr) + Vth , -unif(rand_gen) *  (Vth-Vr) + Vth } ; */
  /* double sigma_I[2] = { unif(rand_gen) *  (Vth-Vr), unif(rand_gen) *  (Vth-Vr) } ; */

  /* double mean, sigma ; */
  
  /* for(i=0;i<n_pop;i++) { */
  /*   mean = unif(rand_gen) * abs(Vth-Vr) ; */
  /*   sigma = unif(rand_gen) * abs(Vth-Vr) ; // divide by 4 so that the distribution is between 0 and 1 at best */
    
  /*   for(i_neuron=0; i_neuron<n_neurons; i_neuron++) { */
  /*     if(i==0) */
  /* 	inputs[i][i_neuron] = max( mean + sqrt(sigma) * white_noise(rand_gen), 0.0 ) ; */
  /*     else */
  /* 	inputs[i][i_neuron] = min( -mean + sqrt(sigma) * white_noise(rand_gen), 0.0 ) ; */
  /*   } */
  /* } */
  
  if(IF_STP) { 
    for(i=0;i<n_neurons;i++) { 
      /* x_stp[i] = unif(rand_gen) ;  */
      /* u_stp[i] = unif(rand_gen) ; */

      x_stp[i] = 1.0 ;
      u_stp[i] = USE ;
    } 
  } 
  
  /* double mean_V[2] = { -unif(rand_gen) * (Vth-Vr) + Vth , -unif(rand_gen) * (Vth-Vr) + Vth } ;  */
  /* double sigma_V[2] = { unif(rand_gen) * (Vth-Vr), unif(rand_gen) * (Vth-Vr) } ; */
  
  for(i_neuron=0; i_neuron<n_neurons; i_neuron++) {  
    pre_pop = which_pop[i_neuron] ;
    /* volt[i_neuron] = mean_V[pre_pop] + sqrt(sigma_V[pre_pop]) * white_noise(rand_gen) ; */
    volt[i_neuron] = unif(rand_gen) * (Vth-Vr) + Vr ;
    /* volt[i_neuron] = mean_V[pre_pop] + sqrt(sigma_V[pre_pop]) * white_noise(rand_gen) ; */
  }    
  /*   if(i_neuron<=5 || i_neuron>=n_neurons-5) */
  /*     cout << volt[i_neuron] << " " ; */
    
  /*   if(volt[i_neuron]>=Vth) { // if spike */
      
  /*     t_spike[i_neuron] = t_time ; */
  /*     volt[i_neuron] = Vr ;  */
  /*     ISI[i_neuron] = t_time - t_spike[i_neuron] ;  */
      
  /*     /\* if(IF_LOW_RANK) *\/ */
  /*     /\* 	update_postsyn_currents_LR() ; *\/ */
  /*     /\* else *\/ */
      
  /*     update_postsyn_currents() ; */
  /*     if(IF_NMDA) */
  /*   	update_postsyn_currents_nmda() ; */
      
  /*     if(IF_STP) */
  /*     	update_stp_variables_lif() ; */

  /*     ISI[i_neuron] = 0.0 ; */
      
  /*   } //endif spike */
  /*   else */
  /*     release_stp() ; */
  /* } */
  
  /* cout << endl ; */
  
  /* update_net_inputs() ; */
  
  /* double kappa_ini = unif(rand_gen) ; */
  /* double sigma_ini = 4.0 * unif(rand_gen) ;  */
  /* double phi_ini = unif(rand_gen) * 2.0 * M_PI ;  */
  
  /* cout << "kappa_ini " << kappa_ini << " phi ini " << phi_ini * 180.0 / M_PI << endl ;  */
  
  /* if(!IF_STRUCTURE)  */
  /*   for(i=0;i<n_per_pop[0]; i++)  */
  /*     ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + unif(rand_gen) ) ;  */

  /* if(IF_LOW_RANK) { */
  /*   init_ksi_init( unif(rand_gen), 4.0 * unif(rand_gen) ) ;  */
  /*   for(i=0;i<n_per_pop[0]; i++) */
  /*     ff_inputs[i] = ext_inputs_scaled[0] * (1.0 + kappa_ini * ksi_init[i] / sqrt_K ) ; */
  /* } */
  
  /* if(IF_SPEC)  */
  /*   for(i=0;i<n_per_pop[0]; i++)  */
  /*     ff_inputs[i] = ext_inputs_scaled[0] * 1.0 + ( kappa_ini / sqrt_K  */
  /* 						    * cos( theta[i] + phi_ini) ) ; */
  
  /* if(IF_RING)  */
  /*   ff_inputs[i] = ext_inputs_scaled[0] * 1.0 + ( kappa_ini * cos( theta[i] + phi_ini) ) ;  */
  
  /* for(t_time=0.; t_time<TIME_INI; t_time+=DT) { */
    
  /*   /\* if(t_time>=TIME_INI/2.0)  */
  /*   /\*   if(kappa_ini>0) *\/  */
  /*   /\* 	kappa_ini -= kappa_ini_0 * 2.0 * (DT/TIME_INI) ; *\/  */
    
  /*   /\* cout << "time " << fixed << setprecision(3) << t_time ; *\/ */
  /*   /\* cout << " kappa_ini " << fixed << setprecision(3) << kappa_ini << "\r" ;  *\/ */

  /*   for(i=0;i<n_per_pop[0]; i++)  */
  /*     ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + sqrt(sigma_ini) * white_noise(rand_gen) ) ;  */
    
  /*   for(i_neuron=0; i_neuron<n_neurons; i_neuron++) {  */
      
  /*     pre_pop = which_pop[i_neuron] ;  */
  /*     vold = volt[i] ;  */
      
  /*     integrate_mem_volt() ;  */
      
  /*     if(volt[i_neuron]>=Vth) { // if spike		 */
  /* 	if(IF_RK2) { */
  /* 	  /\* t_spike[i_neuron] = t_time + DT * (Vth-vold) / (volt[i_neuron]-vold) ;  *\/ */
  /* 	  t_spike[i_neuron] = t_time + DT * (Vth-volt[i_neuron]) / (volt[i_neuron]-vold) ;  */
  /* 	  ISI[i_neuron] = t_time - t_spike[i_neuron] ;  */
  /* 	  volt[i_neuron] = (volt[i_neuron]-Vth) */
  /* 	    * ( 1.0 + dt_over_tau_mem[pre_pop] * (vold-Vl) / (volt[i_neuron]-vold) ) + Vr ;  */
  /* 	}  */
  /* 	else {  */
  /* 	  ISI[i_neuron] = t_time - t_spike[i_neuron] ; */
  /* 	  t_spike[i_neuron] = t_time ;	 */
  /* 	  volt[i_neuron] = Vr ; 	   */
  /* 	} */
	
  /* 	if(IF_STP)  */
  /* 	  update_stp_variables_lif() ;  */
	
  /* 	/\* if(IF_LOW_RANK)  *\/ */
  /* 	/\*   update_postsyn_currents_LR() ;  *\/ */
  /* 	/\* else  *\/	 */
  /* 	update_postsyn_currents() ;  */
  /* 	if(IF_NMDA)  */
  /* 	  update_postsyn_currents_nmda() ;  */
		
  /*     } //endif spike  */
  /*     /\* else *\/ */
  /*     /\* 	release_stp() ; *\/ */
      
  /*   } //endfor neurons */
    
  /*   update_net_inputs() ;  */
    
  /* } //endfor time */
  
  /* if(IF_TUNED_FF) */
  /*   for(i=0; i<n_neurons; i++)  */
  /*     ff_inputs[i] = ext_inputs_scaled[which_pop[i]] * ( 1.0 +  KAPPA_EXT/sqrt_Ka[0] * cos( theta[i] - PHI_EXT ) ) ;  */
  /* else  */
  /*   for(i=0;i<n_per_pop[0]; i++)  */
  /*     ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ;  */
  
  /* for(i_neuron=0; i_neuron<n_neurons; i_neuron++)  */
  /*   t_spike[i_neuron] -= TIME_INI ;  */
  
  /* cout << endl ;  */
  
  /* cout << "kappa_ini " << kappa_ini << " phi ini " << phi_ini << endl ; */
  
}

void run_single_neuron() {

  open_files() ; 
  
  init_lif_globals() ; 
  open_lif_files() ; 
  
  if(IF_STP) { 
    open_stp_files() ; 
    init_stp_globals() ; 
  }

  i_neuron = 0 ;

  double sigma_ext = 10.0 * unif(rand_gen) ; 
    
  for(t_time=0.; t_time<=duration; t_time+=DT) { 
    
    ff_inputs[i_neuron] = 0.25 * ext_inputs_scaled[i_neuron] * ( 1.0 + sqrt(sigma_ext) * white_noise(rand_gen) ) ; 
    
    percentage = t_time/duration ; 
    
    pre_pop = which_pop[i_neuron] ; 
    
    integrate_mem_volt() ; 
    
    if(volt[i_neuron]>=Vth) { // if spike 
      
      ISI[i_neuron] = t_time - t_spike[i_neuron] ; 
      t_spike[i_neuron] = t_time ; 
      volt[i_neuron] = Vr ; 
      
      if(t_time>=time_steady) {	  
	mean_rates[pre_pop] += 1.0 ; 
	filter_rates[i_neuron] += 1.0 ; 
	
	if(t_time<time_rec_spikes) // save spike times 
	  file_spike_times << fixed << setprecision(1) << (float) (i_neuron) 
			   << " " << (float) (t_spike[i_neuron]-time_steady) << endl ; 
      } 
      
      if(IF_STP) 
	update_stp_variables_lif() ; 
      
    } // end_spike 
    
    update_single_neuron_inputs() ;

    print_rates() ; 
    
    if(t_window>=TIME_WINDOW) {
      
      if(t_time<=time_rec)
	save_to_file() ; 
            
      if(IF_HYSTERESIS) 
	hyst_update() ; 
      
      if(IF_STP) 
	save_xy_to_file() ; 
      
      t_window=0. ; 
    } 
    
    if(t_time >= time_steady) {
      t_window += DT ; 
      if(t_time<=time_rec && IF_SAVE_VOLT) 
	save_volt() ; 
    }
    
  } //endfor time
  
  cout << endl ; 
  
  delete_globals() ; 
  close_files() ; 
  
  delete_lif_globals() ; 
  close_lif_files() ; 
  
  if(IF_STP) {
    delete_stp_globals() ; 
    close_stp_files() ;
  }  
  
}

void run_sim_lif() { 

  open_files() ;
  
  init_lif_globals() ; 
  open_lif_files() ; 

  if(IF_STP) {
    open_stp_files() ;
    init_stp_globals() ;
  }
  
  if(IF_HYSTERESIS) 
    hyst_init() ; 
  
  initial_conditions() ;
  
  print_sim_info() ; 
  
  for(t_time=0.; t_time<=duration; t_time+=DT) { 
    
    percentage = t_time/duration ; 

    if(IF_STIM) 
      tasks_inputs() ; 
    if(IF_TRACK) 
      track_input() ; 
    
    for(i_neuron=0; i_neuron<n_neurons; i_neuron++) {
      
      pre_pop = which_pop[i_neuron] ; 
      integrate_mem_volt() ; 
      
      if(volt[i_neuron]>=Vth) { // if spike 
	
	if(IF_RK2) { 
	  /* t_spike[i_neuron] = t_time + DT * (Vth-vold) / (volt[i_neuron]-vold) ;  this is what meunier wrote */ 
	  t_spike[i_neuron] = t_time + DT * (Vth-volt[i_neuron]) / (volt[i_neuron]-vold) ; // this is what hansel is using
	  ISI[i_neuron] = t_time - t_spike[i_neuron] ; 
	  volt[i_neuron] = (volt[i_neuron]-Vth) * (one_plus_dt_over_tau_mem[pre_pop] * (vold-Vl) /(volt[i_neuron]-vold) ) + Vl ; 
	}
	else { 
	  ISI[i_neuron] = t_time - t_spike[i_neuron] ; 
	  t_spike[i_neuron] = t_time ; 
	  volt[i_neuron] = Vr ; 
	}
	
      	if(t_time>=time_steady) { 
      	  mean_rates[pre_pop] += 1.0 ; 
      	  filter_rates[i_neuron] += 1.0 ; 
	  
	  if(IF_LOW_RANK) 
	    overlaps[pre_pop] += ksi[i_neuron] ; 
	  
      	  if(t_time<time_rec_spikes) // save spike times 
      	    file_spike_times << fixed << setprecision(1) << (float) (i_neuron) 
      	  		     << " " << (float) (t_spike[i_neuron]-time_steady) << endl ; 
      	}
	
	if(IF_STP) 
      	  update_stp_variables_lif() ; 
	
	/* if(IF_LOW_RANK) */
	/*   update_postsyn_currents_LR() ; */
	/* else  */
	update_postsyn_currents() ; 
	if(IF_NMDA) 
	  update_postsyn_currents_nmda() ; 
	
      } //endif spike 
      /* else */
      /* 	release_stp() ; */
      
    } //endfor neurons 
        
    update_net_inputs() ; 
    
    print_rates() ; 
    
    if(t_window>=TIME_WINDOW) {
      
      if(t_time<=time_rec)
	save_to_file() ; 
            
      if(IF_HYSTERESIS) 
	hyst_update() ; 
      
      if(IF_STP) 
	save_xy_to_file() ; 
      
      t_window=0. ; 
    } 
    
    if(t_time >= time_steady) {
      t_window += DT ; 
      if(t_time<=time_rec && IF_SAVE_VOLT) 
	save_volt() ; 
    }
    
  } //endfor time 
  cout << endl ; 
  
  delete_globals() ; 
  close_files() ; 
  
  delete_lif_globals() ; 
  close_lif_files() ; 
  
  if(IF_STP) {
    delete_stp_globals() ; 
    close_stp_files() ;
  }
}

#endif
