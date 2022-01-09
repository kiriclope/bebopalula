#ifndef __VOLTUTILS__ 
#define __VOLTUTILS__ 

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

void init_ksi_init(double mean, double var) { 
    
  for(i=0; i<n_pop; i++) 
    for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) 
      ksi_init[j] = mean + sqrt(var) * white_noise(ksi_gen) ; 
} 

void init_lif_globals() { 
  
  volt = new double [n_neurons]() ; // instantaneous individual membrane voltage 
  t_spike = new double [n_neurons]() ; // time of spike emission 

  cout << "Scaling J and ext_inputs ... " << endl ; 
  for(i=0;i<n_pop;i++) { // postsynaptic 
    ext_inputs_scaled[i] = GAIN * ext_inputs[i] * sqrt_Ka[0] * (Vth-Vr) * m0 ; 
    for(j=0;j<n_pop;j++) { // presynaptic 
      if(IF_SYN_DYN) 
	/* J_scaled[j+i*n_pop] = GAIN * J[j+i*n_pop] / TAU_SYN[j+i*n_pop] * (Vth-Vr) / sqrt_Ka[j] * sqrt_Ka[0] / sqrt_Ka[j] ; */
	J_scaled[j+i*n_pop] = GAIN * J[j+i*n_pop] / TAU_SYN[j+i*n_pop] * (Vth-Vr) / sqrt_Ka[j] ; 
      else 
  	J_scaled[j+i*n_pop] = GAIN * J[j+i*n_pop] * (Vth-Vr) / sqrt_Ka[j] * sqrt_Ka[0] / sqrt_Ka[j] ;
      
      cout << "pre " << j << " post " << i << " Jij " << J[j+i*n_pop] << " Jij_scaled " << J_scaled[j+i*n_pop] << endl ; 
      
    } 
  }
  
  cout << "Done" << endl ;
  
  for(i=0; i<n_neurons; i++) 
    ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ;

  double counter=0 ;
  
  if(IF_LOW_RANK) {
    cout << "Scaling low_rank synapses... " << endl ;
    if(RANK==1)
      for(i=0; i<n_per_pop[0]; i++) 
	for(j=0; j<n_per_pop[0]; j++) {
	  ksi_scaled[i+j*n_per_pop[0]] = GAIN * cut_LR( J[0] / sqrt_K + KAPPA * ksi[i]*ksi[j] / K ) * (Vth-Vr) / TAU_SYN[0] ; 
	  if(ksi_scaled[i+j*n_per_pop[0]]==0)
	    counter++ ;
	  if(i<3 && j<3)
	    cout << ksi_scaled[i+j*n_per_pop[0]] << " " ; 	  
	}
    
    cout << endl ;

    cout << "zeros: " << counter << endl ;
    
    if(RANK==2) {
      for(i=0; i<n_per_pop[0]; i++) 
	for(j=0; j<n_per_pop[0]; j++) {
	  ksi_scaled[i+j*n_per_pop[0]] = GAIN * cut_LR( J[0] / sqrt_K + KAPPA * ksi[i]*ksi[j] / K +  KAPPA_1 * ksi_1[i]*ksi_1[j] / K )
	    * (Vth-Vr) / TAU_SYN[0] ; 
	} 
    } 
    
    cout << "Done" << endl ;     
  }
  
  duration = DURATION + TIME_STEADY + TIME_WINDOW ; 
  time_rec = TIME_REC + TIME_STEADY + TIME_WINDOW ; 
  time_rec_spikes = TIME_REC_SPIKES + TIME_STEADY ; 
  time_steady = TIME_STEADY ; 
  
} 

void delete_lif_globals() { 
  delete [] volt ; 
  delete [] t_spike ; 
}

void integrate_mem_volt() {
  volt[i_neuron] = one_minus_dt_over_tau_mem[pre_pop] * volt[i_neuron] 
    + DT * ( net_inputs[i_neuron] + Vl / TAU_MEM[pre_pop] ) ; // V(t+dt) 
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
    cout << " t " << t_time << " ms |";
        
    if(HYST_J_EE!=0)
      cout << " J_EE " << J[0] ;
    if(HYST_M0!=0 || IF_LOOP_M0) 
      cout << " m0 " << m0 ; 
    
    cout << " rates:" << fixed << setprecision(2) ; 
    for(i=0;i<n_pop;i++) 
      cout << " " << mean_rates[i]*1000./TIME_WINDOW/(double)n_per_pop[i] ; 

    if(IF_SPEC) {
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
  }
  
}

void update_postsyn_currents_LR() { 

  j=idx_post[i_neuron] ;
  post_pop = which_pop[id_post[j]] ;

  if(post_pop==0)
    if(IF_STP)
      while(post_pop==0 && j<idx_post[i_neuron] + n_post[i_neuron]) {
	inputs[pre_pop][id_post[j]] += A_u_x_stp[i_neuron] * ksi_scaled[i_neuron+id_post[j]*n_per_pop[0]] ;
	j++ ;
	post_pop = which_pop[id_post[j]] ;
      }
    else
      while(post_pop==0 && j<idx_post[i_neuron] + n_post[i_neuron]) {
	inputs[pre_pop][id_post[j]] += ksi_scaled[i_neuron+id_post[j]*n_per_pop[0]] ; 
	j++ ; 
	post_pop = which_pop[id_post[j]] ; 
      } 
  else
    for(j=idx_post[i_neuron]; j<idx_post[i_neuron] + n_post[i_neuron]; j++) {
      post_pop = which_pop[id_post[j]] ; 
      inputs[pre_pop][id_post[j]] += J_scaled[pre_pop + post_pop * n_pop] ; 
    }
  
}

void update_postsyn_currents() { 
  
  for(j=idx_post[i_neuron]; j<idx_post[i_neuron] + n_post[i_neuron]; j++) { 
    post_pop = which_pop[id_post[j]] ; 
    
    if(IF_STP && stp_synapse[pre_pop + post_pop * n_pop]) 
      /* inputs[pre_pop][id_post[j]] += (1.0 - A_u_x_stp[i_neuron] * KAPPA/sqrt_K*cos(theta[i_neuron]-theta[id_post[j]]) )  */
      /* 	* J_scaled[pre_pop + post_pop * n_pop] ;  */
      inputs[pre_pop][id_post[j]] += A_u_x_stp[i_neuron] * J_scaled[pre_pop + post_pop * n_pop] ; 
    else 
      inputs[pre_pop][id_post[j]] += J_scaled[pre_pop + post_pop * n_pop] ; 
  } 
}

void update_net_inputs() {
  
  // cout << "reseting net inputs" << endl ; 
  for(i=0;i<n_neurons;i++) 
    net_inputs[i] = ff_inputs[i] ; 
    
  if(IF_RK2) 
    for(i=0;i<n_neurons;i++) 
      net_inputs_RK2[i] = ff_inputs[i] ; 
  
  // updating net inputs 
  for(i=0;i<n_pop;i++)  // presynaptic pop 
    for(j=0;j<n_neurons;j++) { // postsynaptic neuron
      
      post_pop = which_pop[j] ; 
      
      if(IF_SYN_DYN) // exponential decay 
	inputs[i][j] *= EXP_DT_TAU_SYN[i+post_pop*n_pop] ; 
      
      net_inputs[j] += inputs[i][j] ; 
      filter_inputs[i][j] += inputs[i][j] ; 
      
      if(!IF_SYN_DYN) // delta synapses 
	inputs[i][j] = 0.0 ; 
    } 
}

void initial_conditions() { 
    
  double mean_I[2] = { -unif(rand_gen) *  (Vth-Vl) + Vth , -unif(rand_gen) *  (Vth-Vl) + Vth } ;
  double sigma_I[2] = { unif(rand_gen) *  (Vth-Vl), unif(rand_gen) *  (Vth-Vl) } ;

  double mean, sigma ;
  
  for(i=0;i<n_pop;i++) {
    mean = unif(rand_gen) * abs(Vth-Vr) ;
    sigma = unif(rand_gen) * abs(Vth-Vr) ; // divide by 4 so that the distribution is between 0 and 1 at best
    
    for(i_neuron=0; i_neuron<n_neurons; i_neuron++) {
      if(i==0)
  	inputs[i][i_neuron] = max( mean + sqrt(sigma) * white_noise(rand_gen), 0.0 ) ;
      else
  	inputs[i][i_neuron] = min( -mean + sqrt(sigma) * white_noise(rand_gen), 0.0 ) ;
    }
  }
  
  if(IF_STP) { 
    for(i=0;i<n_neurons;i++) { 
      x_stp[i] = unif(rand_gen) ; 
      u_stp[i] = unif(rand_gen) ; 
    } 
  } 
  
  double mean_V[2] = { -unif(rand_gen) * (Vth-Vl) + Vth , -unif(rand_gen) * (Vth-Vl) + Vth } ;
  double sigma_V[2] = { unif(rand_gen) * (Vth-Vl), unif(rand_gen) * (Vth-Vl) } ;
  
  for(i_neuron=0; i_neuron<n_neurons; i_neuron++) {
    
    pre_pop = which_pop[i_neuron] ;
    volt[i_neuron] = mean_V[pre_pop] + sqrt(sigma_V[pre_pop]) * white_noise(rand_gen) ;
    
    if(i_neuron<=5 || i_neuron>=n_neurons-5)
      cout << volt[i_neuron] << " " ;
    
    if(volt[i_neuron]>=Vth) { // if spike
      
      ISI = 0 ;
      t_spike[i_neuron] = t_time ;
      volt[i_neuron] = Vr ;

      if(IF_LOW_RANK)
	update_postsyn_currents_LR() ;
      else
	update_postsyn_currents() ; 
      
      if(IF_STP)
      	update_stp_variables_lif(ISI) ;
      
    } //endif spike
    
  }
  
  cout << endl ;
  
  update_net_inputs() ;
  
  double phi_ini = unif(rand_gen) * 2.0 * M_PI ;
  /* double kappa_ini = unif(rand_gen) ; */
  double kappa_ini_0 = unif(rand_gen) ;
  double kappa_ini = kappa_ini_0 ;  
  /* double sigma_ini = unif(rand_gen) * pow(K, 0.25) ;  */
  
  cout << "kappa_ini " << kappa_ini << " phi ini " << phi_ini << endl ;
  
  /* init_ksi_init( unif(rand_gen), 4.0 * unif(rand_gen) ) ;  */
  
  for(i=0;i<n_per_pop[0]; i++) {
    ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + unif(rand_gen) ) ;
    /* ff_inputs[i] = ext_inputs_scaled[0] * ( 1.0 + unif(rand_gen) * pow(K, 0.25) / sqrt_K ) ;  */
    /* ff_inputs[i] = ext_inputs_scaled[0] * 1.0 + ( kappa_ini / sqrt_K  */
    /* 						  * cos( theta[i] + phi_ini) ) ; 	 */
  } 
  
  for(t_time=0.; t_time<=TIME_INI; t_time+=DT) {
    
    /* for(i=0;i<n_per_pop[0]; i++) { */
    /*   ff_inputs[i] = ext_inputs_scaled[0] * 1.0 + ( kappa_ini / sqrt_K  */
    /*  						    * cos( theta[i] + phi_ini) ) ; 	 */
    /* }  */
    /* ff_inputs[i] = ext_inputs_scaled[0] * (1.0 + kappa_ini * ksi_init[i] / sqrt_K ) ;  */ 
    
    /* if(t_time>=TIME_INI/2.0)  */
    /*   if(kappa_ini>0) */
    /* 	kappa_ini -= kappa_ini_0 * 2.0 * (DT/TIME_INI) ; */
    
    /* cout << "time " << fixed << setprecision(3) << t_time ; */
    /* cout << " kappa_ini " << fixed << setprecision(3) << kappa_ini << "\r" ;  */
    
    for(i_neuron=0; i_neuron<n_neurons; i_neuron++) { 

      pre_pop = which_pop[i_neuron] ; 
      
      integrate_mem_volt() ;
      
      if(volt[i_neuron]>=Vth) { // if spike
	
      	ISI = t_time - t_spike[i_neuron] ;
      	t_spike[i_neuron] = t_time ;
      	volt[i_neuron] = Vr ;
	
      	if(IF_STP)
      	  update_stp_variables_lif(ISI) ;
	
	if(IF_LOW_RANK)
	  update_postsyn_currents_LR() ;
	else
	  update_postsyn_currents() ; 
	
      } //endif spike
      
    } //endfor neurons
    
    update_net_inputs() ;
    
  } //endfor time
  
  for(i=0;i<n_per_pop[0]; i++)
    ff_inputs[i] = ext_inputs_scaled[which_pop[i]] ;
  
  for(i_neuron=0; i_neuron<n_neurons; i_neuron++)
    t_spike[i_neuron] -= TIME_INI ;
  
  cout << endl ;
  
  /* cout << "kappa_ini " << kappa_ini << " phi ini " << phi_ini << endl ; */
  
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
    
    for(i_neuron=0; i_neuron<n_neurons; i_neuron++) {
      
      pre_pop = which_pop[i_neuron] ; 
      integrate_mem_volt() ; 
      
      if(volt[i_neuron]>=Vth) { // if spike 
	
      	ISI = t_time - t_spike[i_neuron] ;
      	t_spike[i_neuron] = t_time ; 
      	volt[i_neuron] = Vr ; 
	
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
      	  update_stp_variables_lif(ISI) ; 
	
	if(IF_LOW_RANK)
	  update_postsyn_currents_LR() ;
	else
	  update_postsyn_currents() ; 
	
      } //endif spike 
      
    } //endfor neurons 
    
    if(IF_STIM) 
      tasks_inputs() ; 
    if(IF_TRACK) 
      track_input() ; 
    
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
    
    if(t_time >= time_steady) 
      t_window += DT ; 
    
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
