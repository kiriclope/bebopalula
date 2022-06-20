#include "librairies.h"  
#include "globals.h" 
#include "net_utils.h" 
#include "mat_utils.h" 
#include "low_rank_utils.h" 
#include "lif_utils.h" 
#include "stp_utils.h" 
#include "mean_field.h" 

clock_t t1=clock(); 

int main(int argc , char** argv) { 
  
  unsigned long i,j ; 
  
  // default_random_engine gen  ; 
  random_device rd ; 
  // default_random_engine rand_gen( rd() ) ; 
  mt19937 rand_gen(rd()), rand_gen_xi_left(rd()),
    rand_gen_xi_right(rd()), rand_gen_shared(rd()), rand_gen_ff(rd()) ; 
  
  int pre_pop, post_pop ; 
  double t_window=0., percentage=0. ; 
    
  string dir ;  
  int n_pop ; 
  unsigned long n_neurons ; 
  double K ;
  get_args(argc , argv, dir, n_pop, n_neurons, K) ; 
  
  double *ext_inputs, *J, *tau_syn ; 
  get_param(n_pop, dir, ext_inputs, J, tau_syn) ; 
  
  if(n_pop==1) { 
    J[0] = J0 ; 
    ext_inputs[0] = I0 ; 
    tau_syn[0] = 10.0 ; 
    
    char char_dir[30] ; 
    sprintf(char_dir, "I0_%0.2f_J0_%0.2f", I0, -J0) ; 
    string str_dir = string(char_dir) ; 
    
    dir = str_dir ; 
    cout << dir ; 
  }   
    
  cout << "External Inputs : " ; 
  for(int i=0;i<n_pop;i++) 
    cout << ext_inputs[i] << " "; 
  cout << endl ; 
  
  cout << "J: " ;
  for(int i=0;i<2*n_pop;i++)
    cout << J[i] << " " ; 
  cout << endl ; 
    
  cout << "tau_syn : " ;
  for(i=0;i<n_pop;i++) // postsynaptic pop
    for(j=0;j<n_pop;j++) // presynaptic pop
      cout << tau_syn[j+i*n_pop] << " " ;
  cout << endl ;

  cout << "GAIN: " << GAIN << endl ;
  
  double *mf_rates ; 
  mf_rates = new double [n_pop]() ;
  mean_field_rates(n_pop, ext_inputs, J, mf_rates) ; 
  
  ///////////////////////////////////////////////////////////////////    
  // Random Connectivity, J
  ///////////////////////////////////////////////////////////////////    
  
  int *n_post ;
  unsigned long *id_post, *idx_post ;
  
  string con_path = "/homecentral/alexandre.mahrach/IDIBAPS/connectivity/" ;  
  con_path += to_string(n_pop) + "pop/N" + to_string(n_neurons) + "/K" + to_string((int)K) ; 

  cout << "n_neurons" << n_neurons << endl;
  
  if(n_neurons!=0) 
    get_con_mat_lif(con_path, n_neurons, n_post, id_post, idx_post) ; 
  else { 
    n_post = new int [2]() ;
    idx_post = new unsigned long [2]() ;          
    id_post = (unsigned long *) malloc( (unsigned long) 4 * sizeof(unsigned long) ) ;

    n_post[0] = 2 ; 
    n_post[1] = 0 ;

    idx_post[0] = 0 ; 
    idx_post[1] = 0 ; 
    
    id_post[0] = 0 ; 
    id_post[1] = 1 ; 
    
    id_post[2] = 0 ; 
    id_post[3] = 1 ; 
    
  }
  
  // cout << n_post[0] << " " << id_post[0] << " " << idx_post[0] << endl ; 
  
  ///////////////////////////////////////////////////////////////////    
  // Path
  ///////////////////////////////////////////////////////////////////    
  
  string path = "../" ; 
  create_dir(dir, path, n_pop, n_neurons, K) ; 
  
  if(n_neurons!=0) 
    n_neurons = n_neurons * 10000 ; 
  else 
    n_neurons = 2 ; 
  
  int* n_per_pop ; 
  n_per_pop = new int [n_pop]() ; 
  
  for(i=0;i<n_pop;i++) 
    n_per_pop[i] = n_neurons/n_pop ; 
  
  int* cum_n_per_pop ; 
  cum_n_per_pop = new int [n_pop+1]() ; 
  
  cout <<"cum_n_per_pop=" << " " ; 
  for(i=0;i<n_pop+1;i++) { 
    for(j=0;j<i;j++) 
      cum_n_per_pop[i] += n_per_pop[j] ; 
    cout << cum_n_per_pop[i] << " " ;
  }
  cout << endl ; 
  
  int *which_pop ; 
  which_pop = (int *) malloc( (unsigned long long) n_neurons * sizeof(int) ) ;
  
  for(i=0;i<n_pop;i++) 
    for(j=0; j<n_neurons; j++)
      if(j>=cum_n_per_pop[i] && j<cum_n_per_pop[i+1])
	which_pop[j] = i ; 

  
  ///////////////////////////////////////////////////////////////////    
  // Structured Connectivity, xi
  ///////////////////////////////////////////////////////////////////    

  normal_distribution<double> white_noise(0.0, 1.0) ; 
  
  double *xi, *xi_left, *xi_right ; 
  xi = new double [n_neurons]() ; // gaussian vector 

  double x_left, x_right, x_shared ; 
  
  // if(IF_LOW_RANK) 
  //   create_path_low_rank(path, n_pop) ; 
  
  string sLow = path + "/low_rank_xi.dat" ; 
  ofstream file_low_rank_xi(sLow.c_str(), ios::out | ios::ate) ;

  for(i=0;i<n_neurons;i++) {// xi xj scales as 1/N
    x_left = white_noise(rand_gen_xi_left) ;
    
    xi[i] = (MEAN_XI[which_pop[i]] + sqrt(VAR_XI[which_pop[i]]) * x_left ) / sqrt( (double) n_per_pop[which_pop[i]] ) ; 
  }
  
  cout << "low rank vector: " ;
  for(i=0;i<5;i++)
    cout << xi[i] << " " ;
  cout << endl ; 
  
  /////////////////////////////////////////////////////////////////// 
  // Scaling 
  ///////////////////////////////////////////////////////////////////
  
  for(i=0;i<n_pop;i++) { 
    ext_inputs[i] = sqrt(K) * ext_inputs[i] * (Vth-Vr) * m0 ; 
    for(j=0;j<n_pop;j++) {
      if(IF_SYN_DYN)
  	J[j+i*n_pop] = GAIN * J[j+i*n_pop] / TAU_SYN[j+i*n_pop] * (Vth-Vr) / sqrt(K) ;
      else
  	J[j+i*n_pop] = GAIN * J[j+i*n_pop] * (Vth-Vr) / sqrt(K) ;
    } 
  } 
  
  // for(i=0;i<n_pop;i++) { 
  //   ext_inputs[i] = sqrt(K) * ext_inputs[i] * m0 ; 
  //   for(j=0;j<n_pop;j++) {
  //     if(IF_SYN_DYN)
  // 	J[j+i*n_pop] = GAIN * J[j+i*n_pop] / TAU_SYN[j+i*n_pop] / sqrt(K) ;
  //     else
  // 	J[j+i*n_pop] = GAIN * J[j+i*n_pop] / sqrt(K) ; 
  //   } 
  // } 
  
  double* ext_inputs_BL ; 
  ext_inputs_BL = new double [n_pop]() ; 
  for(i=0; i<n_pop;i++)
    ext_inputs_BL[i] = ext_inputs[i] ; 

  double *ff_inputs ;
  ff_inputs = new double [n_neurons]() ;
  for(i=0; i<n_neurons;i++) 
    ff_inputs[i] = ext_inputs[which_pop[i]] ; 
  
  /////////////////////////////////////////////////////////////////// 
  // Variables 
  ///////////////////////////////////////////////////////////////////    

  double *volt ;
  volt = new double [n_neurons]() ; // instantaneous individual membrane voltage 
  
  double *t_spike ;
  t_spike = new double [n_neurons]() ; // time of spike emission

  double t_spike_old = 0 ;
    
  double *mean_rates ;
  mean_rates = new double [n_pop]() ; // population averaged rate also averaged over TIME_WINDOW
  
  double *filter_rates ;
  filter_rates = new double [n_neurons]() ; // temporal averaged over TIME_WINDOW
    
  // h^(ab)_i=h^b_i, inputs from presynaptic population b to postsynaptic neuron (i,a)
  // here i goes threw all the populations 
  double **inputs ; 
  inputs = new double *[n_pop]() ;
  for(i=0;i<n_pop;i++) // presynaptic population b
    inputs[i] = new double [n_neurons]() ; 
  
  double **filter_inputs ; 
  filter_inputs = new double *[n_pop]() ; 
  for(i=0;i<n_pop;i++) // presynaptic population b
    filter_inputs[i] = new double [n_neurons]() ; 
  
  // htot_i = h^E_i + h^I_i, net input into neuron i
  double *net_inputs ; 
  net_inputs = new double [n_neurons]() ; 
  
  double *net_inputs_RK2 ; 
  net_inputs_RK2 = new double [n_neurons]() ; 
  
  //////////////////////////////////////////////////////////////
  // STP variables
  //////////////////////////////////////////////////////////////
    
  double *u_stp ; // availability variable 
  u_stp = new double [n_neurons]() ; 
  
  double *x_stp ; // resource variable 
  x_stp = new double [n_neurons]() ; 
  
  //////////////////////////////////////////////////////////////
  // Low Rank variables 
  //////////////////////////////////////////////////////////////
  
  double *kappa ; // overlap averaged over TIME_WINDOW
  kappa = new double [n_pop]() ; 
  
  ///////////////////////////////////////////////////////////////////    
  // Files
  ///////////////////////////////////////////////////////////////////    

  string str_volt = path + "/mem_volt.dat" ; 
  ofstream file_volt(str_volt.c_str(), ios::out | ios::ate);

  string str_spike_times = path + "/spike_times.dat" ; 
  ofstream file_spike_times(str_spike_times.c_str(), ios::out | ios::ate);
  
  string str_mean_rates = path + "/mean_rates.dat" ; 
  ofstream file_mean_rates(str_mean_rates.c_str(), ios::out | ios::ate);
  
  string str_filter_rates = path + "/filter_rates.dat" ; 
  ofstream file_filter_rates(str_filter_rates.c_str(), ios::out | ios::ate);

  string str_inputs = path + "/inputs.dat" ; 
  ofstream file_inputs(str_inputs.c_str(), ios::out | ios::ate);

  // STP
  string str_u_stp = path + "/u_stp.dat" ; 
  ofstream file_u_stp(str_u_stp.c_str(), ios::out | ios::ate);
  
  string str_x_stp = path + "/x_stp.dat" ; 
  ofstream file_x_stp(str_x_stp.c_str(), ios::out | ios::ate); 
  
  // Low rank
  string str_overlap = path + "/overlap.dat" ; 
  ofstream file_overlap(str_overlap.c_str(), ios::out | ios::ate); 
  
  ///////////////////////////////////////////////////////////////////    
  // Initial conditions
  ///////////////////////////////////////////////////////////////////
  
  normal_distribution<double> gaussian( 0, (Vth-Vl)/4.0 ) ; 
  
  for(i=0;i<n_neurons;i++) {
    volt[i] = TAU_MEM[which_pop[i]]*Vl ; 
    net_inputs[i] = ext_inputs[which_pop[i]] ; 
  } 
  
  for(i=0;i<n_pop;i++) 
    for(j=0;j<n_neurons;j++) {
      post_pop = which_pop[j] ;
      
      inputs[i][j] = (K * J[i + post_pop * n_pop] * TAU_SYN[i + post_pop * n_pop] * mf_rates[i]
		      / (double) n_per_pop[i] + gaussian(rand_gen) /TAU_MEM[i] ) ; 
      
      net_inputs[j] += inputs[i][j] ; 
    }
  
  ///////////////////////////////////////////////////////////////////    
  // Dynamics of the network : Rate model
  ///////////////////////////////////////////////////////////////////    
  
  cout << "LIF Model " ;
  if(IF_EULER)
    cout << "with EULER integration " ; 
  else 
    if(IF_RK2) 
      cout << "with RK2 and interpolation " ; 
  
  if(IF_STP) 
    cout << "with STP " ; 

  if(IF_LOW_RANK)
    cout << "with added LOW RANK structure " ;
  
  cout << endl ; 
    
  cout << "Main loop :" ; 
  cout << " duration " << DURATION ; 
  cout << " ms | DT " << DT ; 
  cout << " ms | TIME_STEADY " << TIME_STEADY ;
  cout << " ms | TIME_REC " << TIME_REC ; 
  cout << " ms | TIME_WINDOW " << TIME_WINDOW << " ms " << endl ; 

  if(IF_STP==0) {
    cout << "mean_field_rates: " ;
    for(i=0;i<n_pop;i++)
      cout << mf_rates[i]*1000. << " " ;
    cout << endl ;
  }
  
  int dum1=0, dum2=0;
  double vold, RK1, RK2 ; 
  
  for (double time=0.; time<=DURATION; time+=DT) { 
    
    percentage = time/DURATION ; 
    
    if(time>=TIME_STEADY  && time<TIME_STEADY + TIME_REC) 
      file_volt << time-TIME_STEADY ; 
    
    // cout << "updating membrane voltages" << endl ; 
    for(i=0; i<n_neurons; i++) { // presynaptic pop 
      pre_pop = which_pop[i] ; 
      vold = volt[i] ; 
      
      if(IF_STP && ( stp_synapse[pre_pop] || stp_synapse[pre_pop + n_pop] ) ) { 
	
	if(time-TIME_STEADY > 0.30*TIME_REC && time-TIME_STEADY <= 0.40*TIME_REC) { 
	  // ff_inputs[i] = ( 1.0 + KAPPA / sqrt(K) ) * ext_inputs[which_pop[i]] ; 
	  
	  // ff_inputs[i] += ext_inputs[which_pop[i]]*DT/0.1/TIME_REC 
	  //   * cos( 2.0 * (double) i * 3.142 / n_per_pop[pre_pop] ) / sqrt(K) ; 
	  
	  ff_inputs[i] = ext_inputs[which_pop[i]] * ( 2.0 +  KAPPA / sqrt(K) 
	  					      * cos( 2.0 * (double) i * 3.142 / (double) n_per_pop[pre_pop] ) ) ; 
	}
	
	if(time-TIME_STEADY > 0.60*TIME_REC && time-TIME_STEADY <= 0.70*TIME_REC ) { 
	  // ff_inputs[i] -= ext_inputs[which_pop[i]]*DT/0.1/TIME_REC 
	  // * cos( 2.0 * (double) i * 3.142 / (double) n_per_pop[pre_pop] ) / sqrt(K) ; 
	  
	  ff_inputs[i] = ext_inputs[which_pop[i]] ;
	}
      }
      
      if(IF_EULER) { // Standard Euler 
	
	volt[i] = one_minus_dt_over_tau_mem[pre_pop] * volt[i]
	  + DT * ( net_inputs[i] + Vl / TAU_MEM[pre_pop] ) ; // V(t+dt) 
	
	if(volt[i]>=Vth) { // spike
	  
	  if(IF_STP && ( stp_synapse[pre_pop] || stp_synapse[pre_pop + n_pop] ) ) 
	    update_stp_variables_lif(u_stp[i], x_stp[i], time-t_spike[i], pre_pop) ; 
	  
	  t_spike[i] = time ; 
	  volt[i] = Vr ; 
	  
	  if(time>=TIME_STEADY) { 
	    mean_rates[pre_pop] += 1.0 ; 
	    filter_rates[i] += 1.0 ;
	    
	    if(IF_LOW_RANK && low_synapse[pre_pop] ) 
	      kappa[which_pop[i]] += xi[i] ; 
	      
	    if(i<10 || (i>=n_per_pop[0]+10 && i<n_per_pop[0]+20) ) 
	      file_volt << " " << Vpeak ; 
	  } 
	  
	  // update postsynaptic currents 
	  for(j=idx_post[i]; j<idx_post[i] + n_post[i]; j++) { 
	    post_pop = which_pop[id_post[j]] ; 
	    
	    if(IF_STP && stp_synapse[pre_pop + post_pop * n_pop])
	      inputs[pre_pop][id_post[j]] += A_STP * u_stp[i] * x_stp[i] * J[pre_pop + post_pop * n_pop] ; 
	    // inputs[pre_pop][id_post[j]] += J[pre_pop + post_pop * n_pop] 
	    // 	* ( 1.0 + u_stp[i] * x_stp[i] * KAPPA / sqrt(K) * cos( 2.0 * (double) (i - id_post[j] ) * 3.142 / (double) n_per_pop[pre_pop] ) ) ; 
	    else 
	      inputs[pre_pop][id_post[j]] += J[pre_pop + post_pop * n_pop] ; 
	    
	    if(IF_LOW_RANK && low_synapse[pre_pop] ) 
	      if(IF_STP && stp_synapse[pre_pop + post_pop * n_pop]) 
		inputs[pre_pop][id_post[j]] += eps[pre_pop] * xi[i] * xi[i] * u_stp[i] * x_stp[i] ; 
	      else 
		inputs[pre_pop][id_post[j]] += eps[pre_pop] * xi[i] * xi[i] ; 
	    
	  } 
	  
	  // save spike times 
	  if(time>=TIME_STEADY && time<TIME_STEADY + TIME_REC) 
	    file_spike_times << fixed << setprecision(1) << (float) (i) 
			     << " " << (float) (t_spike[i]-TIME_STEADY) << endl ; 
	  
	} // endif spike 
	else
	  if(time>=TIME_STEADY  && time<TIME_STEADY + TIME_REC && (i<10 || (i>=n_per_pop[0]+10 && i<n_per_pop[0]+20) ) ) 
	    file_volt << " " << volt[i] ;       
	
      } // endif EULER 
      else 
	if(IF_RK2) { // RK2 with interpolation as in Hansel et al. (Neltner), 1998. 
	
	  RK1 = -(volt[i]-Vl) / TAU_MEM[pre_pop] + net_inputs[i] ; 
	  RK2 = -(volt[i]-Vl + DT*RK1) / TAU_MEM[pre_pop] + net_inputs_RK2[i] ; 
	  volt[i] = volt[i] + DT/2.0 * ( RK1 + RK2 ) ; 
	  
	  if(volt[i]>=Vth) { // spike
	    t_spike_old = t_spike[i] ; 
	    t_spike[i] = time + DT * (Vth-vold) / (volt[i]-vold) ;
	    
	    if(IF_STP && ( stp_synapse[pre_pop] || stp_synapse[pre_pop + n_pop] ) ) 
	      update_stp_variables_lif(u_stp[i], x_stp[i], t_spike[i]-t_spike_old, pre_pop) ; 
	    
	    volt[i] = (volt[i]-Vth) * (1.0+DT/TAU_MEM[pre_pop] * (vold-Vl) /(volt[i]-vold) ) + Vl ; 
	    
	    if(time>=TIME_STEADY  && time<TIME_STEADY + TIME_REC) { 
	      mean_rates[pre_pop] += 1.0 ; 
	      filter_rates[i] += 1.0 ;
	      if(i<10 || (i>=n_per_pop[0]+10 && i<n_per_pop[0]+20) ) 
		file_volt << " " << Vpeak ; 
	    }
	    
	    // update postsynaptic currents
	    for(j=idx_post[i]; j<idx_post[i] + n_post[i]; j++) { 
	      post_pop = which_pop[id_post[j]] ; 
	      
	      if(IF_STP && stp_synapse[pre_pop + post_pop * n_pop]) 
		inputs[pre_pop][id_post[j]] += u_stp[i] * x_stp[i] * J[pre_pop + post_pop * n_pop] ; 
	      else
		inputs[pre_pop][id_post[j]] += J[pre_pop + post_pop * n_pop] ; 
	      
	      inputs[pre_pop][id_post[j]] *= exp(-(time-t_spike[i])/TAU_SYN[pre_pop+post_pop*n_pop]) ; 
	    } 
	    
	    // save spike times 
	    if( time>=TIME_STEADY && time<TIME_STEADY + TIME_REC )  
	      file_spike_times << fixed << setprecision(1) << (float) (i) 
			       << " " << (float) (t_spike[i]-TIME_STEADY) << endl ; 
	    
	  } //endif spike 
	  else 
	    if(time>=TIME_STEADY && time<TIME_STEADY + TIME_REC && (i<10 || (i>=n_per_pop[0]+10 && i<n_per_pop[0]+20) ) ) 
	      file_volt << " " << volt[i] ; 
	  
	} // endif RK2 
    } // endfor n_neurons 

    if(time>=TIME_STEADY && time<TIME_STEADY + TIME_REC) 
      file_volt << endl ; 
    
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
	if(IF_SYN_DYN)
	  inputs[i][j] *= EXP_DT_TAU_SYN[i+post_pop*n_pop] ; 
	net_inputs[j] += inputs[i][j] ; 
	filter_inputs[i][j] += inputs[i][j] ; 
	if(!IF_SYN_DYN)
	  inputs[i][j] = 0.0 ; 
      }
    
    if(IF_RK2)
    for(i=0;i<n_pop;i++)  // presynaptic pop 
      for(j=0;j<n_neurons;j++) { // postsynaptic neuron 
	post_pop = which_pop[j] ; 
	net_inputs_RK2[j] += inputs[i][j] * EXP_DT_TAU_SYN[i+post_pop*n_pop] ; 
      }
    
    // Writing to files 
    if(t_window<TIME_WINDOW) { 
      cout << int(percentage*100.0) << "% " ; 
      cout << "\r" ; 
      cout.flush() ; 
    }    
    else{ 
      
      cout << int(percentage*100.0) << "% " ; 

      // mean rates
      cout << " t " << time - TIME_STEADY << " ms | ext_input " << ff_inputs[i]  << " << rates:"; 
      for(i=0;i<n_pop;i++) 
      	cout << " " << mean_rates[i]*1000./TIME_WINDOW/(double)n_per_pop[i] ; 

      if(IF_LOW_RANK) {	
	cout << "| overlaps:" ;
	for(i=0;i<n_pop;i++) 
	  cout << " " << kappa[i]*DT/TIME_WINDOW/sqrt( (double)n_per_pop[i]) ; 
      }
      
      cout << "\r" ; 
      cout.flush() ; 
      // cout << endl ; 
      
      file_mean_rates << time - TIME_STEADY ; 
      for(i=0;i<n_pop;i++) {
      	file_mean_rates << " " << mean_rates[i]*1000./TIME_WINDOW/(double)n_per_pop[i]; 
      	mean_rates[i] = 0 ; 
      } 
      file_mean_rates << endl ; 
      
      // filtered rates over tw
      file_filter_rates << time - TIME_STEADY ;
      for(i=0;i<n_neurons;i++) {
  	file_filter_rates << " " << filter_rates[i]*1000./TIME_WINDOW ; 
  	filter_rates[i] = 0 ; 
      } 
      file_filter_rates << endl ; 
      
      // filtered inputs over tw 
      file_inputs << time - TIME_STEADY ; 
      for(i=0;i<n_pop;i++) 
	for(j=0;j<n_neurons;j++) { 
	  file_inputs << " " << filter_inputs[i][j]*DT/TIME_WINDOW ; 
	  filter_inputs[i][j] = 0 ; 
	} 
      file_inputs << endl ; 

      if(IF_STP) { 
	file_u_stp << time-TIME_STEADY ; 
	for(j=0;j<n_neurons;j++) 
	  file_u_stp << " " << u_stp[j] ;  
	file_u_stp << endl ; 
	
	file_x_stp << time-TIME_STEADY ; 
	for(j=0;j<n_neurons;j++) 
	  file_x_stp << " " << x_stp[j] ;  
	file_x_stp << endl ; 
      }

      // overlap kappa
      if(IF_LOW_RANK) { 
	file_overlap << time-TIME_STEADY ; 
	for(i=0;i<n_pop;i++) { 
	  file_overlap << " " << kappa[i]*DT/TIME_WINDOW/sqrt( (double)n_per_pop[i]) ; 
	  kappa[i] = 0.0 ; 
	} 
	file_overlap << endl ; 
      } 
      
      t_window=0. ; 
      
    }//endif 

    
    // printProgress (percentage) ;
    //Defining time window 
    if(time >= TIME_STEADY) 
      t_window += DT ;
    
  } //ENDMAINLOOP
  
    ///////////////////////////////////////////////////////////////////

  delete [] idx_post ; 
  delete [] id_post ; 
  delete [] n_post ; 

  delete [] volt ; 

  delete [] mean_rates ; 
  delete [] filter_rates ; 

  delete [] inputs ; 
  delete [] filter_inputs ; 
  delete [] net_inputs ; 
  delete [] net_inputs_RK2 ; 
  
  delete [] ext_inputs ; 
  delete [] ext_inputs_BL ; 
  delete [] ff_inputs ; 
  delete [] J ; 
  
  delete [] n_per_pop ;
  delete [] cum_n_per_pop ; 
  
  file_volt.close() ; 
  file_spike_times.close() ; 
  
  file_mean_rates.close() ; 
  file_filter_rates.close() ; 
  file_inputs.close() ; 

  //STP
  delete [] u_stp ;
  delete [] x_stp ;

  file_u_stp.close(); 
  file_x_stp.close(); 

  // LOW RANK  
  delete [] xi ;
  delete [] xi_left ; 
  delete [] xi_right ;
  
  delete [] kappa ;
  
  file_overlap.close();
  
  cout << "done" << endl ; 
  
  ///////////////////////////////////////////////////////////////////
  
  cout << "Simulation Done !" << endl ; 
  
  clock_t t2=clock() ; 
  int HOURS=0,MIN=0,SEC=0;
  string str_TIME = path + "/CPU_TIME.txt" ; 
  ofstream TIME(str_TIME.c_str(), ios::out | ios::ate); 

  SEC = (t2-t1)/CLOCKS_PER_SEC ;
  HOURS = SEC/3600 ;
  MIN = SEC/60 ;
  SEC = SEC % 60 ;
  cout << "Elapsed Time = " << HOURS << "h " << MIN << "m " << SEC << "s" << endl;
  TIME << "Elapsed Time = " << HOURS << "h " << MIN << "m " << SEC << "s" << endl;
  TIME.close() ;
  return 0;

}
