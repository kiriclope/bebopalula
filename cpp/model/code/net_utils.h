#ifndef __NETUTILS__
#define __NETUTILS__

string str_mean_rates , str_filter_rates, str_inputs, str_overlaps ; 
ofstream file_mean_rates, file_filter_rates, file_inputs, file_overlaps ;

void print_sim_info() {
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << endl ; 

  if(IF_LIF)
    cout << "LIF model " ;
  if(IF_BIN) 
    cout << "Binary model " ; 
  if(IF_RATE)
    cout << "Rate model " ; 
      
  if(IF_RK2) 
    cout << "with RK2 scheme and interpolation " ; 
  else
    cout << "with EULER integration " ; 
      
  if(IF_STP) 
    cout << "and STP " ; 
  
  cout << endl ; 
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << endl ;

  cout << "path: " << path << endl;
  
  cout << "Main loop :" ; 
  cout << " duration " << DURATION ; 
  cout << " ms | DT " << DT ; 
  cout << " ms | TIME_STEADY " << TIME_STEADY ; 
  cout << " ms | TIME_WINDOW " << TIME_WINDOW ; 
  cout << " ms | TIME_REC " << TIME_REC ; 
  cout << " ms | TIME_REC_SPIKES " << TIME_REC_SPIKES ; 
  cout << endl ;
  
  cout << "mean field rates: " ; 
  if(IF_LIF) { 
    for(i=0; i<n_pop; i++) 
      cout << mf_rates[i] << " " ; 
    cout << endl ;
  }
  else {
    for(i=0; i<n_pop; i++) 
      cout << mf_rates[i] << " " ; 
    cout << endl ;
  } 
  
}

void open_files() {
  str_mean_rates = path + "/mean_rates.dat" ; 
  str_filter_rates = path + "/filter_rates.dat" ; 
  str_inputs = path + "/inputs.dat" ; 
  str_overlaps = path + "/overlaps.dat" ; 
  
  file_mean_rates.open(str_mean_rates.c_str(), ios::out | ios::ate) ; 
  file_filter_rates.open(str_filter_rates.c_str(), ios::out | ios::ate) ; 
  file_inputs.open(str_inputs.c_str(), ios::out | ios::ate) ;   
  file_overlaps.open(str_overlaps.c_str(), ios::out | ios::ate) ; 
}

void close_files() { 
  file_mean_rates.close() ; 
  file_filter_rates.close() ; 
  file_inputs.close() ; 
  file_overlaps.close() ; 
}


void init_globals() { 
  m0 = M0 ; 

  if(n_pop==1) {
    if(IF_SINGLE_NEURON)
      n_neurons = (unsigned long) n_neurons ; 
    else 
      n_neurons = (unsigned long) (n_neurons * 10000) ; 
  }
  else 
    n_neurons = (unsigned long) (n_neurons * 10000) ; 
    
  sqrt_K = sqrt( (double) K) ; 
  
  Ka = new double [n_pop]() ; 
  sqrt_Ka = new double [n_pop]() ;
  
  for(i=0;i<n_pop;i++) 
    if(IF_MATO_K) {
      Ka[i] = K * n_frac[i] ;
      sqrt_Ka[i] = sqrt( Ka[i] ) ;      
    }
    else {
      Ka[i] = K ; 
      sqrt_Ka[i] = sqrt_K ; 
    }
  
  n_per_pop = new unsigned long [n_pop]() ; 
  
  for(i=0;i<n_pop;i++) 
    n_per_pop[i] = (unsigned long) ( n_frac[i] * (double) n_neurons ) ; 
  
  cum_n_per_pop = new unsigned long [n_pop+1]() ; 
  
  cout <<"cum_n_per_pop=" << " " ; 
  for(i=0;i<n_pop+1;i++) { 
    for(j=0;j<i;j++) 
      cum_n_per_pop[i] += n_per_pop[j] ; 
    cout << cum_n_per_pop[i] << " " ;
  }
  cout << endl ; 
  
  which_pop = (int *) malloc( (unsigned long) n_neurons * sizeof(int) ) ; 
  
  for(i=0;i<n_pop;i++) 
    for(j=0; j<n_neurons; j++)
      if(j>=cum_n_per_pop[i] && j<cum_n_per_pop[i+1])
	which_pop[j] = i ; 
  
  mean_rates = new double [n_pop]() ; // population averaged rate also averaged over TIME_WINDOW  
  filter_rates = new double [n_neurons]() ; // temporal averaged over TIME_WINDOW

  ISI = new double [n_neurons]() ; // temporal averaged over TIME_WINDOW 
  
  // h^(ab)_i=h^b_i, inputs from presynaptic population b to postsynaptic neuron (i,a)
  // here i goes threw all the populations 
  inputs = new double *[n_pop]() ; 
  for(i=0;i<n_pop;i++) // presynaptic population b 
    inputs[i] = new double [n_neurons]() ; 

  if(IF_NMDA) {
    inputs_nmda = new double *[n_pop]() ; 
    for(i=0;i<n_pop;i++) // presynaptic population b 
      inputs_nmda[i] = new double [n_neurons]() ; 
  }
  
  filter_inputs = new double *[n_pop]() ; 
  for(i=0;i<n_pop;i++) // presynaptic population b 
    filter_inputs[i] = new double [n_neurons]() ; 
  
  // htot_i = h^E_i + h^I_i, net input into neuron i
  net_inputs = new double [n_neurons]() ;   

  ff_inputs = new double [n_neurons]() ; 
  
  net_inputs_RK2 = new double [n_neurons]() ;   
  
  mf_rates = new double [n_pop]() ; 
  mean_field_rates() ; 

  if(IF_LOW_RANK) {
    overlaps = new double [n_pop]() ; 
    ksi = new double [n_neurons]() ;
    ksi_scaled = new double [n_neurons * n_neurons]() ; 
    
    shared_ksi  = new double [n_neurons]() ; 
    for(j=0; j<n_neurons; j++)
      shared_ksi[j] = white_noise(covar_ksi_gen) ; 

    shared_sample  = new double [n_neurons]() ; 
    for(j=0; j<n_neurons; j++)
      shared_sample[j] = white_noise(ksi_gen) ; 

    shared_dist  = new double [n_neurons]() ; 
    for(j=0; j<n_neurons; j++)
      shared_dist[j] = white_noise(ksi_1_gen) ; 
    
    sample = new double [n_neurons]() ; 
    distractor = new double [n_neurons]() ; 
    
    ksi_init = new double [n_neurons]() ; 
    if(RANK==2) { 
      ksi_1 = new double [n_neurons]() ; 
      ksi_1_scaled = new double [n_neurons * n_neurons]() ; 
    }
    
  }
  
  if(IF_GAUSS)
    prefactor = new double [n_pop*n_neurons]() ; 
  
  if(IF_SPEC)
    if(RANK==2) {
      idx_perm_E = new unsigned long [n_per_pop[0]]() ; 
      /* idx_perm_I = new unsigned long [n_per_pop[1]]() ;  */
      idx_perm = new unsigned long [n_neurons]() ; 
      theta_1 = new double [n_neurons]() ; 
    }
  
  ext_inputs_scaled = new double [n_pop]() ; 
  J_scaled = new double [n_pop*n_pop]() ; // instantaneous individual membrane voltage 
  J_nmda = new double [n_pop*n_pop]() ; // instantaneous individual membrane voltage 
  
}

void delete_globals() { 
  delete [] mean_rates ; 
  delete [] filter_rates ; 
  delete [] ISI ;
  
  delete [] inputs ;
  if(IF_NMDA)
    delete [] inputs_nmda ;
  
  delete [] filter_inputs ; 
  delete [] net_inputs ; 
  delete [] net_inputs_RK2 ; 
  
  delete [] ext_inputs ; 
  delete [] ff_inputs ; 
  delete [] J ; 
  delete [] J_scaled ;
  delete [] J_nmda ;
  
  delete [] n_per_pop ;
  delete [] cum_n_per_pop ;
  free(which_pop) ;
  
  delete [] sqrt_Ka ;
  delete [] Ka ;
  
  if(IF_LOW_RANK) {
    delete [] ksi ; 
    delete [] ksi_init ;
    delete [] ksi_scaled ;

    delete [] shared_ksi ; 
    delete [] shared_sample ; 
    delete [] shared_dist ; 
    
    if(RANK==2) {
      delete [] ksi_1 ;
      delete [] ksi_1_scaled ; 
    }
    delete [] overlaps ;
  }

  if(IF_GAUSS) 
    delete [] prefactor ; 
  if(IF_SPEC) {
    delete [] theta ; 
    if(RANK==2) {
      delete [] idx_perm ;
      delete [] idx_perm_E ;
      delete [] idx_perm_I ;
      delete [] theta_1 ;       
    }
  }
}

void make_dir(string path) {
  string mkdirp = "mkdir -p " ; 
  mkdirp += path ; 
  
  const char * cmd = mkdirp.c_str() ; 
  const int dir_err = system(cmd) ; 
  
  cout << path << endl ; 
}

#define VarToString(name) var_to_string(#name, (name))

template <class T>
const char* var_to_string(const char *name, const T& val) {
  // cout << name << " = " << val << endl;
  return name ;
}

void get_args(int argc , char** argv) { 
  if(argv[1] != NULL) { 
    n_pop = (int) atoi(argv[1]) ; 
    n_neurons = (unsigned long) atoi(argv[2]) ; 
    K = (double) atof(argv[3]) ;
    if(n_pop!=1)
      dir = argv[4] ;
  }
  else {
    cout << "n_pop ? " ;
    cin >> n_pop ;
    cout << "n_neurons ? " ;
    cin >> n_neurons ;
    cout << "K ? " ;
    cin >> K ;
    cout << "Directory ? " ;
    cin >> dir ;
  } 
}

///////////////////////////////////////////////////////////////////////

void get_param() {

  cout << "reading parameters from : " ;
  string file_name = "../parameters/" + to_string(n_pop) + "pop/" + dir +".txt" ; 
  cout << file_name << endl; 
  
  ext_inputs = new double [n_pop]() ;
  J = new double [n_pop*n_pop]() ;
  
  double *Tsyn ;
  Tsyn = new double [n_pop*n_pop]() ; 
  
  string token ;
  string::size_type sz;
  ifstream file(file_name.c_str()) ; 
  int i,j;
  
  i=0 ; 
  while(getline(file, token)) {
    j=-1 ; 
    istringstream line(token);
    while(line >> token) {
      if(i==0)
	if(j!=-1) ext_inputs[j] = stod(token, &sz) ; 
      
      if(i==1)
	if(j!=-1) J[j] = stod(token, &sz) ; 

      if(i==2)
	if(j!=-1) Tsyn[j] = stod(token, &sz) ; 
      
      j++ ;
    }
    
    if(file.unget().get() == '\n') {
      i++ ;
    }
  }

  if(n_pop==1) { 
    J[0] = J0 ; 
    ext_inputs[0] = I0 ; 
  }
  
  cout << "ext_inputs" << endl ;
  for(int i=0;i<n_pop;i++)
    cout << ext_inputs[i] << " " ;
  cout << endl ;
  
  cout << "J" << endl ;
  for(int i=0;i<n_pop;i++) {
    for(int j=0;j<n_pop;j++)
      cout << J[j+i*n_pop] << " " ;
    cout << endl ;
  }
  
  cout << "Tsyn" << endl ;
  for(int i=0;i<n_pop;i++) {
    for(int j=0;j<n_pop;j++)
      cout << Tsyn[j+i*n_pop] << " " ;
    cout << endl ;
  }
    
}

double Phi(double x) { // Gaussian CDF
  return 0.5 * ( 1.0 + erf(x/sqrt(2.0) ) ) ; 
}

double threshold_linear(double x) {
  if(x>0.0) 
    return x ; 
  else 
    return 0. ; 
}

double cut_LR(double x) {
  if(eps[2-n_pop]*x>0.0 && x<=1.0) 
    return x ; 
  else 
    return 0. ; 
}

void create_dir() { 
  
  path += "simulations/" ; 
  
  if(IF_LIF) 
    path += "lif/" ; 
  if(IF_BIN) 
    path += "binary/" ; 
    
  ostringstream str_I0, str_J0 ;
  str_I0 << fixed << setprecision(2) << I0 ;
  str_J0 << fixed << setprecision(2) << abs(J0) ;
  
  if(n_pop==1)
    dir = "I0_" + str_I0.str()  + "_J0_" + str_J0.str() ;
  
  path += to_string(n_pop)+"pop/" + dir  ;
  /* path += "/N" + to_string(n_neurons) ; */
  
  if(n_pop==1)
    path += "/N" + to_string(n_per_pop[0]/1000) ; 
  else 
    path += "/NE_" + to_string(n_per_pop[0]/1000) +  "_NI_" + to_string(n_per_pop[1]/1000) ; 
  
  path += "/K" + to_string((int)K) ; 
  
  ostringstream str_tau_fac, str_tau_rec, str_use ;
  str_tau_fac << fixed << setprecision(0) << TAU_FAC ; 
  str_tau_rec << fixed << setprecision(0) << TAU_REC ; 
  str_use << fixed << setprecision(2) << USE ; 

  if(IF_STP) 
    path += "/STP/Tf_" + str_tau_fac.str() + "_Tr_" +  str_tau_rec.str() + "_U_" +  str_use.str() ;  
  
  ostringstream str_kappa, str_kappa_var, str_kappa_1 ;
  str_kappa << fixed << setprecision(2) << KAPPA ; 
  str_kappa_var << fixed << setprecision(2) << KAPPA_VAR ; 
  str_kappa_1 << fixed << setprecision(2) << KAPPA_1 ; 
  
  ostringstream str_ksi, str_ksi_var, str_ksi_1 , str_ksi_var_1 ; 
  str_ksi << fixed << setprecision(2) << MEAN_KSI ; 
  str_ksi_var << fixed << setprecision(2) << VAR_KSI ; 
  str_ksi_1 << fixed << setprecision(2) << MEAN_KSI ; 
  str_ksi_var_1 << fixed << setprecision(2) << VAR_KSI ; 
 
  ostringstream str_map_seed ; 
  str_map_seed << fixed << setprecision(0) << MAP_SEED ; 

  ostringstream str_EE, str_EI, str_IE, str_II ;  
  str_EE << fixed << setprecision(0) << SIGMA[0] ; 
  str_EI << fixed << setprecision(0) << SIGMA[1] ; 
  str_IE << fixed << setprecision(0) << SIGMA[2] ; 
  str_II << fixed << setprecision(0) << SIGMA[3] ; 
    
  if(IF_STRUCTURE) { 
    if(IF_SPEC) {
      if(RANK==1)
	path += "/spec/kappa_" + str_kappa.str() ;
      if(RANK==2) {
	path += "/spec/kappa_" + str_kappa.str() + "_kappa_1_" + str_kappa_1.str() ;
	if(FIX_MAP_SEED) 
	  path += "/seed_" + str_map_seed.str() ; 
      }
    }
    
    if(IF_LOW_RANK) {
      if(RANK==1){
	path += "/low_rank/kappa_" + str_kappa.str() ; 
	/* path += "/mean_" + str_ksi.str() + "_var_" + str_ksi_var.str() ; */ 
      }
      if(RANK==2) { 
	path += "/low_rank/kappa_" + str_kappa.str() + "_kappa_1_" + str_kappa_1.str() ; 
	/* path += "/mean_" + str_ksi.str() + "_var_" + str_ksi_var.str() ;  */
	/* path += "/mean_" + str_ksi_1.str() + "_var_" + str_ksi_var_1.str() ;  */
      }
    }
  
    if(IF_RING)
      path += "/ring/kappa_" + str_kappa.str() ;
        
    if(IF_GAUSS)
      path += "/gauss/EE_" + str_EE.str() + "_EI_" + str_EI.str() +"_IE_" + str_IE.str() + "_II_" + str_II.str() ;     
  } 
  
  ostringstream str_kappa_ext, str_phi_ext ; 
  str_kappa_ext << fixed << setprecision(2) << KAPPA_EXT ; 
  str_phi_ext << fixed << setprecision(2) << PHI_EXT ; 
  
  if(IF_DPA) 
    path += "/DPA/kappa_" + str_kappa_ext.str() + "_phi_" + str_phi_ext.str() ; 
  if(IF_DUAL) 
    path += "/dual_task/kappa_" + str_kappa_ext.str() + "_phi_" + str_phi_ext.str() ; 
  if(IF_DRT) 
    path += "/DRT/kappa_" + str_kappa_ext.str() + "_phi_" + str_phi_ext.str() ; 
  
  if(IF_HYSTERESIS) {
    if(HYST_J_EE==1) 
      path += "/hysteresis/Jee_up"; 
    if(HYST_J_EE==-1) 
      path += "/hysteresis/Jee_down"; 

    if(HYST_M0==1) 
      path += "/hysteresis/m0_up"; 
    if(HYST_M0==-1) 
      path += "/hysteresis/m0_down"; 
  }

  ostringstream str_m0 ; 
  str_m0 << fixed << setprecision(4) << M0 ; 
  if(IF_LOOP_M0) 
    path += "/m0_" + str_m0.str() ; 

  ostringstream str_gain ; 
  str_gain << fixed << setprecision(1) << GAIN ; 
  if(IF_LOOP_GAIN) 
    path += "/gain_" + str_gain.str() ; 
  
  if(IF_TRIALS) 
    path += "/trial_" + to_string( (int) TRIAL_ID ) ;
  
  if(IF_INI_COND) 
    path += "/ini_cond_" + to_string( (int) INI_COND_ID ) ; 
    
  make_dir(path) ; 
  
  cout << "Created directory : " ; 
  cout << path << endl ; 

}

template <class T> 
void read_from_file(string path, string file_name, T * &array, size_t array_size) {

  int dum ; 
  
  string file_path ; 
  file_path = path + "/" + file_name + ".dat" ; 
  cout << "reading from: " << file_path << endl ; 
  
  struct stat buffer ; 
  FILE *file ; 
  
  if (stat (file_path.c_str(), &buffer) == 0) { 
    file = fopen(file_path.c_str(), "rb") ; 
    dum = fread(&array[0], sizeof array[0], array_size, file) ; 
    fclose(file) ; 
  } 
  else { 
    cout << "ERROR: " << file_name << ".dat not found" << endl ; 
    exit(-1) ; 
  }
}

template <class T> 
void write_to_file(string path, string file_name, T * &array, size_t array_size) {

  int dum ; 
  
  string file_path ; 
  file_path = path + "/" + file_name + ".dat" ; 
  cout << "writing to: " << file_path << endl ; 
  
  struct stat buffer ; 
  FILE *file ; 
  
  file = fopen(file_path.c_str(), "wb") ; 
  dum = fwrite(&array[0], sizeof array[0], array_size, file) ; 
  fclose(file) ;
  
}

void save_to_file() { 

  if(!IF_HYSTERESIS) 
    file_mean_rates << t_time - TIME_STEADY ; 
  else {
    if(HYST_J_EE!=0)
      file_mean_rates << J[0] ;
    if(HYST_M0!=0)
      file_mean_rates << m0 ;
  }
  
  for(i=0;i<n_pop;i++) { 
    file_mean_rates << " " << mean_rates[i]*1000./TIME_WINDOW/(double)n_per_pop[i] ; 
    mean_rates[i] = 0 ; 
  }
  
  file_mean_rates << endl ; 
  
  if(IF_LOW_RANK) {
    file_overlaps << t_time - TIME_STEADY ; 
    for(i=0;i<n_pop;i++) { 
      file_overlaps << " " << overlaps[i] * 1000. / TIME_WINDOW * IS_STRUCT_SYN[i] ; 
      /* file_overlaps << " " << overlaps[i]*1000./TIME_WINDOW/(double)n_per_pop[i] ;  */
      overlaps[i] = 0 ; 
    } 
    file_overlaps << endl ; 
  }
  
  // filtered rates over tw

  if(!IF_HYSTERESIS)
    file_filter_rates << t_time - TIME_STEADY ; 
  else {
    if(HYST_J_EE!=0)
      file_filter_rates << J[0] ;
    if(HYST_M0!=0)
      file_filter_rates << m0 ;     
  }
  
  for(i=0;i<n_neurons;i++) { 
    file_filter_rates << " " << filter_rates[i]*1000./TIME_WINDOW ; 
    filter_rates[i] = 0 ; 
  } 
  file_filter_rates << endl ; 
  
  // filtered inputs over tw 
  if(!IF_HYSTERESIS)
    file_inputs << t_time - TIME_STEADY ;
  else {
    if(HYST_J_EE!=0)
      file_inputs << J[0] ;
    if(HYST_M0!=0)
      file_inputs << m0 ; 
  }
  for(i=0;i<n_pop;i++) 
    for(j=0;j<n_neurons;j++) { 
      file_inputs << " " << filter_inputs[i][j]*DT/TIME_WINDOW ; 
      filter_inputs[i][j] = 0 ; 
    } 
  file_inputs << endl ;   
}

void get_m1_phase() { 
  
  double dPhi = 0 ; 
  double xCord = 0, yCord = 0 ; 
  
  for(int i_pop=0; i_pop<n_pop; i_pop++) {
    
    dPhi = M_PI / (double) n_per_pop[i_pop] ; 
    xCord = 0;
    yCord = 0 ; 
    
    for(unsigned long j=cum_n_per_pop[i_pop]; j < cum_n_per_pop[i_pop+1]; j++) {
      xCord += filter_rates[j] * cos(2.0 * j * dPhi) / TIME_WINDOW ; 
      yCord += filter_rates[j] * sin(2.0 * j * dPhi) / TIME_WINDOW ; 
    }
    
    m1[i_pop] = ( 2.0 / (double) n_per_pop[i_pop]) * sqrt(xCord * xCord + yCord * yCord) ; 
    phase[i_pop] = 0.5 * atan2(yCord, xCord) ; 
    
    if(phase[i] < 0)
      phase[i] = phase[i] + M_PI ;
    
    phase[i] *= 180.0/M_PI ; 
  }
}

#endif
