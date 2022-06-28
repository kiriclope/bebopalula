#ifndef _UTILS_ 
#define _UTILS_

double norm_array(double *array, size_t n) { 
  
  double norm = 0 ; 
  for(i=0; i<n; i++) 
    norm += array[i]*array[i] ;
  
  return sqrt(norm) ;
}

double cos_array(double *a, double *b, size_t n) {
  double cosine_ab = 0, norm_a=0, norm_b=0 ; 
  for(i=0; i<n; i++) 
    cosine_ab += a[i]*b[i] ; 

  norm_a = norm_array(a, n) ;
  norm_b = norm_array(b, n) ;

  return cosine_ab / norm_a / norm_b ; 
}

void outer_product(double *a, double *b, size_t n, double *&outer) { 
  
  outer = new double [n*n]() ; 
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      outer[i+j*n] = a[i]*b[j] ; 
  
}

void normalize_array(double *&a, size_t n) {
  double norm = 0 ;
  norm = norm_array(a,n) ;
  
  for(i=0; i<n; i++)
    a[i] /= norm ; 
}
  
void rotate_ab(double* a, double* b, size_t n, double angle, double *&result) {

  normalize_array(a, n) ; 
  normalize_array(b, n) ; 
  
  double *Id;
  Id = new double [n*n]() ; 
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      Id[i+j*n] = 1.0 ; // Identity Matrix 

  double *outer_ab, *outer_ba, *outer_aa, *outer_bb; 
  outer_product(a, a, n, outer_aa) ;
  outer_product(b, b, n, outer_bb) ;  
  outer_product(a, b, n, outer_ab) ;
  outer_product(b, a, n, outer_ba) ;
  
  /* R = I + ( np.outer(n2,n1) - np.outer(n1,n2) ) * np.sin(a) + ( np.outer(n1,n1) + np.outer(n2,n2) ) * (np.cos(a)-1) ; */
  
  double *R;
  R = new double [n*n]() ;
  
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      R[i+j*n] = Id[i+j*n] + ( outer_ba[i+j*n] - outer_ab[i+j*n] ) * sin(angle) + ( outer_aa[i+j*n] + outer_bb[i+j*n] ) * (cos(angle) - 1.0) ;

  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
      result[i] += R[i+j*n] * b[j] ; 
  
  delete [] Id ;
  delete [] R ; 
}

void angle_ksi() { 
   
  double cos_ksi_0 = 0.0 ;
  cos_ksi_0 = cos_array(ksi, ksi_1, n_per_pop[0]) ; 

  cout << "cos_ksi_0: " << cos_ksi_0 ; 
  cout << " angle_ksi_0: " << acos(cos_ksi_0) * 180.0 / M_PI << "°" << endl ; 
  
  /* double norm_ksi=0, norm_ksi_1=0, dot=0 ;  */
  /* norm_ksi = norm_array(ksi, n_per_pop[0]) ;  */
  /* norm_ksi_1 = norm_array(ksi, n_per_pop[0]) ;  */
  
  /* for(i=0; i<n_per_pop[0]; i++)  */
  /*   dot += ksi[i]*ksi_1[i] ;  */
  
  /* for(i=0; i<n_per_pop[0]; i++)  */
  /*   ksi_1[i] = ksi_1[i] - dot/norm_ksi/norm_ksi * ksi[i] ; // Gram Schmidt ortho  */
  /*   /\* ksi_1[i] = ksi_1[i] - dot/norm_ksi_1/norm_ksi_1 * ksi_1[i] ; // Gram Schmidt ortho  *\/  */
  
  /* double cos_ksi = 0.0 ;  */
  /* cos_ksi = cos_array(ksi, ksi_1, n_per_pop[0] ) ;  */
  
  /* cout << "cos_ksi: " << cos_ksi ;  */
  /* cout << " angle_ksi: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ; */
  
  /* double *ksi_rotate; */
  /* ksi_rotate = new double [n_per_pop[0]]() ;  */
  
  /* rotate_ab( ksi, ksi_1, n_per_pop[0], ANGLE_KSI, ksi_rotate) ;  */
  
  /* cos_ksi = cos_array(ksi, ksi_rotate, n_per_pop[0] ) ;  */
  
  /* cout << "cos_ksi: " << cos_ksi ;  */
  /* cout << " angle_ksi: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ; */
  
}

void angle_maps() { 

  cout << "theta: " ; 
  for(i=0; i<10; i++) 
    cout << theta[i] << " " ; 
  cout << endl ;
  
  cout << "theta_1: " ; 
  for(i=0; i<10; i++) 
    cout << theta_1[i] << " " ; 
  cout << endl ; 
  
  double sum_theta=0.0, norm_theta ; 
  
  for(i=0; i<n_per_pop[0]; i++) 
    sum_theta += theta[i] ; 
  
  cout << "sum theta: " << sum_theta << endl ;
  double sum_theta_1=0.0, norm_theta_1 ;
  
  for(i=0; i<n_per_pop[0]; i++) { 
    theta_1[i] = theta_1[i] - theta[i] * theta_1[i] / sum_theta ; 
    theta_1[i] = theta_1[i] + cos(MAP_ANGLE) / sum_theta ; 
    /* sum_theta_1 += theta_1[i] ;  */ 
  }

  cout << "theta_1 ortho: " ; 
  for(i=0; i<10; i++) 
    cout << theta_1[i] << " " ; 
  cout << endl ; 
  
  double cos_maps = 0.0 ; 
  for(i=0; i<n_per_pop[0]; i++) {
    /* theta_1[i] += cos(MAP_ANGLE) / sum_theta_1 ;  */ 
    cos_maps += theta[i]*theta_1[i] ; 
  }
  
  norm_theta = norm_array(theta, n_per_pop[0]) ;
  norm_theta_1 = norm_array(theta_1, n_per_pop[0]) ; 
  
  cos_maps = cos_maps / norm_theta / norm_theta_1 ; 
  
  cout<< "norm_theta " << norm_theta << " norm_theta_1 " << norm_theta_1 << endl ; 
  cout << "cos_maps: " << cos_maps << ", angle_maps: " << acos(cos_maps) * 180.0 / M_PI << "°" << endl ;
}

double Gaussian1D(double mu, double sigma) {

  mu = min(abs(mu), 2.0*M_PI-abs(mu)) ;  // mato et al.
  sigma *= DEG_TO_RAD ; 
  
  if(sigma!=0.) 
    return exp(-mu*mu/2./sigma/sigma)/sqrt(2.*M_PI)/sigma ; 
  else 
    return 1. ; 
}

void my_shuffle(void *base, size_t n, size_t size) { 
  
  const gsl_rng_type * T ; 
  gsl_rng * r ; 
  T = gsl_rng_default ; 
  r = gsl_rng_alloc (T) ;
  
  if(!FIX_MAP_SEED)
    gsl_rng_set(r, clock());
  else {
    cout << "MAP SEED " << exp( (double) MAP_SEED) << endl ; 
    gsl_rng_set(r, exp( (double) MAP_SEED) ); 
  }
  gsl_ran_shuffle (r, base, (size_t) n, size) ; 
  
}

size_t* permutation(size_t N) {

  FILE *pFile ; 
  string file_path ; 
  
  /* file_path = con_path + "idx_perm.dat" ;  */
  /* cout << file_path << endl ;  */
  
  /* pFile = fopen(file_path.c_str(), "wb") ;  */
  
  const gsl_rng_type * T ; 
  gsl_rng * r ; 
  gsl_permutation * p = gsl_permutation_alloc (N) ; 
  
  gsl_rng_env_setup(); 
  T = gsl_rng_default; 
  r = gsl_rng_alloc (T); 
  gsl_rng_set(r, clock()); 
  
  gsl_permutation_init (p); 
  gsl_ran_shuffle (r, p->data, N, sizeof(size_t)) ; 
  /* gsl_permutation_fprintf (pFile, p, " %u");  */
  
  /* fclose(pFile) ;  */
  return p->data ; 
}

void init_con_globals() { 
  cout << "initialize globals" << endl ; 
  
  IF_STRUCTURE = IF_RING || IF_SPEC || IF_LOW_RANK || IF_GAUSS ; 
  
  n_per_pop = new unsigned long [n_pop]() ; 
  
  for(i=0;i<n_pop;i++) 
    n_per_pop[i] = (unsigned long) ( n_frac[i] * (double) n_neurons ) ; 
  
  cum_n_per_pop = new unsigned long [n_pop+1]() ; 
  
  for(i=0;i<n_pop+1;i++) 
    for(j=0;j<i;j++) 
      cum_n_per_pop[i] += n_per_pop[j] ; 
  
  which_pop = new int [n_neurons] ; 
  
  for(i=0;i<n_pop;i++) 
    for(j=0; j<n_neurons; j++) 
      if(j>=cum_n_per_pop[i] && j<cum_n_per_pop[i+1]) 
	which_pop[j] = i ; 
  
  prefactor = new double [n_pop*n_neurons]() ;
  
}

void delete_con_globals() { 
  cout << "delete globals" << endl ; 
  
  delete [] n_per_pop ; 
  delete [] cum_n_per_pop ; 
  delete [] which_pop ; 
  
  /* delete [] con_prob ;  */
  /* delete [] con_vec ;  */
  
  delete [] n_post ; 
  delete [] idx_post ;
  free(id_post) ; 
  
  delete [] prefactor ;
  
  delete [] X ; 
    
}

void init_X() { 
  
  X = new double [n_neurons]() ; 
  
  for(i=0; i<n_pop; i++) 
    for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++)      
      X[j] = fmod( (double) (j-cum_n_per_pop[i]) , (double) n_per_pop[i] ) * 2.0 * M_PI / (double) n_per_pop[i] ;     
  
  cout << "X: " ; 
  for(i=0; i<10; i++) 
    cout << X[i] << " " ; 
  cout << endl ; 
  
}

void init_theta() {
  
  theta = new double [n_neurons]() ; 
  
  for(i=0; i<n_pop; i++) { 
    for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) 
      theta[j] = ( 2.0 * M_PI * (double) (j-cum_n_per_pop[i]) ) / (double) n_per_pop[i] ; 
  } 
  
  cout << "theta: " ; 
  for(i=0; i<10; i++) 
    cout << theta[i] << " " ; 
  cout << endl ; 
  
}

double distance_between_maps() {
  double sum = 0 ;
  for(i=0; i<n_per_pop[0]; i++) 
    if( fmod(2.0 * M_PI * (double) (idx_perm_E[i]-i) / (double) n_per_pop[0], 2.0*M_PI) >= M_PI/4.0 ) 
      sum += 1 ; 
  
  return sum / (double) n_per_pop[0] ; 
}

void init_theta_1() {
  
  /* idx_perm = permutation( (size_t) n_per_pop[0] ) ; */
  
  for(i=0; i<n_per_pop[0]; i++) 
    idx_perm_E[i] = i ; 
  
  my_shuffle(idx_perm_E, n_per_pop[0], sizeof(unsigned long) ) ; 
  
  /* double distance = 0 ; */
  /* distance = distance_between_maps() ; */
  /* cout << "distance between maps " << distance << endl ; */
  
  cout << "idx_perm_E: " ; 
  for(i=0;i<10;i++) 
    cout << idx_perm_E[i] << " " ; 
  cout << endl ; 
  
  /* for(i=0; i<n_per_pop[1]; i++)  */
  /*   idx_perm_I[i] = i + n_per_pop[0] ;  */
  
  /* my_shuffle(idx_perm_I, n_per_pop[1], sizeof(unsigned long) ) ;  */
  
  /* cout << "idx_perm_I: " ;  */
  /* for(i=0;i<10;i++)  */
  /*   cout << idx_perm_I[i] << " " ;  */
  /* cout << endl ;  */
  
  for(i=0; i<n_pop; i++)
    if(i==0)
      for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) 
	idx_perm[j] = idx_perm_E[j] ; 
    /* else  */
    /*   for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++)  */
    /* 	idx_perm[j] = idx_perm_I[j-cum_n_per_pop[i]] ;  */
  
  cout << "idx_perm: " ; 
  for(i=0;i<10;i++) 
    cout << idx_perm[i] << " " ; 
  cout << endl ; 
  
  for(i=0; i<n_pop; i++) 
    for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++)
      theta_1[j] = theta[idx_perm[j]] ; 
  
  /* if(i==0)  */
  /* 	theta_1[j] = 2.0 * M_PI * (double) (idx_perm[j]-cum_n_per_pop[i]) / (double) n_per_pop[i] ;  */
  /* else  */
  /* 	theta_1[j] = 2.0 * M_PI * (double) (j-cum_n_per_pop[i]) / (double) n_per_pop[i] ;  */
  
  /* angle_maps() ;  */


  cout << "theta_1: " ; 
  for(i=0; i<10; i++) 
    cout << theta_1[i] << " " ; 
  cout << endl ; 
  
}

void init_ksi() { 
  
  double norm ;
  
  if(FIX_KSI_SEED) 
    for(i=0; i<n_pop; i++) { 
      norm = 0 ; 
      for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) { 
	ksi[j] = MEAN_KSI ;
	  /* + sqrt(abs(VAR_KSI-COVAR_KSI-COVAR_KSI_SAMPLE-COVAR_KSI_DIST)) * white_noise(ksi_gen)  */
	  /* + sqrt(COVAR_KSI) * shared_ksi[j]  */
	  /* + sqrt(abs(VAR_KSI-COVAR_KSI_SAMPLE-COVAR_KSI_DIST)) * white_noise(ksi_gen)  */
	  /* + sqrt(COVAR_KSI_SAMPLE) * shared_sample[j]  */
	  /* - sqrt(COVAR_KSI_DIST) * shared_dist[j] ;  */
	norm += ksi[j]*ksi[j] ; 
      } 
      
      if(IF_NORM_KSI) 
	for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) 
	  ksi[j] /= sqrt(norm) ; 
    } 
  else
    for(i=0; i<n_pop; i++) {
      norm = 0 ;
      for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) { 
	ksi[j] = MEAN_KSI + sqrt(VAR_KSI) * white_noise(rand_gen) ; 
	norm += ksi[j]*ksi[j] ; 
      } 
      
      if(IF_NORM_KSI) 
	for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) 
	  ksi[j] /= sqrt(norm) ; 
    }
  
  cout << "ksi: " ;
  for(i=0; i<10; i++)
    cout << ksi[i] << " " ;
  cout << endl ;
  
} 

void init_ksi_1() { 

  double norm ;
  
  if(FIX_KSI_SEED)
    for(i=0; i<n_pop; i++) {
      norm = 0 ;
      for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) { 
	ksi_1[j] = MEAN_KSI_1 ;
	  /* + sqrt(abs(VAR_KSI_1-COVAR_KSI-COVAR_KSI1_SAMPLE-COVAR_KSI1_DIST)) * white_noise(ksi_1_gen)  */
	  /* + sqrt(COVAR_KSI) * shared_ksi[j]  */
	  /* + sqrt(abs(COVAR_KSI1_SAMPLE-COVAR_KSI1_DIST)) * shared_sample[j]  */
	  /* + sqrt(COVAR_KSI1_DIST) * shared_dist[j] ;  */
	norm += ksi_1[j]*ksi_1[j] ; 
      } 
      
      if(IF_NORM_KSI_1)
	for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) 
	  ksi_1[j] /= sqrt(norm) ; 
    } 
  else
    for(i=0; i<n_pop; i++) {
      norm = 0 ;
      for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) { 
	ksi_1[j] = MEAN_KSI_1 + sqrt(VAR_KSI_1) * white_noise(rand_gen) ; 
	norm += ksi_1[j]*ksi_1[j] ; 
      } 
      
      if(IF_NORM_KSI_1)
	for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i+1]; j++) 
	  ksi_1[j] /= sqrt(norm) ; 
    } 

  angle_ksi() ;
  
  cout << "ksi_1: " ;
  for(i=0; i<10; i++)
    cout << ksi_1[i] << " " ;
  cout << endl ;
  
} 

void func_con_prob() { 

  double kappa, kappa_K_N ; 
  double kappa_1, kappa_1_K_N ; 
  double cos_D_theta=0.0 ; 
    
  cout << "generate probability of connection: " ; 
  
  con_prob = new double [n_neurons*n_neurons]() ; 
  K_over_Na = new double [n_pop]() ; 
  
  for(i=0; i<n_pop; i++) { 
    /* K_over_Na[i] = ( K / (double) n_per_pop[i] ) ;  */ 
    if(IF_MATO_K) 
      K_over_Na[i] = (double) ( K * n_frac[i] / (double) n_per_pop[i] ) ;
    else
      K_over_Na[i] = (double) ( K / (double) n_per_pop[i] ) ; 
    cout << K_over_Na[i] << " " ; 
  } 
  cout << endl ; 
  
  if(IF_STRUCTURE) { 
    if(IF_MATO_K)
      kappa = KAPPA / sqrt(K * n_frac[0]) ; 
    else 
      kappa = KAPPA / sqrt_K ;
    
    if(RANK==2) 
      if(IF_MATO_K)
	kappa_1 = KAPPA / sqrt(K * n_frac[0]) ; 
    else
      kappa_1 = KAPPA / sqrt_K ; 
    
    if(IF_SPEC) { 
      cout << "random network with specific connections: rank, " << RANK << ", " ; 
      
      if(RANK==1) 
	cout << "kappa, " << KAPPA << endl ; 
      
      init_theta() ; 
      
      if(RANK==2) {
	cout << "kappa " << KAPPA << " kappa_1 " << KAPPA_1 << endl ;	
	init_theta_1() ; 
      }
      
      /* double phases[6] = {0, M_PI/4.0, M_PI/2.0, 3.0*M_PI/4.0, M_PI, 3.0*M_PI/2.0} ;  */
      
      for(j=0; j<n_neurons; j++) { 
	pre_pop = which_pop[j] ; 
	
	if(IF_KAPPA_DIST) { 
	  /* kappa = ( KAPPA * unif(rand_gen) ) / sqrt_K ;  */ 
	  kappa = ( KAPPA + sqrt(KAPPA_VAR) * white_noise(rand_gen) ) / sqrt_K ; 
	  kappa = log_normal(rand_gen) / sqrt_K ; 
	} 
	
	kappa_K_N = K_over_Na[pre_pop] * kappa ; 
	
	if(RANK==2) 
	  kappa_1_K_N = K_over_Na[pre_pop] * kappa_1 ; 
	
	for(i=0; i<n_neurons; i++) { 
	  post_pop = which_pop[i] ; 
	  
	  con_prob[j + i * n_neurons] += K_over_Na[pre_pop] ; 
	  
	  /* for(int i_phase=0; i_phase<6; i_phase++)  */
	  if(IS_STRUCT_SYN[pre_pop + post_pop * n_pop]) { 
	    cos_D_theta = cos(theta[i]-theta[j]) ; 
	    /* cos_D_theta = cos(theta[i]-theta[j] + phases[i_phase] ) ;  */ 
	    con_prob[j + i * n_neurons] += kappa_K_N * cos_D_theta ; 
	    if(RANK==2) 
	      con_prob[idx_perm[j] + idx_perm[i] * n_neurons] += kappa_1_K_N * cos_D_theta ; 
	    
	    if(con_prob[j + i * n_neurons]<=0 || con_prob[j + i * n_neurons]>=1) 
	      cout << "error con_prob>1 or <0"  << endl ; 
	  } 
	} 
      } 
      
    } //endif spec 
    
    if(IF_LOW_RANK) { 
      cout << "random network with low-rank specific connections: rank, " << RANK << ", " ;
      
      if(RANK==1)
	cout << "kappa, " << KAPPA << endl ;
      
      /* init_ksi() ; */  
      /* /\* kappa = KAPPA / sqrt_K  ; *\/ */
      
      if(RANK==2) {
	cout << "kappa, " << KAPPA << " kappa_1, " << KAPPA_1 << endl ;
    	/* init_ksi_1() ; */
    	/* /\* kappa_1 = KAPPA_1 / sqrt_K  ;  *\/ */ 	
      }
      
      for(j=0; j<n_neurons; j++) {
	pre_pop = which_pop[j] ; 
	
	kappa_K_N = K_over_Na[pre_pop] * kappa ; 
	if(RANK==2) 
	  kappa_1_K_N = K_over_Na[pre_pop] * kappa_1 ; 
	
	for(i=0; i<n_neurons; i++) { 
	  post_pop = which_pop[i] ; 
	  con_prob[j + i * n_neurons] += K_over_Na[pre_pop] ; 
	  if(IS_STRUCT_SYN[pre_pop + post_pop * n_pop]) { 
	    con_prob[j + i * n_neurons] += kappa_K_N * ksi[i] * ksi[j] ; 
	    if(RANK==2) 
	      con_prob[j + i * n_neurons] += kappa_1_K_N * ksi_1[i] * ksi_1[j] ; 
	    con_prob[j + i * n_neurons]	= cut_LR(con_prob[j + i * n_neurons]) ;	    
	  }
	  
	  /* if(con_prob[j + i * n_neurons]<=0 || con_prob[j + i * n_neurons]>=1)  */
	  /*   cout << "error con_prob>1 or <0"  << endl ;  */
	} 
      } 
    } // end low rank 
    
    if(IF_RING) { 
      cout << "ring: " ; 
      cout << "kappa, " << KAPPA << endl ;
      
      if(KAPPA>0.5) 
	cout << "ERROR: KAPPA TOO LARGE" << endl ; 
      
      init_theta() ;
      if(RANK==2) {
	cout << "kappa, " << KAPPA_1 << endl ; 
	init_theta_1() ; 
      }
      
      for(j=0; j<n_neurons; j++) {  
	pre_pop = which_pop[j] ;
	
	for(i=0; i<n_neurons; i++) { 
	  post_pop = which_pop[i] ;
	  
	  con_prob[j + i * n_neurons] += K_over_Na[pre_pop] ; 
	  /* if(IS_STRUCT_SYN[pre_pop + post_pop * n_pop])  */
	  con_prob[j + i * n_neurons] += 2.0 * K_over_Na[pre_pop] * kappas[pre_pop] * cos( theta[i] - theta[j] ) ; 
	  
	}
      }
    } // end ring

    if(IF_GAUSS) { 
      cout << "gaussian structure: " ; 
      /* cout << " L, " << L <<", sigma, " ; */
      cout << SIGMA[0] <<" "<< SIGMA[1] ; 
      cout <<" "<< SIGMA[2] <<" "<< SIGMA[3] <<  endl ; 
      
      /* init_X() ;  */
      init_theta() ;
      
      for(i=0; i<n_neurons; i++) { // need loop over post before pre to compute prefactor
	post_pop = which_pop[i] ; 
	
	for(j=0; j<n_neurons; j++) { // Jij: j (pre) to i (post)
	  pre_pop = which_pop[j] ;
	  
	  con_prob[j + i * n_neurons] = Gaussian1D(theta[i]-theta[j], SIGMA[pre_pop+post_pop*n_pop]) ; 
	  // Pij = Zb Cij with Zb = K / sum_j Cij then sum_j Pij = K 
	  prefactor[i + pre_pop*n_neurons] += con_prob[j + i * n_neurons] ; // sum over presynaptic neurons j 
	}
	
	for(j=0; j<n_neurons; j++) {
	  pre_pop = which_pop[j] ; 
	  con_prob[j + i * n_neurons] *= Ka[pre_pop] / prefactor[i+pre_pop*n_neurons] ; 
	}
	/* prefactor[pre_pop+post_pop*n_pop] = 0.0 ; // reinit prefactor for each post  */
	
      } // end loop post
      
    }// endif Gauss
    
  }// endif structure 
  
  else {
    cout << "random network" << endl ; 
    
    for(j=0; j<n_neurons; j++) { // Pij -> j (cols, pre) to i (rows, post) 
      pre_pop = which_pop[j] ; 
      for(i=0; i<n_neurons; i++) 
	con_prob[j + i * n_neurons] = K_over_Na[pre_pop] ; // Pij -> P[j + i * n_cols] (index = indexX * arrayWidth + indexY) ; 
    } 
  } 
  
  delete [] K_over_Na ; 
  
}

void func_con_vec() {
  
  cout << "generate connectivity vector" << endl ; 
  
  con_vec = new int [n_neurons*n_neurons]() ; 
  
  for(i=0; i<n_neurons; i++) { 
    for(j=0; j<n_neurons; j++) { // Jij -> j (cols, pre) to i (rows, post) 
      if(con_prob[j + i * n_neurons] >= unif(con_gen) ) { 
  	con_vec[j + i * n_neurons] = 1 ; 
  	total_n_post++ ; 
      } 
    } 
  } 
  
  delete [] con_prob ; 
  
  if(IF_SAVE_CON_VEC && !IF_TRIALS) { 
    cout << "saving connectivity vector to: "; 
    write_to_file(con_path, "con_vec", con_vec , n_neurons*n_neurons) ; 
  } 
} 

void func_con_sparse_rep() { 
  
  unsigned long counter = 0 ;  
  
  cout << "generate sparse representation" << endl ;
  cout << "total_n_post " << total_n_post << endl ; 
  
  idx_post = new unsigned long [n_neurons]() ; 
  id_post = (unsigned long *) malloc( (unsigned long) total_n_post * sizeof(unsigned long) ) ; 
  
  n_post = new int [n_neurons]() ;
    
  avg_n_post = new int [n_pop*n_pop]() ;
  
  /* n_pre = new int *[n_pop]() ; */
  /* for(i=0;i<n_pop;i++) // presynaptic population b */
  /*   n_pre[i] = new int [n_neurons]() ; */
  
  /* avg_n_pre = new int [n_pop*n_pop]() ;  */
  
  for(j=0; j<n_neurons; j++) { // Jij -> j (cols, pre) to i (rows, post) 
    pre_pop = which_pop[j] ; 
    
    for(i=0; i<n_neurons; i++) { 
      post_pop = which_pop[i] ;
      
      if(con_vec[j + i * n_neurons]) { // j->i = 1 with proba Kj/Nj 
	id_post[counter] = i ; 
	n_post[j]++ ; // Kb to post pop a, K/N * N = K 
	/* n_pre[pre_pop][i]++ ; // Ka from pre pop a, K/N * Nj = Kj */ 
	avg_n_post[pre_pop + post_pop*n_pop]++ ; 
	counter++ ; 
      } 
    } 
  } 
  
  delete [] con_vec ; 
  
  cout << "average number of postsynaptic connections per pop: " ;
  for(i=0; i<n_pop; i++)
    for(j=0; j<n_pop; j++)
      cout << i << j << " " << avg_n_post[j+i*n_pop] / n_per_pop[i] << " " ;
  cout << endl ;
  
  /* for(i=0; i<2*n_pop; i++) */
  /*   cout << avg_n_post[i] / n_per_pop[i%2] << " " ; */
  /* cout << endl ; */
  
  delete [] avg_n_post ;
  
  /* cout << "number of presynaptic connections per neuron: " << endl ; */

  /* for(i=0; i<n_pop; i++) { //post */
  /*   for(j=0; j<n_pop; j++) { //pre */
  /*     cout << "neuron in " << i << " n_pre in " << j << " | " ; */
  /*     for(k=n_per_pop[i]-1; k>n_per_pop[i]-10; k--)  */
  /* 	cout << n_pre[j][k] << " " ; */
  /*     cout << endl ; */
  /*   } */
  /* } */
  
  /* delete [] n_pre ; */
  
  /* cout << "number of postsynaptic connections per neuron: " << endl ;  */
  
  /* for(i=0; i<n_pop; i++) { //post  */
  /*   for(j=0; j<n_pop; j++) { //pre  */
  /*     cout << "neuron in " << j << " n_post in " << i << " " ;  */
  /*     for(k=cum_n_per_pop[j]; k<cum_n_per_pop[j]+10; k++)  */
  /* 	cout << n_post2[i][k] << " " ;  */
  /*     cout << endl ;  */
  /*   }  */
  /* }  */
   
  cout << "total number of postsynaptic connections per neuron: " << endl ;
  
  for(j=0; j<n_pop; j++) { //pre
    cout << j << " " ;
    for(k=n_per_pop[j]; k>n_per_pop[j]-10; k--)
      cout << n_post[k] << " " ; 
    cout << endl ;
  }
  
  idx_post[0] = 0 ; 
  for(i=1; i<n_neurons; i++) 
    idx_post[i] = idx_post[i-1] + n_post[i-1] ; 
  
}

void check_sparse_rep() { 
  
  cout << "checking sparse representation" << endl ; 
  con_vec = new int [n_neurons*n_neurons]() ; 
  
  for(j=0; j<n_neurons; j++) // Jij -> j (cols, pre) to i (rows, post) 
    for(i=idx_post[j]; i<idx_post[j] + n_post[j]; i++) 
      con_vec[j + id_post[i] * n_neurons] = 1 ; 
  
  cout << "saving rebuild connectivity vector to: "; 
  write_to_file(con_path, "con_vec", con_vec , n_neurons*n_neurons) ; 
  
  delete [] con_vec ;
}

void create_con_dir() { 
  
  con_path += to_string(n_pop)+"pop"; 
  
  if(n_pop==1)
    con_path += "/N" + to_string(n_per_pop[0]/1000) ; 
  else 
    con_path += "/NE_" + to_string(n_per_pop[0]/1000) +  "_NI_" + to_string(n_per_pop[1]/1000) ; 
  
  con_path += "/K" + to_string((int)K) ;
  
  ostringstream str_kappa, str_kappa_var, str_kappa_1 ;
  str_kappa << fixed << setprecision(2) << KAPPA ; 
  str_kappa_var << fixed << setprecision(2) << KAPPA_VAR ; 
  str_kappa_1 << fixed << setprecision(2) << KAPPA_1 ; 
  
  ostringstream str_ksi, str_ksi_var, str_ksi_1 , str_ksi_var_1 ; 
  str_ksi << fixed << setprecision(2) << MEAN_KSI ; 
  str_ksi_var << fixed << setprecision(2) << VAR_KSI ; 
  str_ksi_1 << fixed << setprecision(2) << MEAN_KSI ; 
  str_ksi_var_1 << fixed << setprecision(2) << VAR_KSI ; 
  
  ostringstream str_seed_ksi ; 
  str_seed_ksi << fixed << setprecision(0) << SEED_KSI ; 

  ostringstream str_map_seed ; 
  str_map_seed << fixed << setprecision(0) << MAP_SEED ; 

  ostringstream str_EE, str_EI, str_IE, str_II ;  
  str_EE << fixed << setprecision(0) << SIGMA[0] ; 
  str_EI << fixed << setprecision(0) << SIGMA[1] ; 
  str_IE << fixed << setprecision(0) << SIGMA[2] ; 
  str_II << fixed << setprecision(0) << SIGMA[3] ; 
  
  if(IF_STRUCTURE) {
    if(IF_RING)
      con_path += "/ring/kappa_" + str_kappa.str() ; 
    
    if(IF_SPEC) { 
      if(RANK==1) 
	con_path += "/spec/kappa_" + str_kappa.str() ; 
      if(RANK==2) {
	con_path += "/spec/kappa_" + str_kappa.str() + "_kappa_1_" + str_kappa_1.str() ; 
	if(FIX_MAP_SEED) 
	  con_path += "/seed_" + str_map_seed.str() ; 
      }
        
    }

    if(IF_GAUSS)
      con_path += "/gauss/EE_" + str_EE.str() + "_EI_" + str_EI.str() +"_IE_" + str_IE.str() + "_II_" + str_II.str() ; 
    
    ksi_path = con_path ;
    
    if(IF_LOW_RANK) {
      if(RANK==1){ 
	ksi_path += "/low_rank/rank_1" ; 
	if(FIX_KSI_SEED) 
	  ksi_path += "/seed_ksi_" + str_seed_ksi.str() ; 
	/* con_path += "/low_rank/kappa_" + str_kappa.str() ; */
	/* con_path += "/mean_" + str_ksi.str() + "_var_" + str_ksi_var.str() ; */
      }
      if(RANK==2) { 
	ksi_path += "/low_rank/rank_2" ; 
	if(FIX_KSI_SEED) 
	  ksi_path += "/seed_ksi_" + str_seed_ksi.str() ; 
	/* con_path += "/low_rank/kappa_" + str_kappa.str() + "_kappa_1_" + str_kappa_1.str() ;  */
	/* con_path += "/mean_" + str_ksi.str() + "_var_" + str_ksi_var.str() ;  */
	/* con_path += "/mean_" + str_ksi_1.str() + "_var_" + str_ksi_var_1.str() ;  */
      }
    }
    make_dir(ksi_path) ;
  } 
  
  /* if(IF_INI_COND)  */
  /*   con_path += "/ini_cond_" + to_string( (int) INI_COND_ID ) ;  */
  
  /* if(IF_TRIALS)  */
  /*   con_path += "/trial_" + to_string( (int) TRIAL_ID ) ;  */
  
  make_dir(con_path) ; 
  
  cout << "created directory: " ; 
  cout << con_path << endl ; 
  
}

void save_to_con_file() {

  cout << "save to file" << endl ;
    
  write_to_file(con_path, "id_post", id_post, total_n_post) ; 
  write_to_file(con_path, "idx_post", idx_post, n_neurons) ; 
  write_to_file(con_path, "n_post", n_post, n_neurons) ; 
  
}

#endif
