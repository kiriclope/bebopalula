#ifndef __MATRIXUTILS__ 
#define __MATRIXUTILS__ 

void get_con_sparse_vec() { 

  if(IF_SPEC || IF_RING || IF_GAUSS) 
    init_theta() ; 
  
  create_con_dir() ; 
  
  cout << "###############################################" << endl ; 
  cout << "getting connectivity sparse vectors from:" ; 
  cout << con_path << endl ; 
  
  n_post = new int [n_neurons]() ; 
  read_from_file(con_path, "n_post", n_post, n_neurons) ; 
  
  cout << "average number of postsynaptic connections per neuron: " << endl ; 
  for(i=0; i<n_pop; i++) { 
    cout << "pop " << i << " " ; 
    for(j=cum_n_per_pop[i]; j<cum_n_per_pop[i]+10; j++) 
      cout << n_post[j] << " " ; 
    cout << endl ;     
  }
  
  idx_post = new unsigned long [n_neurons]() ; 
  read_from_file(con_path, "idx_post", idx_post, n_neurons) ; 
  
  cout << "postsynaptic neurons idx: " ; 
  for(i=0; i<10; i++) 
    cout << idx_post[i] << " " ; 
  cout << endl ; 
  
  for(i=0 ; i<n_neurons; i++) 
      total_n_post += n_post[i] ; 
  
  cout << "total number of connections: " << total_n_post << endl ; 
  id_post = (unsigned long *) malloc( (unsigned long) total_n_post * sizeof(unsigned long) ) ;
  /* id_post = new unsigned long [total_n_post]() ; // for some reason this is not working ... */ 
  read_from_file(con_path, "id_post", id_post, total_n_post) ;
  
  cout << "postsynaptic neurons id: " ; 
  for(int i=0; i<10; i++) 
    cout << id_post[i] << " " ; 
  cout << endl ; 
  
  if(IF_LOW_RANK) {
    read_from_file(ksi_path, "ksi", ksi, n_neurons) ; 
    
    cout << "ksi " << endl ;
    for(i=0;i<10;i++)
      cout << ksi[i] << " " ;
    cout << endl ;

    if(IF_SAMPLE)
      read_from_file(ksi_path, "dist", sample, n_neurons) ;
    else
      read_from_file(ksi_path, "sample", sample, n_neurons) ;
    
    read_from_file(ksi_path, "ksi_1", distractor, n_neurons) ; 

    cout << "sample " << endl ;
    for(i=0;i<10;i++)
      cout << sample[i] << " " ;
    cout << endl ;

    cout << "distractor " << endl ;
    for(i=0;i<10;i++)
      cout << distractor[i] << " " ;
    cout << endl ;

    if(RANK==2) {
      read_from_file(ksi_path, "ksi_1", ksi_1, n_neurons) ;
      
      cout << "ksi_1 " << endl ;
      for(i=0;i<10;i++)
  	cout << ksi_1[i] << " " ;
      cout << endl ;
    }
  }
  
  if(IF_SPEC && RANK==2) {
    read_from_file(con_path, "idx_perm", idx_perm, n_neurons) ; 
    read_from_file(con_path, "theta_1", theta_1, n_neurons) ; 
    
    cout << "idx_perm: " ; 
    for(i=0;i<5;i++) 
      cout << idx_perm[i] << " " ; 
    cout << endl ;
    
    cout << "theta_1: " ; 
    for(i=0;i<5;i++)
      cout << theta_1[i] << " " ; 
    cout << endl ;    
  } 
} 

void gen_ksi() {

  cout << "Generate ksi " << endl ; 

  create_con_dir() ; 

  float *array ; 
  random_normal_multivariate(COV_MAT, array, 4, n_neurons) ;
  
  for(i=0; i<10; i++) { 
    for(int j=0; j<4; j++)
      printf("%f ", array[j+i*4] ) ; 
    printf("\n") ; 
  }
  
  for(i=0; i<n_neurons; i++) { 
    sample[i] = array[0+i*4] ; 
    distractor[i] = array[1+i*4] ; 
    ksi[i] = array[2+i*4] ;
    ksi_1[i] = array[3+i*4] ;      
  }

  cout << "sample: " ;
  cout << "mean " << mean_array(sample, n_per_pop[0]) << " var " << var_array(sample, n_per_pop[0]) << endl ; 
  
  cout << "distractor: " ;
  cout << "mean " << mean_array(distractor, n_per_pop[0]) << " var " << var_array(distractor, n_per_pop[0]) << endl ; 

  cout << "covar sample dist " << covar_array(sample, distractor, n_per_pop[0]) << endl ; 
  
  cout << "ksi: " ;
  cout << "mean " << mean_array(ksi, n_per_pop[0]) << " var " << var_array(ksi, n_per_pop[0]) << endl ; 
  cout << "covar sample " << covar_array(ksi, sample, n_per_pop[0]) << " covar dist " << covar_array(ksi, distractor, n_per_pop[0]) ; 
  cout << "covar ksi_1 " << covar_array(ksi, ksi_1, n_per_pop[0]) << endl ; 
  
  /* for(i=0; i<10; i++) */
  /*   cout << ksi[i] << " " ; */
  /* cout << endl ; */
  
  cout << "ksi_1: " ;
  cout << "mean " << mean_array(ksi_1, n_per_pop[0]) << " var " << var_array(ksi_1, n_per_pop[0]) << endl ;
  cout << "covar sample " << covar_array(ksi_1, sample, n_per_pop[0]) << " covar dist " << covar_array(ksi_1, distractor, n_per_pop[0]) ; 
  cout << " covar ksi " << covar_array(ksi_1, ksi, n_per_pop[0]) << endl ;
  
  /* for(i=0; i<10; i++) */
  /*   cout << ksi_1[i] << " " ; */
  /* cout << endl ; */
  
  /* double cos_ksi ;  */
  /* cos_ksi = cos_array(sample, distractor, n_per_pop[0]) ;  */

  /* cout << " angle_sample_dist: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ;  */

  /* cos_ksi = cos_array(ksi, ksi_1, n_per_pop[0]) ;  */
  
  /* // cout << "cos_ksi_ksi_1: " << cos_ksi ;  */
  /* cout << " angle_ksi_ksi_1: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ;  */
  
  /* cos_ksi = cos_array(ksi, sample, n_per_pop[0]) ;  */
  
  /* // cout << "cos_ksi_sample: " << cos_ksi ;  */
  /* cout << " angle_ksi_sample: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ;  */

  /* cos_ksi = cos_array(ksi, distractor, n_per_pop[0]) ;  */

  /* // cout << "cos_ksi_dist: " << cos_ksi ;  */
  /* cout << " angle_ksi_dist: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ;  */

  /* cos_ksi = cos_array(ksi_1, sample, n_per_pop[0]) ;  */

  /* // cout << "cos_ksi_sample: " << cos_ksi ;  */
  /* cout << " angle_ksi1_sample: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ;  */
  
  /* cos_ksi = cos_array(ksi_1, distractor, n_per_pop[0]) ;  */

  /* // cout << "cos_ksi_dist: " << cos_ksi ;  */
  /* cout << " angle_ksi1_dist: " << acos(cos_ksi) * 180.0 / M_PI << "°" << endl ;  */
  
  // init_ksi() ; 
  cout << "###############################################" << endl ;
  
  write_to_file(ksi_path, "ksi", ksi , n_neurons) ; 
  // write_to_file(con_path, "ksi_scaled", ksi_scaled , n_per_pop[0]) ; 
  write_to_file(ksi_path, "sample", sample , n_neurons) ; 
  write_to_file(ksi_path, "dist", distractor , n_neurons) ; 
    
  if(RANK==2) {
    // init_ksi_1() ; 
    write_to_file(ksi_path, "ksi_1", ksi_1 , n_neurons) ; 
    // write_to_file(con_path, "ksi_1_scaled", ksi_1_scaled , n_per_pop[0]) ; 
  }
  
  /* init_ksi() ;  */
  /* write_to_file(ksi_path, "ksi", ksi , n_neurons) ;  */
  /* // write_to_file(con_path, "ksi_scaled", ksi_scaled , n_per_pop[0]) ;  */
    
  /* if(RANK==2) { */
  /*   init_ksi_1() ;  */
  /*   write_to_file(ksi_path, "ksi_1", ksi_1 , n_neurons) ;  */
  /*   // write_to_file(con_path, "ksi_1_scaled", ksi_1_scaled , n_per_pop[0]) ;  */
  /* } */
}

void gen_con_sparse_vec() { 
  cout << "###############################################" << endl ;
  
  /* if( (IF_SAVE_CON_VEC || IF_SAVE_SPARSE_REP) && !IF_TRIALS)  */
  create_con_dir() ; 

  /* if(IF_GEN_KSI) */
  /*   gen_ksi() ;     */
  /* else */
  if(IF_LOW_RANK) {
    cout << ksi_path << endl ; 
    read_from_file(ksi_path, "ksi", ksi, n_neurons) ; 
    
    if(IF_SAMPLE)
      read_from_file(ksi_path, "dist", sample, n_neurons) ;
    else
      read_from_file(ksi_path, "sample", sample, n_neurons) ;
    
    read_from_file(ksi_path, "ksi_1", distractor, n_neurons) ; 
    
    cout << "ksi " << endl ;
    for(i=0;i<10;i++)
      cout << ksi[i] << " " ;
    cout << endl ;
    
    if(RANK==2) {
      read_from_file(ksi_path, "ksi_1", ksi_1, n_neurons) ;
      
      cout << "ksi_1 " << endl ;
      for(i=0;i<10;i++)
	cout << ksi_1[i] << " " ;
      cout << endl ;
    }
  }
  
  func_con_prob() ; 
  func_con_vec() ; 
  
  func_con_sparse_rep() ;
  
  if(IF_CHECK_SPARSE_REP) 
    check_sparse_rep() ; 
  
  if(IF_SPEC) 
    if(RANK==2) { 
      write_to_file(con_path, "idx_perm", idx_perm, n_neurons) ; 
      write_to_file(con_path, "theta_1", theta_1, n_neurons) ; 
    }
    
  /* if(IF_LOW_RANK) { */
  /*   write_to_file(ksi_path, "ksi", ksi , n_neurons) ;  */
  /*   if(RANK==2)  */
  /*     write_to_file(ksi_path, "ksi_1", ksi_1 , n_neurons) ;  */
  /* } */
  
  if(IF_SAVE_SPARSE_REP && !IF_TRIALS) 
    save_to_con_file() ; 
  
}

#endif
