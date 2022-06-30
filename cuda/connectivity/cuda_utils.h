#ifndef __CUDA_UTILS__
#define __CUDA_UTILS__

__global__ void setup_kernel() { 
  unsigned long id = threadIdx.x + blockIdx.x * blockDim.x ; 
  unsigned long i_neuron = (unsigned long) ( id + dev_i_chunck * N_NEURONS_PER_CHUNCK ) ;
  /* Each thread gets different seed, a different sequence number, no offset */ 
  /* if(id < N_NEURONS_PER_CHUNCK) */
  /*   curand_init(clock64(), id, 0, &dev_states[id]) ; */
  if(id < N_NEURONS_PER_CHUNCK & i_neuron < N_NEURONS) 
    curand_init(clock64(), id, 0, &dev_states[i_neuron]) ; 
  
  /* if(id<2) {  */
  /*   cuPrintf("currand init: ") ;  */
  /*   cuPrintf("%d", dev_states[id]) ;  */
  /* }   */
}

__device__ double unif_dist(unsigned long i_neuron) { 
  /*RETURNS ONE SAMPLE FROM UNIFORM DISTRIBUTION*/ 
  double randNumber= 0.0 ; 
  if(i_neuron < N_NEURONS) { 
    /* save state in local memory for efficiency */ 
    curandState localState = dev_states[i_neuron] ;
    randNumber = curand_uniform(&localState) ; 
    dev_states[i_neuron] = localState ; 
  } 
  return randNumber ; 
}

__device__ void init_dev_theta() { 
  
  for(int i=0; i<n_pop; i++) 
    for(unsigned long j = dev_cum_n_per_pop[i]; j<dev_cum_n_per_pop[i+1]; j++) 
      dev_theta[j] = 2.0 * (double) M_PI * (double) (j-dev_cum_n_per_pop[i]) / (double) dev_n_per_pop[i] ; 
  
  /* unsigned long id = (unsigned long) threadIdx.x + blockIdx.x * blockDim.x ;  */
  
  /* if(id<1) {  */
  /*   cuPrintf("theta:") ; */
  /*   for(int i=0; i<5; i++)  */
  /*     cuPrintf("%.2f", dev_theta[i]) ;  */
  /*   cuPrintf("\n") ;  */
  /* } */
  
} 

__global__ void kernel_gen_con_prob() { 
  
  unsigned long id = (unsigned long) threadIdx.x + blockIdx.x * blockDim.x ; 
  unsigned long i_neuron = (unsigned long) ( id + dev_i_chunck * N_NEURONS_PER_CHUNCK ) ; 

  double kappa = 2.0 * KAPPA_E ; 

  DEV_IF_STRUCTURE = IF_SPEC || IF_RING ; 
  
  init_dev_theta() ; 
  
  if(IF_SPEC) 
    kappa /= sqrt_K ; 
  
  if(id < N_NEURONS_PER_CHUNCK & i_neuron < N_NEURONS) {     
    dev_pre_pop = dev_which_pop[i_neuron] ; 

    if(dev_pre_pop==1) 
	kappa = 2.0 * KAPPA_I ; 
        
    for(unsigned long i=0; i<N_NEURONS; i++) { // id (pre) -> i (post)       
      dev_post_pop = dev_which_pop[i] ; 
      /* dev_con_prob_chunck[id + i * N_NEURONS_PER_CHUNCK ] = dev_K_over_Na[dev_pre_pop] ; */
      dev_con_prob_chunck[i + id * N_NEURONS] = dev_K_over_Na[dev_pre_pop] ;
      
      if(DEV_IF_STRUCTURE)
	if(DEV_IS_STRUCT_SYN[dev_pre_pop + dev_post_pop * n_pop]) 
	  /* dev_con_prob_chunck[id + i * N_NEURONS_PER_CHUNCK ] *= ( 1.0 + kappa * cos( dev_theta[i_neuron] - dev_theta[i] ) ) ; */
	  dev_con_prob_chunck[i + id * N_NEURONS] *= ( 1.0 + kappa * cos( dev_theta[i_neuron] - dev_theta[i] ) ) ; 
      
    } 
  } 
} 


__global__ void kernel_gen_con_vec() { 

  unsigned long id = (unsigned long) threadIdx.x + blockIdx.x * blockDim.x ; 
  unsigned long i_neuron = (unsigned long) ( id + dev_i_chunck * N_NEURONS_PER_CHUNCK ) ; 
    
  /* if(id<1) {  */
  /*   cuPrintf("DEV_IF_STRUCTURE %d \n", DEV_IF_STRUCTURE) ;     */ 
  /*   cuPrintf("dev_K_over_Na:") ;  */
  /*   for(int i=0; i<n_pop; i++)  */
  /*     cuPrintf("%.2f", dev_K_over_Na[i]) ; */
  /*   cuPrintf("\n") ; */
  /* } */

  /* if(id<2) {    */
  /*   cuPrintf("check random seed: ") ; */
    
  /*   for(int i=0; i<5; i++) */
  /*     cuPrintf("%.2f", unif_dist(i) ) ;  */
  /*   cuPrintf("\n") ; */
    
  /* } */
  
  if(id < N_NEURONS_PER_CHUNCK && i_neuron < N_NEURONS) // presynaptic neuron (columns) 
    for(unsigned long i=0; i<N_NEURONS; i++) { // post
      /* if( dev_con_prob_chunck[id + i * N_NEURONS_PER_CHUNCK] >= unif_dist(id) ) { // WARNING must be id inside unif otherwise problems */
      /* 	dev_con_vec_chunck[id + i * N_NEURONS_PER_CHUNCK] = 1 ;  */
      /* 	dev_total_n_post++ ;  */
      /* 	dev_n_post[i_neuron]++ ;  */
      /* } */

      if( dev_con_prob_chunck[i + id * N_NEURONS] >= unif_dist(id) ) { 
	dev_con_vec_chunck[i + id * N_NEURONS] = 1 ; 
	dev_total_n_post++ ; 
	dev_n_post[i_neuron]++ ; 
      } 
      
      /* else  */ 
      /* 	dev_con_vec_chunck[id + i * N_NEURONS_PER_CHUNCK] = 0 ;  */
    }  
  
  /* if(id<1) {  */
  /*   cuPrintf("dev_n_post: ") ; */
  /*   for(int i=0; i<5; i++)  */
  /*     cuPrintf("%d ", dev_n_post[i]) ; */
  /*   cuPrintf("\n") ;  */
  /* }  */
} 

#endif 
