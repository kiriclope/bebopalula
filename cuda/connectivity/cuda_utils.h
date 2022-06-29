#ifndef __CUDA_UTILS__
#define __CUDA_UTILS__

__global__ void setup_kernel() { 
  unsigned long id = threadIdx.x + blockIdx.x * blockDim.x ; 
  /* Each thread gets different seed, a different sequence number, no offset */ 
  if(id < N_NEURONS) 
    curand_init(clock64(), id, 0, &dev_states[id]) ; 

  if(id<2) { 
    cuPrintf("currand init: ") ; 
    cuPrintf("%d", dev_states[id]) ; 
  }  
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

  double kappa = KAPPA_E ;
  
  init_dev_theta() ;
  
  if(IF_SPEC) 
    if(id < N_NEURONS_PER_CHUNCK & i_neuron < N_NEURONS) { 
      
      dev_pre_pop = dev_which_pop[i_neuron] ; 
      
      for(unsigned long j=0; j<N_NEURONS; j++) { 
	
	dev_post_pop = dev_which_pop[j] ; 
	
	dev_con_prob_chunck[j + id * N_NEURONS ] = dev_K_over_Na[dev_post_pop] 
	  * ( 1.0 + DEV_IS_STRUCT_SYN[dev_pre_pop + dev_post_pop * n_pop] / sqrt_K * KAPPA * cos( dev_theta[i_neuron] - dev_theta[j] ) ) ; 
      } 
    }
  
  if(IF_RING) 
    if(id < N_NEURONS_PER_CHUNCK & i_neuron < N_NEURONS) { 
    
      dev_pre_pop = dev_which_pop[i_neuron] ;
      
      if(dev_pre_pop==1)
	kappa = KAPPA_I ;
	
      for(unsigned long j=0; j<N_NEURONS; j++) { 
	
	dev_post_pop = dev_which_pop[j] ; 
	
	dev_con_prob_chunck[j + id * N_NEURONS ] = dev_K_over_Na[dev_post_pop] 
	  * ( 1.0 + DEV_IS_STRUCT_SYN[dev_pre_pop + dev_post_pop * n_pop] * kappa * cos( dev_theta[i_neuron] - dev_theta[j] ) ) ; 
      } 
    } 
} 


__global__ void kernel_gen_con_vec() { 

  unsigned long id = (unsigned long) threadIdx.x + blockIdx.x * blockDim.x ; 
  unsigned long i_neuron = (unsigned long) ( id + dev_i_chunck * N_NEURONS_PER_CHUNCK ) ; 
  
  DEV_IF_STRUCTURE = IF_SPEC || IF_RING ; 
  
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
  
  if(DEV_IF_STRUCTURE) { 
    if(id < N_NEURONS_PER_CHUNCK && i_neuron < N_NEURONS) 
      for(unsigned long i=0; i<N_NEURONS; i++) { 
	if( dev_con_prob_chunck[i + id * N_NEURONS] >= unif_dist(i_neuron) ) { 
	  dev_con_vec_chunck[i + id * N_NEURONS] = 1 ; 
	  dev_total_n_post++ ; 
	  dev_n_post[i_neuron]++ ; 
	} 
	else 
	  dev_con_vec_chunck[i + id * N_NEURONS] = 0 ; 
      } 
  } 
  else { 
    if(id < N_NEURONS_PER_CHUNCK && i_neuron < N_NEURONS) { // presynaptic neuron (columns) 
      
      dev_pre_pop = dev_which_pop[i_neuron] ; 
      
      for(unsigned long i=0; i<N_NEURONS; i++) { // candidates for postsynaptic neurons (rows) 
	
	dev_post_pop = dev_which_pop[i] ; 
	
	if( dev_K_over_Na[dev_post_pop] >= unif_dist(i_neuron) ) { 
	  dev_con_vec_chunck[i + id * N_NEURONS] = 1 ; 
	  dev_total_n_post++ ; 
	  dev_n_post[i_neuron]++ ;
	} 
	else 
	  dev_con_vec_chunck[i + id * N_NEURONS] = 0 ; 
      } //endfor 
    } //endif 
    
  } //endelse

  if(id<1) { 
    cuPrintf("dev_n_post: ") ;
    for(int i=0; i<5; i++) 
      cuPrintf("%d ", dev_n_post[i]) ;
    cuPrintf("\n") ; 
  }   
} 

#endif 
