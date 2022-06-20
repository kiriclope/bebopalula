#ifndef __STPUTILS__ 
#define __STPUTILS__ 

//STP Parameters 

double *u_stp ; // availability variable 
double *x_stp, *x_stp_old ; // resource variable 
double *A_u_x_stp ; // A*u*x variable 

double *u_stp_mongillo, *x_stp_mongillo ; // availability variable 

double *unbind_rate, *refill_rate, *bind_rate ;

string str_u_stp ; 
ofstream file_u_stp ; 

string str_x_stp ; 
ofstream file_x_stp ; 

const int stp_synapse[4] = {1, 0, 0, 0} ; 

void init_stp_globals(){ 
  u_stp = new double [n_neurons]() ; 
  x_stp = new double [n_neurons]() ;

  u_stp_mongillo = new double [n_neurons*n_neurons]() ; 
  x_stp_mongillo = new double [n_neurons*n_neurons]() ; 
  
  x_stp_old = new double [n_neurons]() ; 
  A_u_x_stp = new double [n_neurons]() ;
  
  unbind_rate = new double [n_neurons]() ; 
  refill_rate = new double [n_neurons]() ; 
  bind_rate = new double [n_neurons]() ; 
  
}

void delete_stp_globals(){ 
  delete [] u_stp ; 
  delete [] x_stp ;
  
  delete [] x_stp_old ; 

  delete [] u_stp_mongillo ; 
  delete [] x_stp_mongillo ;
  
  delete [] unbind_rate ;
  delete [] refill_rate ;
  delete [] bind_rate ; 
}

void open_stp_files(){ 
  str_u_stp = path + "/u_stp.dat" ; 
  file_u_stp.open(str_u_stp.c_str(), ios::out | ios::ate) ;
  
  str_x_stp = path + "/x_stp.dat" ; 
  file_x_stp.open(str_x_stp.c_str(), ios::out | ios::ate) ; 
}

void close_stp_files(){
  file_u_stp.close() ; 
  file_x_stp.close() ; 
}

void markram() {
  
  x_stp[i_neuron] -= u_stp[i_neuron] * x_stp[i_neuron] ; 
  u_stp[i_neuron] += USE * (1.0 - u_stp[i_neuron]) ; 
  A_u_x_stp[i_neuron] = u_stp[i_neuron] * x_stp[i_neuron] ;
  
}

/* Markram98, Mato  */
void mato() { 
  u_stp[i_neuron] = u_stp[i_neuron] * exp(- ISI[i_neuron] / TAU_FAC ) + USE * ( 1. - u_stp[i_neuron] * exp(- ISI[i_neuron] / TAU_FAC ) ) ; 
  x_stp[i_neuron] = x_stp[i_neuron] * ( 1. - u_stp[i_neuron] ) * exp(- ISI[i_neuron] / TAU_REC ) + 1. - exp(- ISI[i_neuron] / TAU_REC ) ; 
  A_u_x_stp[i_neuron] = u_stp[i_neuron] * x_stp[i_neuron] ; 
}

void release_stp() {
  
  if(IF_MARKRAM) { 
    u_stp[i_neuron] *= EXP_DT_TAU_FAC ; 
    /* x_stp_old[i_neuron] = x_stp[i_neuron] ;  */
    x_stp[i_neuron] += DT_OVER_TAU_REC * (1.0 - x_stp[i_neuron]) ; 
  }
  
  if(IF_MONGILLO==1) {
    if( u_stp[i_neuron] == 1.0 && unif(rand_gen) < 1.0 / TAU_FAC) {// u = 1 -> 0 , calcium unbinds with rate 1/Tfac 
      u_stp[i_neuron] = 0.0 ; 
      unbind_rate[i_neuron] += 1.0 ; 
    }
    if( x_stp[i_neuron] == 0.0 && unif(rand_gen) < 1.0 / TAU_REC) {// x = 0 -> 1 , neurotransmitter refills with rate 1/Trec 
      x_stp[i_neuron] = 1.0 ;
      refill_rate[i_neuron] += 1.0 ; 
    }
  }
}

void mongillo() { // mongillo 
  
  // Updating STP variables
  if( u_stp[i_neuron]==0.0 && unif(rand_gen) < USE) { 
    u_stp[i_neuron] = 1.0 ; // calcium binds with probability Use u = 0->1
    /* bind_rate[i_neuron] += 1.0 ;  */
    
    if(x_stp[i_neuron]==1.0) {
      A_u_x_stp[i_neuron] = 1.0 ; x_stp[i_neuron] = 0.0 ;
    } // neurotransmitter release if available x = 1->0    
    else
      A_u_x_stp[i_neuron] = 0.0 ;
  }
  else
    A_u_x_stp[i_neuron] = 0.0 ;
}

void mongillo_alt() { // mongillo
  
  u_stp[i_neuron] = u_stp_mongillo[i_neuron+id_post[j]*n_neurons] ; 
  x_stp[i_neuron] = x_stp_mongillo[i_neuron+id_post[j]*n_neurons] ; 
  
  if (u_stp[i_neuron]==1.0)
    if ( unif(rand_gen)<(1.0-exp(-ISI[i_neuron]/TAU_FAC)))
      u_stp[i_neuron]=0.0;
  
  if (u_stp[i_neuron]==0.0)
    if (unif(rand_gen)<USE)
      u_stp[i_neuron]=1.0;
  
  if (x_stp[i_neuron]==0.0)
    if ( unif(rand_gen)<(1.0-exp(-ISI[i_neuron]/TAU_REC)))
      x_stp[i_neuron]=1.0 ; 
  
  A_u_x_stp[i_neuron] = u_stp[i_neuron] * x_stp[i_neuron] ;

  if (u_stp[i_neuron]==1.0)
    u_stp_mongillo[i_neuron+id_post[j]*n_neurons] = 1.0 ;
  else
    x_stp_mongillo[i_neuron+id_post[j]*n_neurons] = x_stp[i_neuron] ; 
}

void update_stp_variables_lif() {
  
  if(IF_MARKRAM)
    markram() ; 
  if(IF_MATO)
    mato() ; 
  if(IF_MONGILLO==1)
    mongillo() ;
  if(IF_MONGILLO==2)
    mongillo_alt() ; 
  
}

void save_xy_to_file() { 
  
  /* cout << "unbind_rate" ;  */
  /* for(j=0;j<5;j++) {  */
  /*   cout << unbind_rate[j] *1000./TIME_WINDOW << " " ;  */
  /*   unbind_rate[j] = 0 ;  */
  /* } */
  /* cout << endl ;  */
  
  /* cout << "refill rate" ; */
  /* for(j=0;j<5;j++) { */
  /*   cout << refill_rate[j] *1000./TIME_WINDOW << " " ; */
  /*   refill_rate[j] = 0 ;     */
  /* } */
  /* cout << endl ;  */
  
  /* cout << "bind rate" ;  */
  /* for(j=0;j<5;j++) {  */
  /*   cout << bind_rate[j] *1000./TIME_WINDOW << " " ;  */
  /*   refill_rate[j] = 0 ;  */
  /* } */
  /* cout << endl ;  */

  file_u_stp << t_time-TIME_STEADY ; 
  for(j=0;j<n_neurons;j++) 
    file_u_stp << " " << u_stp[j] ;  
  file_u_stp << endl ; 
  
  file_x_stp << t_time-TIME_STEADY ; 
  for(j=0;j<n_neurons;j++) 
    file_x_stp << " " << x_stp[j] ;  
  file_x_stp << endl ; 
}


#endif
