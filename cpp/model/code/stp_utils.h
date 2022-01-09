#ifndef __STPUTILS__ 
#define __STPUTILS__ 

//STP Parameters 
#define A_STP (double) 1.0 

double *u_stp ; // availability variable 
double *x_stp ; // resource variable 
double *A_u_x_stp ; // A*u*x variable 

string str_u_stp ; 
ofstream file_u_stp ; 

string str_x_stp ; 
ofstream file_x_stp ; 

const int stp_synapse[4] = {1, 0, 0, 0} ; 

void init_stp_globals(){ 
  u_stp = new double [n_neurons]() ; 
  x_stp = new double [n_neurons]() ; 
  A_u_x_stp = new double [n_neurons]() ; 
}

void delete_stp_globals(){ 
  delete [] u_stp ; 
  delete [] x_stp ; 
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

// Markram98, Mato 
void update_stp_variables_lif(double ISI) { 
  u_stp[i_neuron] = u_stp[i_neuron] * exp(- ISI / TAU_FAC ) + USE * ( 1. - u_stp[i_neuron] * exp(- ISI / TAU_FAC ) ) ; 
  
  x_stp[i_neuron] = x_stp[i_neuron] * ( 1. - u_stp[i_neuron] ) * exp(- ISI / TAU_REC ) + 1. - exp(- ISI / TAU_REC ) ; 
  
  /* A_u_x_stp[i_neuron] = A_STP * u_stp[i_neuron] * x_stp[i_neuron] ;  */
  A_u_x_stp[i_neuron] = u_stp[i_neuron] * x_stp[i_neuron] ; 
} 

/* void mongillo() { */
  
/*   if( u_stp[i_neuron] = 1.0 && unif(rand_gen) < DT / TAU_FAC ) // u = 1 -> 0 , calcium unbinds with rate 1/Tfac  */
/*     u_stp[i_neuron] = 0.0 ;  */
/*   if( x_stp[i_neuron] == 0.0 && unif(rand_gen) < DT / TAU_REC ) // x = 0 -> 1 , neurotransmitter refills with rate 1/Trec  */
/*     x_stp[i_neuron] = 1.0 ;  */
  
/*   /\* // Updating STP variables *\/ */
/*   /\* if( u_stp[i_neuron]==0.0 && unif(rand_gen) < DT * USE ) { u_stp[i_neuron] = 1.0 ; // calcium binds with probability Use u = 0->1  *\/ */
/*   /\*   if(x_stp[i_neuron]==1.0) { Jstp[i_neuron] = 1.0 ; x[i_neuron] = 0.0 ; } ; // neurotransmitter release if available x = 1->0  *\/ */
/*   /\* } *\/ */
  
/*   A_u_x_stp[i_neuron] = u_stp[i_neuron] * x_stp[i_neuron] ;  */
  
/* } */


void save_xy_to_file() { 
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
