#ifndef __GLOBALS__
#define __GLOBALS__

#include<iostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace:: std ; 

#define n_pop 2 
#define N_NEURONS (unsigned long) 40000 
#define K (double) 2000 
#define sqrt_K (double) sqrt(K) 

#define E_frac (double) 0.8

const double n_frac[2] = { E_frac, round( (1.0 - E_frac)*100.0) / 100.0 } ; 

// string path = "/homecentral/alexandre.mahrach/IDIBAPS/connectivity/" ;
string path = "../../cpp/model/connectivity/" ;

#define IF_CON_DIR 1 
#define CUE (double) 0.35 

unsigned long i, j, i_neuron ;
int pre_pop, post_pop ; 
int n_pref ; 

unsigned long n_per_pop[n_pop] ; 
unsigned long cum_n_per_pop[n_pop+1] ; 

int which_pop[N_NEURONS] ; 
double K_over_Na[n_pop] ; 

#define IF_SAVE_CON_VEC 0 // save Cij matrix 
#define IF_SAVE_SPARSE_REP 1 

////////////////////////////////// 
// structure 
////////////////////////////////// 

int IF_STRUCTURE ; 
const double IS_STRUCT_SYN[4] = {1.0, 1.0, 1.0, 1.0} ; // WARNING check that it is the same in cuda globals

#define IF_RING 1 
#define IF_SPEC 0 
#define IF_GAUSS 0 

#define KAPPA (double) 0.25
#define KAPPA_E (double) 0.25
#define KAPPA_I (double) 0.125

#define SIGMA (double) 1.0

#define IF_LOW_RANK 0 
#define MEAN_XI 0.0 
#define VAR_XI 1.0 

#endif
