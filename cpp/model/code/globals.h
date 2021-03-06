#ifndef __GLOBALVARS__ 
#define __GLOBALVARS__ 

using namespace:: std ; 

random_device rd ; 
mt19937 rand_gen(rd()) ;

uniform_real_distribution<double> unif(0.0, 1.0) ; 
normal_distribution<double> white_noise(0.0, 1.0) ;  

#define IF_BIN 0 
#define IF_LIF 1 
#define IF_RATE 0 

#define IF_SINGLE_NEURON 0 

#define E_frac (double) 0.8 
const double n_frac[2] = { E_frac, round( (1.0 - E_frac)*100.0) / 100.0 } ; 

#define IF_MATO_K 1 
#define IF_RESCALE 0 

////////////////////////////////////////// 
// Simulation globals 
////////////////////////////////////////// 
#define DT (double) .1 

#define TIME_INI (double) 0.000E3 // 1.00E3 // 10.E3 // 
#define TIME_STEADY (double) 10.000E3 //  2.0E3 // 10.E3 //

#define DURATION (double) 30.00E3 // 10.E3 // 

#define TIME_WINDOW (double) 0.2500E3 // 1.00E3 // 10.E3 // 
#define TIME_REC (double) 60.00E3 // 1.00E3 // 10.E3 // 
#define TIME_REC_SPIKES (double) 10.00E3 // 1.00E3 // 10.E3 // 

#define IF_RK2 0 

#define N_PREF 10000 

#define IF_TRIALS 0 
#define TRIAL_ID 0 

#define IF_INI_COND 0 
#define INI_COND_ID 0 

#define IF_SAVE_VOLT 0 
////////////////////////////////////////// 
// Network globals 
////////////////////////////////////////// 
const double eps[2] = {1.0,-1.0} ; 
const double Trate[2] = {10.0, 10.0} ; 

#define J0 -1.0 // for one population 
#define I0 1.0 // for one population 
#define Tsyn0 1.0 // for one population 

#define GAIN (double) 1.0 
#define M0 (double) 0.001

// 0.05 if with 0.01 without stp 
#define IF_LOOP_M0 0 
#define IF_LOOP_GAIN 0 

double m0, m1[2], phase[2] ; 
double duration, time_rec, time_rec_spikes, time_steady ; 

string dir ; 
string path = "../" ; 

unsigned long i, j, k, l, i_neuron ; 
double t_time, percentage, t_window=0. ; 
int pre_pop, post_pop ; 

int n_pop ; 
unsigned long n_neurons ; 
double K, sqrt_K, *sqrt_Ka, *Ka ; 

double *ext_inputs, *ff_inputs, *J, *tau_syn , *J_scaled, *J_nmda, *ext_inputs_scaled ; 

unsigned long *n_per_pop, *cum_n_per_pop ; 
int *which_pop ; 

double *mf_rates ; 

double *volt, vold ;
double RK1, RK2 ; 
double *t_spike ; 

double *mean_rates ; 
double *filter_rates ; 

double **inputs, **inputs_nmda ; 
double **filter_inputs ; 

double *net_inputs ; 
double *net_inputs_RK2 ; 

double *overlaps ; 

////////////////////////////////////////// 
// Connectivity globals 
//////////////////////////////////////////

#define IF_GEN_CON 0
#define IF_GEN_KSI 0
#define IF_SAVE_CON_VEC 0 
#define IF_SAVE_SPARSE_REP 0 
#define IF_CHECK_SPARSE_REP 0 

#define IF_RING 1
#define IF_SPEC 0 

#define KAPPA (double) .25 // 4 // 3.5 // 12.0 14 // rank 1: 4.0 
#define KAPPA_1 (double) 4 // 3.5 // 8.0 12 
#define KAPPA_VAR (double) 0.0

#define IF_KAPPA_DIST 0 
lognormal_distribution<double> log_normal(KAPPA, KAPPA_VAR) ; 

const double kappas[2] = {KAPPA, KAPPA*.5} ; 

#define KAPPA_CORR (double) 0 

#define IF_GAUSS 0 
/* #define L 2.0*M_PI  */
#define DEG_TO_RAD (double) M_PI/180.0 
const double SIGMA[4] = {60.0, 60.0, 70.0, 60.0} ; 

#define SIGMA_FF 100.0 

double *X ; 

#define PHASE M_PI/2.0 
#define N_PHASES 0 
#define RANK 1 
#define MAP_ANGLE (double) M_PI/4.0 
// in radians 
#define FIX_MAP_SEED 1 
#define MAP_SEED 4 

#define IF_LOW_RANK 0
#define FIX_KSI_SEED 1 
#define SEED_KSI (double) 3

#define SEED_CON (double) 1 // rd()

mt19937 con_gen(exp(SEED_CON)) ; 
mt19937 ksi_gen(exp(SEED_KSI)) ; 
mt19937 ksi_1_gen(sqrt(SEED_KSI)) ; 
mt19937 covar_ksi_gen(10.0) ; 

#define IF_NORM_KSI 0 
#define MEAN_KSI (double) 0.0 
#define VAR_KSI (double) 1.0 

#define IF_NORM_KSI_1 0 
#define MEAN_KSI_1 (double) 0.0 
#define VAR_KSI_1 (double) 1.0 

#define VAR_SAMPLE (double) 1.0
#define VAR_DIST (double) 1.0

double COV_MAT[16] = { // sample, dist, ksi, ksi1
  1.0, 0.0, 0.6, 0.2, 
  0.0, 1.0, -0.6, 0.2, 
  0.6, -0.6, 1.0, -0.1/2000.,
  0.2, 0.2, -0.1/2000., 1.0, 
} ;

#define ANGLE_KSI (double) M_PI/4.0 

double *ksi, *ksi_1, *ksi_init, *ksi_scaled, *ksi_1_scaled ;
double *sample, *distractor ;

#define IF_SAMPLE 0

string con_path ; 
string ksi_path ; 

int *n_post, *avg_n_post, *con_vec, **n_pre, *avg_n_pre ; 
unsigned long *id_post, *idx_post ; 
unsigned long *idx_perm, *idx_perm_E, *idx_perm_I ; 

unsigned long total_n_post = 0 ; 

double *K_over_Na ; 
double *con_prob ; 
double *prefactor ; 

const double IS_STRUCT_SYN[4] = {1.0, 0.0, 0.0, 0.0} ; 
double *theta, *theta_1 ; 

int IF_STRUCTURE ; 

//////////////////////////////////////
// Stimulus / Task globals 
//////////////////////////////////////

int IF_STIM ;

int SWITCH_ON = 0 ; 
int SWITCH_OFF = 0 ;

#define IF_TRACK 0 

#define KAPPA_EXT (double) 1.0 
#define PHI_EXT (double) 0.25 

#define KAPPA_DIST (double) 0.5
#define PHI_DIST (double) M_PI 

#define KAPPA_CUE (double) 0.25 
#define KAPPA_TEST (double) 2.0

#define IF_TUNED_FF 0

#define IF_STEP 1
#define T_STEP_ON (double) 2000.0
#define T_STEP_OFF (double) 3000.0
const double A_STEP[2] = {1.55, 1.25} ; 

//////////////////////////////////////
// Christos
//////////////////////////////////////

#define IF_CHRISTOS 0
#define T_CUE_ON (double) 2000 
#define T_CUE_OFF (double) 3000 

#define T_ERASE_ON (double) 600000
#define T_ERASE_OFF (double) 700000 

const double A_CUE[2] = {0.25, 0.0} ; // {2.4, 1.0} 
const double EPS_CUE[2] =  {0.17, 0.0} ; // {.17 , 0.0} 

const double A_ERASE[2] = {5.2, 3.7} ; 
const double EPS_ERASE[2] = {0.23, 0.28} ; 

//////////////////////////////////////
// DUAL TASK
//////////////////////////////////////

#define IF_DPA 0
#define IF_DUAL 0 
#define IF_DRT 0 

#define T_SAMPLE_ON (double) 2000 
#define T_SAMPLE_OFF (double) 3000 

#define T_DIST_ON (double) 4500 
#define T_DIST_OFF (double) 5500 

#define T_RWD_ON (double) 6500 
#define T_RWD_OFF (double) 7500 

#define T_TEST_ON (double) 9000 
#define T_TEST_OFF (double) 10000 

//////////////////////////////////////
// LIF globals 
//////////////////////////////////////

#define IF_SYN_DYN 1 // 1 exponential synapses, 0 delta synapses 

// mongillo 
/* #define Vl 0.75 // -70. // -55. // -70. // Membrane leak Potential  */
/* #define Vr 0.75 //-70. // -55. // -70. // Membrane leak Potential  */
/* #define Vth 1.0 //-50. // -40. // -50. // Voltage Threshold  */
/* #define Vpeak 1.2 // 20. // 20. // Spikes Peak  */

/* #define Vl -70. // -55. // -70. // Membrane leak Potential */
/* #define Vr -70. // -55. // -70. // Membrane leak Potential */
/* #define Vth -50. // -40. // -50. // Voltage Threshold */
/* #define Vpeak 20. // 20. // Spikes Peak */

/* /\* // mato *\/ */
#define Vl 0. // -55. // -70. // Membrane leak Potential
#define Vr -3.33 // -55. // -70. // Membrane leak Potential
#define Vth 20. // -40. // -50. // Voltage Threshold
#define Vpeak 20. // 20. // Spikes Peak

/* #define Vl 0. // -55. // -70. // Membrane leak Potential */
/* #define Vr -5.0 // -55. // -70. // Membrane leak Potential  */
/* #define Vth 30.0 // -40. // -50. // Voltage Threshold */
/* #define Vpeak 20. // 20. // Spikes Peak */

const double TAU_MEM[2] = {20.0, 10.0} ; 
const double dt_over_tau_mem[2] = { DT/TAU_MEM[0], DT/TAU_MEM[1] } ; 
const double one_minus_dt_over_tau_mem[2]={ 1.0-dt_over_tau_mem[0], 1.0-dt_over_tau_mem[1] } ; 
const double one_plus_dt_over_tau_mem[2]={ 1.0+dt_over_tau_mem[0], 1.0+dt_over_tau_mem[1] } ; 
const double dt_over_dt_tau_mem[2]={ DT/(1.0+dt_over_tau_mem[0]), DT/(1.0+dt_over_tau_mem[1]) } ; 

const double EXP_DT_TAU_MEM[2] = { exp(-DT/TAU_MEM[0]) , exp(-DT/TAU_MEM[1]) } ; 

const double TAU_SYN[4] = {3.0, 2.0, 3.0, 2.0} ; 
const double EXP_DT_TAU_SYN[4] = { exp(-DT/TAU_SYN[0]) , exp(-DT/TAU_SYN[1]), exp(-DT/TAU_SYN[2]), exp(-DT/TAU_SYN[3])} ; 

#define IF_NMDA 1
const double TAU_NMDA[4] = {40.0, 20.0, 40.0, 20.0} ; 
const double EXP_DT_TAU_NMDA[4] = { exp(-DT/TAU_NMDA[0]) , exp(-DT/TAU_NMDA[1]), exp(-DT/TAU_NMDA[2]), exp(-DT/TAU_NMDA[3])} ; 
const double R_NMDA[2] = {1.0, 9.0} ; 

double *ISI ; 

string str_volt ; 
ofstream file_volt ; 

string str_spike_times ; 
ofstream file_spike_times ; 

////////////////////////////////////// 
// STP globals 
////////////////////////////////////// 
#define IF_STP 1 

#define IF_MARKRAM 0
#define IF_MATO 1 
#define IF_MONGILLO 0 

#define TAU_FAC (double) 450 
#define TAU_REC (double) 200 
#define USE (double) 0.03 

const double EXP_DT_TAU_FAC = exp(-DT/TAU_FAC) ; 
const double DT_OVER_TAU_REC = DT/TAU_REC ; 

//////////////////////////////////////
// other globals 
//////////////////////////////////////

#define IF_HYSTERESIS 0 
#define HYST_J_EE 0 
#define HYST_M0 1 
#define HYST_DX (double) 0.001 
#define HYST_X_MIN (double) 0.00 
#define HYST_X_MAX (double) 0.02 

#endif 
