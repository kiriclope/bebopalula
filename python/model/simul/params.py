import warnings 
warnings.filterwarnings("ignore") 
import os 
import numpy as np 
import seaborn as sns

from plot_settings import *
SetPlotParams()

global n_pop, n_neurons, K, m0, folder, model 
# model = 'binary' 
model = 'lif' 
n_pop = 2
n_neurons = 4
K = 2000
# folder='shrisha'
folder = 'christos_on_2'
# folder = 'christos_Ie_840.00_Jee_800.00' # bifurcation for N4K2 n_frac 0.8 no STP
# folder = 'christos_on' # bifurcation for N4K2 n_frac 0.8 no STP 
m0 = .001 
nu0 = .05 

global IF_LOOP_J0 
IF_LOOP_J0 = 0 

global E_frac, n_frac, n_size 
if n_pop==2:
    E_frac = 0.8
    n_frac = [E_frac, np.round(1.0-E_frac,2)]
    n_size = [int(n_frac[0] * n_neurons * 10000.0) , int(n_frac[1] * n_neurons * 10000.0)] 
else:
    E_frac = 1.0
    n_frac = [E_frac, np.round(1.0-E_frac,2)]
    n_size = [int(n_frac[0] * n_neurons * 10000.0)]
    
global IF_STRUCTURE, IF_SPEC, IF_RING, IF_LOW_RANK 
IF_SPEC = 0 
IF_RING = 1
IF_LOW_RANK = 0
IF_GAUSS = 0 

IF_STRUCTURE = IF_SPEC or IF_RING or IF_LOW_RANK or IF_GAUSS

global SIGMA
SIGMA = [60, 60, 70, 60] 

global RANK, MAP 
RANK = 1
MAP = 0 

global KAPPA, KAPPA_1 
KAPPA = .25 # 4 # 12 #14 
KAPPA_1 = 4 # 8 # 12 

KAPPA_VAR = 2
FIX_SEED = 1 
MAP_SEED = 4 

global MEAN_KSI, VAR_KSI 
MEAN_KSI = 0 
VAR_KSI = 1 
SEED_KSI = 3

global MEAN_KSI_1, VAR_KSI_1 
MEAN_KSI_1 = 0 
VAR_KSI_1 = 1 

global IF_STP 
IF_STP = 1 
TAU_FAC = 450 
TAU_REC = 200 
USE = .03 

global IF_TRIALS, TRIAL_ID 
IF_TRIALS = 0 
TRIAL_ID = 1
N_TRIALS = 10

global IF_INI_COND, INI_COND_ID 
IF_INI_COND = 0
INI_COND_ID = 1

global IF_HYSTERESIS, HYST_JEE, HYST_M0 
IF_HYSTERESIS = 0 
HYST_JEE = -1 
HYST_M0 = 0 

global IF_DPA, IF_DUAL, IF_ODR, IF_DRT 
IF_DPA = 0
IF_DUAL = 0
IF_DRT = 0
IF_ODR = 0 

global KAPPA_EXT, PHI_EXT 
KAPPA_EXT = 1.0 
PHI_EXT = 0.25 

global T_SAMPLE_ON, T_SAMPLE_OFF
T_SAMPLE_ON = 2
T_SAMPLE_OFF = 3

global T_DIST_ON, T_DIST_OFF
T_DIST_ON = 4.5
T_DIST_OFF = 5.5

global T_TEST_ON, T_TEST_OFF
T_TEST_ON = 9
T_TEST_OFF = 10

global T_CUE_ON, T_CUE_OFF
T_CUE_ON = 6.5
T_CUE_OFF = 7.5

global MEAN_XI, VAR_XI
MEAN_XI = 0.0 
VAR_XI = 1.0 

global ext_inputs, J, J2, Tsyn, JEE, JEE2
ext_inputs = [] 
J = []
J2 = []
Tsyn = [] 
JEE = 0
JEE2 = 0 

global pal
try:
    folder.index('on')
    pal = ['r', 'b']
except:
    pal = [sns.color_palette('colorblind')[2],
        sns.color_palette('colorblind')[0],
        sns.color_palette('colorblind')[1],
        sns.color_palette('colorblind')[3]
        ] 

global label
label = ['E', 'I'] 
 
global J0, I0 
J0 = -1.0
I0 = 1.0 

if n_pop==1 :
    folder = 'I0_%.2f_J0_%.2f' % (I0, J0) 

global filename 
filename = 'inputs.dat'

global IF_LEFT_RIGHT, MEAN_XI_LEFT, MEAN_XI_RIGHT, VAR_XI_LEFT, VAR_XI_RIGHT, RHO
IF_LEFT_RIGHT = 0
MEAN_XI_LEFT = -0.0
VAR_XI_LEFT = 5.0 
MEAN_XI_RIGHT = -0.0
VAR_XI_RIGHT = 5.0 
RHO = 0.5

global IF_FF, MEAN_FF, VAR_FF, VAR_ORTHO, IF_RHO_FF, RH0_FF_XI
IF_FF = 0 
MEAN_FF = 1.0
VAR_FF = 1.0
VAR_ORTHO = 0.0

IF_RHO_FF = 1
RHO_FF_XI = 1.0 

global path, con_path, fig_path, ksi_path
path=''
con_path=''
ksi_path=''
fig_path = ''

def init_param():
    global folder, n_pop, n_neurons, K, ext_inputs, J, J2, m0, nu0, Tsyn, path, con_path, KAPPA, KAPPA_1 , fig_path, MAP 
    global TRIAL_ID, INI_COND_ID , HYST_JEE, IF_HYSTERESIS, HYST_M0, IF_INI_COND , JEE, JEE2, ksi_path, con_path
    
    if(n_pop!=1): 
        file_name = "/homecentral/alexandre.mahrach/IDIBAPS/cpp/model/parameters/%dpop/%s.txt" % (n_pop, folder) 
        print("reading parameters from:", file_name) 
        
        i=0 
        with open(file_name, 'r') as file:  # Open file for read 
            for line in file: # Read line-by-line
                line = line.strip().split()  # Strip the leading/trailing whitespaces and newline
                line.pop(0)
                if i==0:
                    ext_inputs = np.asarray([float(j) for j in line])
                    ext_inputs *= nu0
                    # print(ext_inputs)
                if i==1:
                    J = np.asarray([float(j) for j in line])
                    J = J.reshape(n_pop, n_pop)
                    JEE = J[0][0] 
                    J2 = J*J 
                    JEE2 = J2[0][0] 
                    # print(J)
                if i==2:
                    Tsyn = np.asarray([float(j) for j in line])
                    Tsyn = Tsyn.reshape(n_pop, n_pop)
                    # print(Tsyn)
                i=i+1
    else:
        ext_inputs = I0
        J = J0 
        J2 = J0*J0
        Tsyn = 2 ;
    
    path = '/homecentral/alexandre.mahrach/IDIBAPS/cpp/model/simulations/%s/%dpop/%s' % (model, n_pop, folder) 
    fig_path = '/homecentral/alexandre.mahrach/IDIBAPS/python/model/simul/figures/%s/%dpop/%s' % (model, n_pop, folder) 
    con_path = '/homecentral/alexandre.mahrach/IDIBAPS/connectivity/%dpop' % n_pop 
    
    if(K!=np.inf):
        
        if(n_pop==1):
            path+= '/N%d/K%d' % (n_size[0]/1000, K) 
            con_path+= '/N%d/K%d' % (n_size[0]/1000, K) 
            fig_path+= '/N%d/K%d' % (n_size[0]/1000, K) 
        else:
            path+= '/NE_%d_NI_%d/K%d' % (n_size[0]/1000, n_size[1]/1000, K) 
            fig_path+= '/NE_%d_NI_%d/K%d' % (n_size[0]/1000, n_size[1]/1000, K) 
            con_path+= '/NE_%d_NI_%d/K%d' % (n_size[0]/1000, n_size[1]/1000, K)
            
    else:
        path+= '/N_inf/K_inf' 
        con_path+= '/N_inf/K_inf'  
        fig_path+= '/N_inf/K_inf'  
        
    if(IF_STP):
        path += '/STP/Tf_%d_Tr_%d_U_%.2f' % (TAU_FAC, TAU_REC, USE) 
    
    if(IF_STRUCTURE):
        if(IF_SPEC):
            if(RANK==1):
                path += '/spec/kappa_%.2f' % KAPPA 
                fig_path += '/spec/kappa_%.2f' % KAPPA 
                con_path += '/spec/kappa_%.2f' % KAPPA 
            elif(RANK==2):
                path += '/spec/kappa_%.2f_kappa_1_%.2f' % (KAPPA, KAPPA_1)
                fig_path += '/spec/kappa_%.2f_kappa_1_%.2f' % (KAPPA, KAPPA_1)
                con_path += '/spec/kappa_%.2f_kappa_1_%.2f' % (KAPPA, KAPPA_1)
                if(FIX_SEED):
                    path += '/seed_%d' % (MAP_SEED) 
                    fig_path += '/seed_%d' % (MAP_SEED) 
                    con_path += '/seed_%d' % (MAP_SEED) 
                
        if(IF_RING):
            if(RANK==1):
                path += '/ring/kappa_%.2f' % KAPPA                
                con_path += '/ring/kappa_%.2f' % KAPPA                
                fig_path += '/ring/kappa_%.2f' % KAPPA                
            elif(RANK==2):
                path += '/ring/kappa_%.2f_kappa_1_%.2f' % (KAPPA, KAPPA_1)
                con_path += '/ring/kappa_%.2f_kappa_1_%.2f' % (KAPPA, KAPPA_1)
                fig_path += '/ring/kappa_%.2f_kappa_1_%.2f' % (KAPPA, KAPPA_1)

        if(IF_GAUSS):
            path += '/gauss/EE_%d_EI_%d_IE_%d_II_%d' % (SIGMA[0], SIGMA[1], SIGMA[2], SIGMA[3]) 
            con_path += '/gauss/EE_%d_EI_%d_IE_%d_II_%d' % (SIGMA[0], SIGMA[1], SIGMA[2], SIGMA[3]) 
            fig_path += '/gauss/EE_%d_EI_%d_IE_%d_II_%d' % (SIGMA[0], SIGMA[1], SIGMA[2], SIGMA[3]) 
            
        if(IF_LOW_RANK):
            ksi_path = con_path
            
            if(RANK==1):
                path += '/low_rank/kappa_%.2f' % KAPPA 
                fig_path += '/low_rank/kappa_%.2f' % KAPPA 
                con_path += '/low_rank/kappa_%.2f' % KAPPA
                ksi_path += '/low_rank/rank_1/seed_ksi_%d' % SEED_KSI
                
                # path += '/low_rank/kappa_%.2f/mean_%.2f_var_%.2f' % (KAPPA,MEAN_KSI, VAR_KSI) 
                # fig_path += '/low_rank/kappa_%.2f/mean_%.2f_var_%.2f' % (KAPPA,MEAN_KSI, VAR_KSI) 
                # con_path += '/low_rank/kappa_%.2f/mean_%.2f_var_%.2f' % (KAPPA,MEAN_KSI, VAR_KSI) 
            elif(RANK==2):
                path += '/low_rank/kappa_%.2f_kappa_1_%.2f' % (KAPPA, KAPPA_1) 
                fig_path += '/low_rank/kappa_%.2f_kappa_1_%.2f' % (KAPPA, KAPPA_1) 
                con_path += '/low_rank/kappa_%.2f_kappa_1_%.2f' % (KAPPA, KAPPA_1) 
                ksi_path += '/low_rank/rank_2/seed_ksi_%d' % SEED_KSI 
                # path += '/low_rank/kappa_%.2f/mean_%.2f_var_%.2f/mean_%.2f_var_%.2f' % (KAPPA, MEAN_KSI, VAR_KSI, MEAN_KSI, VAR_KSI) 
                # fig_path += '/low_rank/kappa_%.2f/mean_%.2f_var_%.2f/mean_%.2f_var_%.2f' % (KAPPA, MEAN_KSI, VAR_KSI, MEAN_KSI, VAR_KSI) 
                # con_path += '/low_rank/kappa_%.2f/mean_%.2f_var_%.2f/mean_%.2f_var_%.2f' % (KAPPA, MEAN_KSI, VAR_KSI, MEAN_KSI, VAR_KSI) 
        
    if(IF_DPA) : 
        path += '/DPA/kappa_%.2f_phi_%.2f' % (KAPPA_EXT, PHI_EXT)  
        fig_path += '/DPA/kappa_%.2f_phi_%.2f' % (KAPPA_EXT, PHI_EXT) 
    elif(IF_DUAL) :
        path += '/dual_task/kappa_%.2f_phi_%.2f' % (KAPPA_EXT, PHI_EXT)
        fig_path += '/dual_task/kappa_%.2f_phi_%.2f' % (KAPPA_EXT, PHI_EXT)
    elif(IF_DRT) :
        path += '/DRT/kappa_%.2f_phi_%.2f' % (KAPPA_EXT, PHI_EXT) 
        fig_path += '/DRT/kappa_%.2f_phi_%.2f' % (KAPPA_EXT, PHI_EXT) 
    elif(IF_ODR) :
        path += '/ODR'
        fig_path += '/ODR'
    
    if(IF_HYSTERESIS):
        fig_path += '/hysteresis'; 
        if(HYST_JEE==1):
            path += '/hysteresis/Jee_up'; 
        if(HYST_JEE==-1):
            path += '/hysteresis/Jee_down';
            
        if(HYST_M0==1):
            path += '/hysteresis/m0_up'; 
        if(HYST_M0==-1):
            path += '/hysteresis/m0_down'; 

    if(IF_LOOP_J0):
        path += '/m0_%.3f' % m0 ;
                
    if(IF_TRIALS):
        path += '/trial_%d' % TRIAL_ID ;
        con_path += '/trial_%d' % TRIAL_ID ;
        fig_path += '/trial_%d' % TRIAL_ID ;
        
    if(IF_INI_COND):
        path += '/ini_cond_%d' % INI_COND_ID ; 
        fig_path += '/ini_cond_%d' % INI_COND_ID ;     
            
    if(~IF_TRIALS): 
        print('reading simulations data from:', path) 
    
    if not os.path.isdir(fig_path):
        os.makedirs(fig_path)

