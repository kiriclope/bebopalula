import numpy as np
from write import replace_global

for Ie in np.arange(2, 3., 0.1): 
    # for Jee in np.arange(1., 1.1, 0.1):
    Jee = 1 
    folder = "'L23_Ie_%.1f_Jee_%.1f'" % (Ie, Jee)
    replace_global('folder', folder)
    
    exec(open("./rates.py").read())
