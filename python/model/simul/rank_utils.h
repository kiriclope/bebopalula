import numpy as np

def get_idx_perm():
dtype = np.dtype("L") 
  
    try: 
        with open(gv.con_path + 'idx_perm.dat', 'rb') as f: 
            idx_perm = np.fromfile(f, dtype)            
        f.close()
        
    except EOFError:
        pass
    
print('idx_perm', idx_perm.shape, idx_perm[:10], idx_perm[-10:]) 
