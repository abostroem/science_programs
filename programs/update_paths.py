import os
import glob
import numpy as np
import subprocess as sp

if __name__ == "__main__":
    paths = np.genfromtxt('path_config.dat', dtype = 'str', delimiter = ',')
    for old_path, new_path in paths:
        sed_str = "sed -i 's:%s:%s:' *.py" %(old_path, new_path)
        sp.call(sed_str, shell = True)
