import numpy as np
import sys
import pdb

def a(x, lag):

#for i in lag:
    #    assert type(i) == int, print 'all elements of lag must be integer'
    p= np.empty((0,))
    x = np.array(x)
    lag = np.array(lag)
    nx = len(x)
    data = x - np.mean(x)
    lag = np.abs(lag)
    for k, l in enumerate(lag):
        #pdb.set_trace()
        psum_top =  np.sum(data[0:nx - l] * data[l:])
        psum_bottom =  np.sum(data**2)
        #pdb.set_trace()
        p = np.append(p, psum_top/psum_bottom)
        
    return p

if __name__ == "__main__":
    x = sys.argv[1]
    lag = sys.argv[2]
