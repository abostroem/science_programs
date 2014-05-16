'''
This program was written to match IDLs C_CORRELATE function

The C_CORRELATE function computes the cross correlation Pxy(L) or cross covariance Rxy(L) of two sample populations X and Y as a function of the lag L
'''
import numpy as np
import pdb
def c_corr(x, y, lag = False):
    #check for np array
    if not lag:
        lag = np.arange(np.shape(x)[0] + np.shape(y)[0] - 1) - np.shape(x)[0] + 1
    p = np.float_(np.ones(len(lag)))
    n = np.shape(x)[0]
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    for i, l in enumerate(lag):
        #pdb.set_trace()
        denom1 = np.sum((x - mean_x)**2)
        denom2 = np.sum((y - mean_y)**2)
        if l < 0:
            numerator = np.sum((x[abs(l):n] - mean_x)*(y[:n-abs(l)] - mean_y))
        if l >= 0:
            numerator = np.sum((x[:n-l] - mean_x)*(y[l:] - mean_y))
        #pdb.set_trace()
        p[i] = numerator/np.sqrt(denom1*denom2)
    return lag, p
