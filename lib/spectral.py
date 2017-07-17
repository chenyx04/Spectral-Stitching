#!/usr/bin/python
from numpy.linalg import eig
import numpy as np
"""
Functions for running pure spectral algorithm on small blocks
"""
def spectral(We, Q):
    eigenv, V = eig(We)
    n = We.shape[0]
    I = sorted(range(len(eigenv)), key=lambda k: eigenv[k], reverse=True)
    V1 = V[:,I]
    evec = V1[:,0:(Q-1)].flatten()
    idxt = np.ones((n,1))
    idxt[evec<0] = -1
    idxt[evec==0] = 0
    return idxt