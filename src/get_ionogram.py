import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import constants


# Physical constants
qe = constants.e
e0 = constants.epsilon_0
me = constants.electron_mass
pi = constants.pi
eps = 1e-20


def ne_to_fp2(ne): # input: m^{-3} output: Hz^{2}
    return (qe**2/(e0*me)/(2*pi)**2) * ne

def fp2_to_ne(fp2): # input: Hz^{2} output: m^{-3}
    return ((2*pi)**2)*(e0*me)/(qe**2) * fp2

def get_ionogram(h,fp,range_f):
    """
    INPUT:
    - h: real heights
    - fp: plasma frequency
    - range_f: frequency range in which we are interested in
    IMPORTANT: `fp` and `range_f` should have the same units

    OUTPUT:
    - f : frequency range | same values as in `range_f`
    - hv: virtual heights | same units as in `h`
    """

    ne = fp2_to_ne(fp**2)

    # some modifications for numerical stability
    ne[0]=0.00
    epsilon = 1e-13
    range_f = range_f - epsilon
    # end of modifications

    if h[0]!=0:
        h  = np.append(0,h)
        ne = np.append(0,ne)

    n_elements = len(h)
    dif_h = h[1:len(h)] - h[0:len(h)-1]

    # arrays where we will save the integral results per frequency
    f  = np.array([])
    hv = np.array([])

    fp2 = ne_to_fp2(ne)
    max_possible_fp = max(fp2)

    for f_i in range_f:
        f_probe = f_i**2
        if f_probe<max_possible_fp:
            n2l = 1 - fp2[0]/f_probe

            integral = 0
            if n2l>0:
                sqrt_n2l = np.sqrt(n2l)
                for i in range(n_elements-1):
                    n2r = 1 - fp2[i+1]/f_probe
                    if n2r>0:
                        sqrt_n2r = np.sqrt(n2r)
                        integral += 2*dif_h[i]/(sqrt_n2l+sqrt_n2r)
                        n2l = n2r
                        sqrt_n2l=sqrt_n2r
                    else:
                        integral += 2*dif_h[i]*sqrt_n2r/(n2l-n2r)
                        break
            f  = np.append(f,f_i)       # Hz
            hv = np.append(hv,integral) # m
    return f, hv # (Hz,m)