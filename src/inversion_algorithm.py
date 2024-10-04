# INPUT
'''
author: Renzo Kenyi Takagui Perez 
date  : 10 - 02 - 2024
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
from scipy import constants
from scipy import integrate
from scipy.interpolate import PchipInterpolator
import numpy.polynomial.chebyshev as cheb
import inspect
import time
import tkinter as tk  # or use pyautogui or wx

class InversionAlgorithm:
    # Physical constants
    qe = constants.e
    e0 = constants.epsilon_0
    me = constants.electron_mass
    pi = constants.pi
    eps = 1e-6

    def __init__(self, data_r, data_f, foE):
        """
        Initialization of the InversionAlgorithm class.
        - data_r: should receive in [meters]
        - data_f: should receive in [Hz]
        - foE   : should receive in [MHz]
        """
        # Convert inputs to appropriate units
        self.data_r = np.array(data_r)       # Convert to Km  (except if they from .sao files)
        self.data_f = np.array(data_f)       # Convert to MHz (except if they from .sao files)
        self.foE = foE                       # MHz

        # Initialize QP dictionary for storing results
        self.QP = {}
        
        # Initialize variables storing the final reconstructed profile
        self.plasma_frequency = None 
        self.real_height      = None
        self.frequency        = None
        self.virtual_height   = None

        # Run the inversion algorithm on initialization
        self.frequency, self.virtual_height     = self.run_inversion()
        self.plasma_frequency, self.real_height = self.QP['plasma_frequency'], self.QP['real_height']


    def run_inversion(self):
        """
        Executes the inversion algorithm to generate the ionogram and plasma frequency profile.
        """
        self._compute_E_layer()
        self._compute_F_layer()

        # Find the produced final ionogram
        fp_ans, hv_ans = self.get_ionogram(self.QP['real_height'], self.QP['plasma_frequency'], self.QP['plasma_frequency'])
        
        return fp_ans, hv_ans
        # OPTIONAL: Plot final comparison
        # self.plot_ionogram(fp_ans, hv_ans)


    def E_ionogram(self, f,r_bE,y_mE,f_c):
        return r_bE + 0.5*y_mE*(f/f_c)*np.log((f_c + self.eps + f)/(f_c + self.eps - f))

    def delta_he(self, fk,zE,yE,fE):
        """
        E layer height contribution for fk < fE
        """
        return  zE-yE+(1/2)*(yE)*(fk/fE)*np.log(((fE + self.eps)+fk)/((fE+self.eps)-fk))

    def fp_E_layer(self, rme, yme, fce):
        rbe = rme - yme
        r = np.linspace(rbe,rme,num=100)
        inside_sqrt = 1 - (((r-rme)/yme)**2) * ((rbe/r)**2)
        inside_sqrt[inside_sqrt<0] = 0
        fp = np.sqrt( (fce**2)*( inside_sqrt ) )
        return r,fp

    def compute_ab(self, fc, rb, ym):
        """
        need f_ci, r_bi, y_mi
        """
        a = fc**2
        b = a*((rb/ym)**2)
        return a,b
    
    def find_Q0(self, E_data_f, E_data_r, fE):
        data_zE = np.linspace(1.0,250.0,num=400)
        data_yE = np.linspace(1.0,200.0,num=400)
        process = {'error': 1e9, 'rm': -1, 'rb': -1, 'ym': -1}
        for ze in data_zE:
            for ye in data_yE:

                if ye > ze:
                    continue

                P = self.delta_he(E_data_f,zE=ze,yE=ye,fE=fE)
                dif = E_data_r - P
                error_i = (np.sum(dif**2))/len(dif)

                if error_i<process['error']:
                    process['error'] = error_i
                    process['rm'], process['rb'], process['ym'] = (ze, ze-ye, ye)

        return process['rm'], process['rb'], process['ym']

    def _compute_E_layer(self):
        """
        Computes the E layer based on the initial ionogram data.
        """
        index_fE = np.searchsorted(self.data_f, self.foE, side='right')

        self.QP['f_c0'] = self.data_f[index_fE]
        self.QP['f_c0'] = self.foE
        self.QP['QPs'] = 0

        # Data related to the E layer
        E_data_f = self.data_f[self.data_f <= self.foE]
        E_data_r = self.data_r[self.data_f <= self.foE]

        self.QP['numt'] = len(E_data_f)
        self.QP['r_m0'], self.QP['r_b0'], self.QP['y_m0'] = self.find_Q0(E_data_f, E_data_r, self.foE)
        self.QP['a_0'], self.QP['b_0'] = self.compute_ab(fc=self.QP['f_c0'], rb=self.QP['r_b0'], ym=self.QP['y_m0'])

        # Interpolation to get the E layer plasma frequency profile
        r_E, fp_E = self.fp_E_layer(rme=self.QP['r_m0'], yme=self.QP['y_m0'], fce=self.foE)
        pchip = PchipInterpolator(fp_E, r_E)
        xnew = E_data_f
        ynew = pchip(xnew)

        self.QP['plasma_frequency'] = xnew
        self.QP['real_height'] = ynew

    def ne_to_fp2(self, ne): # input: m^{-3} output: Hz^{2}
        return (self.qe**2/(self.e0*self.me)/(2*self.pi)**2) * ne

    def fp2_to_ne(self, fp2): # input: Hz^{2} output: m^{-3}
        return ((2*self.pi)**2)*(self.e0*self.me)/(self.qe**2) * fp2

    def _compute_F_layer(self):
        """
        Computes the F layer based on the ionogram data.
        """
        numt_min = 2
        numt_max = 6
        start_time = time.time()
        for tmr in range(30):
            self.get_qp(qp_number=tmr + 1, data_f=self.data_f, data_r=self.data_r, numt_min=numt_min, numt_max=numt_max)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time for computing QP layers: {elapsed_time} seconds \n")

    def get_ionogram(self,h,fp,range_f):
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

        ne = self.fp2_to_ne(fp**2)

        # some modifications for numerical stability
        ne[0]=0.00
        epsilon = 1e-7
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

        fp2 = self.ne_to_fp2(ne)
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

    def get_qp(self, qp_number, data_f, data_r, numt_min, numt_max):
        # Ensure that 'plasma_frequency' is initialized
        if self.QP.get('plasma_frequency') is None:
            self.QP['plasma_frequency'] = np.array([])

        # Check if we already have enough points
        numt = len(self.QP['plasma_frequency'])
        if numt >= len(data_f): return

        # Set a lowerbound and upperbound depending on if it is a QP or anti-QP layer
        tam_without_grid = 1000
        fc_lowerbound, fc_upperbound, qp_type = None, None, None

        if qp_number % 2 == 1:
            fc_lowerbound = self.QP['f_c'+str(qp_number-1)] - 2.0
            fc_upperbound = self.QP['f_c'+str(qp_number-1)] - self.eps
            qp_type = 'neg_to_pos'
        else:
            fc_lowerbound = self.QP['plasma_frequency'][-1] + self.eps
            fc_upperbound = self.QP['plasma_frequency'][-1] + 2.0
            qp_type = 'pos_to_neg'

        possible_fc = np.linspace(fc_lowerbound, fc_upperbound, num=tam_without_grid)

        # Variable to store the best set of parameters based on error L2
        store = {'error': 1e12, 'a_1': -1, 'b_0': -1, 'r_m1': -1, 'data_f': np.array([]), 'data_r': np.array([])}

        for fi in possible_fc:
            # Compute current QP parameters
            fi = fi * 1e6  # Hz
            ai = fi ** 2   # Hz^2
            prev_real = self.QP['real_height'][-1] * 1e3  # m

            rmi = self.find_rm(ri=prev_real, a=ai, ai1=self.QP['a_'+str(qp_number-1)] * 1e12, bi1=self.QP['b_'+str(qp_number-1)] * 1e12, rmi1=self.QP['r_m'+str(qp_number-1)] * 1e3, tipo=qp_type)  # m
            bi = self.find_b(ri=prev_real, a=ai, ai1=self.QP['a_'+str(qp_number-1)] * 1e12, bi1=self.QP['b_'+str(qp_number-1)] * 1e12, rmi1=self.QP['r_m'+str(qp_number-1)] * 1e3, tipo=qp_type)  # Hz^2

            # Create a local plasma frequency profile to toy with it
            tmp_f, tmp_r = self.QP['plasma_frequency'] * 1e6, self.QP['real_height'] * 1e3

            # Data to be added to the final reconstructed real height profile IF its error is the smallest
            i_data_r = np.array([])
            i_data_f = np.array([])

            numt_for = numt_min if qp_number + 1 <= 10 else numt_max

            for idx_numt in range(numt_for):
                if numt + idx_numt >= len(data_f): break

                frq = data_f[numt + idx_numt] * 1e6  # Hz
                rh_frq = self.get_real_height(frq, ai, bi, rmi, prev_real, qp_type=qp_type)  # m
                prev_real = rh_frq

                tmp_f, tmp_r = np.append(tmp_f, frq), np.append(tmp_r, rh_frq)
                i_data_f, i_data_r = np.append(i_data_f, frq), np.append(i_data_r, rh_frq)

            fp_ans, hv_ans = self.get_ionogram(h=tmp_r / 1e3, fp=tmp_f / 1e6, range_f=data_f[:len(tmp_f)])

            # Compute initial error L2 | compare with original virtual heights from the ionogram's data
            original_vh = data_r[:len(fp_ans)]
            difference = abs(hv_ans - original_vh)
            errors = np.array([difference])
            L2_error = self.L2_ERROR(errors)

            # If the error with the current f_c1 is lesser than the previous one, update the store
            if L2_error < store['error']:
                store['error'] = L2_error
                store['a_'+str(qp_number)] = ai / 1e12
                store['b_'+str(qp_number)] = bi / 1e12
                store['r_m'+str(qp_number)] = rmi / 1e3
                store['data_r'] = i_data_r / 1e3
                store['data_f'] = i_data_f / 1e6
      
        # Print some info
        print('--------------------------------------------------------------------------------------------------------------')
        print(f'QP NUMBER = {qp_number}')
        print('Error = ',store['error'])
      
        # Update QP dictionary with the best parameters found
        self.QP['a_'+str(qp_number)] = store['a_'+str(qp_number)]
        self.QP['f_c'+str(qp_number)] = np.sqrt(store['a_'+str(qp_number)])
        self.QP['b_'+str(qp_number)] = store['b_'+str(qp_number)]
        self.QP['r_m'+str(qp_number)] = store['r_m'+str(qp_number)]
        self.QP['plasma_frequency'] = np.append(self.QP['plasma_frequency'], store['data_f'])
        self.QP['real_height'] = np.append(self.QP['real_height'], store['data_r'])


    
    def solve_for_r(self, f, a, b, rm, qp_type):
        f = f**2
        discriminant = (f - a) / b if qp_type == 'neg_to_pos' else -(f - a) / b
        if discriminant < 0: return -1, -1

        s = math.sqrt(discriminant)
        return rm / (1 - s), rm / (1 + s)

    def get_real_height(self, frq, a, b, rm, h_ref, qp_type):
        r1, r2 = self.solve_for_r(frq, a, b, rm, qp_type)
        if r1 == -1 and r2 == -1: return -1
        return max(r1, r2) if qp_type == 'neg_to_pos' else min(r1, r2)

    def y_der_negative(self, r, b, rm):
        return -2 * b * rm * (1 - rm / r) * (1 / r**2)

    def y_der_positive(self, r, b, rm):
        return 2 * b * rm * (1 - rm / r) * (1 / r**2)

    def y_negative(self, r, a, b, rm):
        return a - b * ((1 - rm / r)**2)

    def y_positive(self, r, a, b, rm):
        return a + b * ((1 - rm / r)**2)

    def find_rm(self, ri, a, ai1, bi1, rmi1, tipo):
        deriv = self.y_der_negative(ri, bi1, rmi1) if tipo == 'neg_to_pos' else self.y_der_positive(ri, bi1, rmi1)
        y_type = self.y_negative(ri, ai1, bi1, rmi1) if tipo == 'neg_to_pos' else self.y_positive(ri, ai1, bi1, rmi1)

        return (ri**2 * deriv) / (2 * (y_type - a) + ri * deriv)

    def find_b(self, ri, a, ai1, bi1, rmi1, tipo):
        deriv = self.y_der_negative(ri, bi1, rmi1) if tipo == 'neg_to_pos' else self.y_der_positive(ri, bi1, rmi1)
        y_type = self.y_negative(ri, ai1, bi1, rmi1) if tipo == 'neg_to_pos' else self.y_positive(ri, ai1, bi1, rmi1)
        cte = 1.0 if tipo == 'neg_to_pos' else -1.0

        numerator = 2 * (y_type - a) + ri * deriv
        denominator = 4 * (y_type - a)
        return cte * (numerator**2) / denominator

    def L2_ERROR(self, errors):
        return np.sqrt(np.sum(errors**2)) / len(errors)

    def plot_ionogram(self, fp_ans, hv_ans):
        """
        Plots the comparison between the original ionogram and the reconstructed ionogram.
        """
        plt.figure()
        plt.title(f'Final comparison')
        plt.xlabel('f [Mhz]')
        plt.ylabel('r [Km]')
        plt.plot(self.QP['plasma_frequency'], self.QP['real_height'], label='reconstructed fp profile', alpha=0.3, mew=0, mec='none')
        plt.plot(fp_ans, hv_ans, label='reconstructed ionogram', mew=0, mec='none')
        plt.plot(self.data_f, self.data_r, 'o', label='original ionogram', alpha=0.3, mew=0, mec='none')
        plt.legend()
        plt.show()



def check_E_layer_existance(frq,vh,foE):
    index_fE = np.searchsorted(frq, foE, side='right')
    return index_fE >= 2


def plot_inversion_algorithm(ionogram_f_original, ionogram_hv_original, reco_fp, reco_hr, reco_f, reco_hv, path, date, foE, i):
    """
    Plots the comparison between the original ionogram and the reconstructed ionogram.
    """
    plt.figure()
    plt.title(date + f' | foE = {foE}')
    plt.ylabel('Height [km]')
    plt.xlabel('Freq [$MHz$]')
    plt.plot(reco_fp, reco_hr, label='reconstructed fp profile', alpha=0.3, mew=0, mec='none')
    plt.plot(reco_f, reco_hv, label='reconstructed ionogram', mew=0, mec='none')
    plt.plot(ionogram_f_original, ionogram_hv_original, 'o', label='original ionogram', alpha=0.3, mew=0, mec='none')
    plt.legend()
    npath = path + '/' + str(i) + '.png'
    plt.savefig(npath)
    plt.close()