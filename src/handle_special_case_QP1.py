import numpy as np
import matplotlib.pyplot as plt
from src.get_ionogram import get_ionogram
from src.utils import compute_ab, generate_values_exp, generate_fp_profile, L2_ERROR

params = {'mew': 0, 'mec': 'none', 'alpha': 0.8, 'markersize':1}

def find_QP1(original_f, original_vh, QP, num_points):
    """
    Finds suitable parameters for the E layer/0 QP layer
    """
    rm = QP['r_m0']
    possible_ym = np.linspace(1.0,350.0,num=400)
    process = {'error': 1e9, 'fc': -1, 'a':-1, 'rm': -1, 'b': -1, 'numt':-1, 'fp': None, 'rh': None}
    
    fp, rh = QP['plasma_frequency'], QP['real_height']

    epsilon = 1e-10
    numt = QP['numt']-1
    tmp_probe = original_f[numt:numt+num_points]

    for ym in possible_ym:
        if ym >= rm: continue
        
        # generate the real heights for the curated data 
        gen_rh = generate_fp_profile(rm, ym, QP['plasma_frequency'][-1], fp_data_test=tmp_probe, qp_type='neg_to_pos')
        
        # append to the E layer 
        tmp_fp = np.append(fp,tmp_probe)
        tmp_rh = np.append(rh,gen_rh)

        _, hv = get_ionogram(tmp_rh, tmp_fp, tmp_probe-epsilon)

        error_i = L2_ERROR(original_vh[numt:numt+num_points],hv)
        if error_i<process['error']:
            a,b = compute_ab(fc=QP['plasma_frequency'][-1], rb=rm-ym, ym=ym)
            process['error'] = error_i
            process['fc'] = QP['plasma_frequency'][-1]
            process['a'] = a
            process['rm'] = rm
            process['b'] = b
            process['numt'] = numt+num_points
            process['fp'], process['rh'] = (tmp_fp, tmp_rh)
    
    QP['f_c1'] = process['fc']
    QP['a_1']  = process['a']
    QP['r_m1'] = process['rm']
    QP['b_1']  = process['b']
    QP['numt'] = process['numt']
    QP['plasma_frequency'], QP['real_height'] = (process['fp'], process['rh'])
    print('QP NUMBER = 1 | SPECIAL FUNCTION GETS ERROR = ', process['error'])
    return QP, process['error']