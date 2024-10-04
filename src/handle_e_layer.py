import numpy as np
import matplotlib.pyplot as plt
from src.get_ionogram import get_ionogram
from src.utils import compute_ab, generate_values_exp, generate_fp_profile, L2_ERROR
from tqdm import tqdm

def find_Q0(foE, probe_f, original_vh):
    """
    Finds suitable parameters for the E layer/0 QP layer
    """
    possible_rmE = np.linspace(1.0,350.0,num=400)
    possible_ymE = np.linspace(1.0,350.0,num=400)
    process = {'error': 1e9, 'rm': -1, 'ym': -1, 'fp': None, 'rh': None}
    
    epsilon = 1e-13
    num_points = 50
    curated_data_f = generate_values_exp(epsilon,foE-epsilon, num_points, k = 3) 
     
    for rmE in tqdm(possible_rmE):
        for ymE in possible_ymE:
            if ymE >= rmE: continue
            
            rh = generate_fp_profile(rmE, ymE, foE, fp_data_test=curated_data_f, qp_type='pos_to_neg')
            probe_f, hv = get_ionogram(rh, curated_data_f, probe_f-epsilon)

            error_i = L2_ERROR(original_vh,hv)
            if error_i<process['error']:
                process['error'] = error_i
                process['rm'], process['ym'] = (rmE, ymE)
                process['fp'], process['rh'] = (curated_data_f, rh)

    print('E LAYER ERROR = ',process['error'])
    return process['rm'], process['ym'], process['fp'], process['rh']

def handle_e_layer(frq,vh,foE):
    
    data_f   = frq[frq <= foE]
    data_hv  = vh[frq <= foE]
    index = len(data_f)
    epsilon = 1e-13
    
    rmE, ymE, fp_add, rh_add = find_Q0(foE=foE, probe_f=data_f, original_vh=data_hv)
    fp_add += epsilon # for numerical stability when computing ionogram
    a,b = compute_ab(fc=foE, rb=rmE-ymE, ym=ymE)
    QP = {'plasma_frequency': fp_add, 'real_height': rh_add, 
            'numt': index, 'a_0': a, 'r_m0': rmE, 'b_0': b,
            'f_c0': np.sqrt(a)}
    return QP
