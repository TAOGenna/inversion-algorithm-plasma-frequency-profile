import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange, tqdm
from src.plot_results import plot_results
from src.get_ionogram import get_ionogram
from src.utils import generate_values_exp, find_rm, find_b, L2_ERROR, generate_fp_profile
from src.handle_special_case_QP1 import find_QP1


def handle_f_layer(QP, frq, vh):
    """
    Computes the F layer based on the ionogram data.
    """
    numt_max = 3
    print(f'AMOUNT OF POINTS PER QP = {numt_max}')

    # Handle special case of QP1 
    tmp1_QP, error1 = find_QP1(original_f=frq, original_vh=vh, QP=QP.copy(), num_points=numt_max)
    tmp2_QP, error2 = get_qp(qp_number = 1, QP=QP.copy(), data_f=frq, data_r=vh, numt_max=numt_max)
    QP = tmp1_QP if error1<error2 else tmp2_QP
    plot_results(frq, vh, QP['plasma_frequency'], QP['real_height'], filename='QP1 layer')

    # Handle the rest of cases 
    for tmr in trange(2,30):
        QP, _ = get_qp(qp_number=tmr, QP=QP.copy(), data_f=frq, data_r=vh, numt_max=numt_max)
    
    return QP

# 
def get_qp(qp_number, QP, data_f, data_r, numt_max):

    # Check if we already have enough points
    numt = QP['numt']
    if numt >= len(data_f): return QP

    # Set a lowerbound and upperbound depending on if it is a QP or anti-QP layer
    search_length = 100
    fc_lowerbound, fc_upperbound, qp_type = None, None, None
    epsilon = 1e-10

    if qp_number % 2 == 1:
        fc_lowerbound = max(epsilon,QP['f_c'+str(qp_number-1)] - 2.0)
        fc_upperbound = QP['f_c'+str(qp_number-1)] - epsilon
        qp_type = 'neg_to_pos'
    else:
        fc_lowerbound = QP['plasma_frequency'][-1] + epsilon
        fc_upperbound = QP['plasma_frequency'][-1] + 2.0
        qp_type = 'pos_to_neg'

    possible_fc = np.linspace(fc_lowerbound, fc_upperbound, num=search_length)

    # Variable to store the best set of parameters based on error L2
    store = {'error': 1e12, 'a_1': -1, 'b_0': -1, 'r_m1': -1, 'data_f': np.array([]), 'data_r': np.array([])}

    for fi in possible_fc:
        # Compute current QP parameters
        fi = fi * 1e6  # Hz
        ai = fi ** 2   # Hz^2
        prev_real = QP['real_height'][-1] * 1e3  # m | r_i: intersection between QP_{i-1} and QP_{i}

        # ai1: a_{i-1} | bi1: b_{i-1} | rmi1: r_{m,i-1}
        rmi = find_rm(ri=prev_real, a=ai, ai1=QP['a_'+str(qp_number-1)] * 1e12, bi1=QP['b_'+str(qp_number-1)] * 1e12, rmi1=QP['r_m'+str(qp_number-1)] * 1e3, tipo=qp_type)  # m
        bi  =  find_b(ri=prev_real, a=ai, ai1=QP['a_'+str(qp_number-1)] * 1e12, bi1=QP['b_'+str(qp_number-1)] * 1e12, rmi1=QP['r_m'+str(qp_number-1)] * 1e3, tipo=qp_type)  # Hz^2

        #print(fi/1e6, ai/1e12, bi/1e12, rmi/1e3)
        #print('---------------------------------------------------------')

        # Create a local plasma frequency profile to toy with it
        tmp_f, tmp_r = QP['plasma_frequency'] * 1e6, QP['real_height'] * 1e3

        # Select the frequencies from the data set whose real heights we want
        test_frq = np.array([]) 
        for idx_numt in range(numt_max):
            if numt + idx_numt >= len(data_f): break
            test_frq = np.append(test_frq, data_f[numt + idx_numt] * 1e6)

        # Use logarithmically separated points to construct the fp profile 
        
        #curated_data_f = generate_values_exp(A=test_frq[0]/1e6, B=test_frq[-1]/1e6-epsilon, N=10, k=3) * 1e6
        curated_data_f = test_frq
        generated_rh   = generate_fp_profile(rm=rmi, ym=-1, fo = fi, fp_data_test=curated_data_f, qp_type=qp_type, b=bi)
        
        tmp_f, tmp_r = np.append(tmp_f,curated_data_f), np.append(tmp_r, generated_rh)

        # Compute a new ionogram with the syntethic fp profile
        fp_ans, hv_ans = get_ionogram(h=tmp_r / 1e3, fp=tmp_f / 1e6, range_f=test_frq/1e6-epsilon)

        # Compute error L2 | compare with original virtual heights from the ionogram's data
        L2_error = L2_ERROR(cmp1 = data_r[numt:numt+len(hv_ans)], cmp2 = hv_ans)

        # If the error with the current f_c is lesser than the previous one, update the store
        if L2_error < store['error']:
            i = str(qp_number)
            store['error'], store['a_'+i], store['b_'+i], store['r_m'+i] = L2_error, ai/1e12, bi/1e12, rmi/1e3
            store['data_r'], store['data_f'] = generated_rh/1e3, curated_data_f/1e6
    
    # Print some info
    final_error = store['error']
    tqdm.write(f'QP NUMBER = {qp_number} | Error = {final_error:.4f}')
    
    # Update QP dictionary with the best parameters found
    i = str(qp_number)
    QP['a_' +i]  = store['a_'+i]
    QP['f_c'+i]  = np.sqrt(store['a_'+i])
    QP['b_' +i]  = store['b_'+i]
    QP['r_m'+i]  = store['r_m'+i]
    QP['numt']   = QP['numt'] + numt_max
    QP['plasma_frequency'] = np.append(QP['plasma_frequency'], store['data_f'])
    QP['real_height'] = np.append(QP['real_height'], store['data_r'])
    return QP, store['error']

