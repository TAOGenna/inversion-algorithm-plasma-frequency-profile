import numpy as np
import os

def generate_values_exp(A, B, N, k):
    """
    Generates N values in the range [A, B] with higher density near B using exponential mapping.

    Parameters:
    A (float): Lower bound of the range.
    B (float): Upper bound of the range.
    N (int): Number of values to generate.
    k (float): Exponential rate parameter (k > 0).

    Returns:
    numpy.array[float]: List of generated values.
    """
    return np.array([A + (B - A) * (1 - np.exp(-k * i / N)) / (1 - np.exp(-k)) for i in range(1, N + 1)])

def generate_fp_profile(rm, ym, fo, fp_data_test, qp_type, b = None):
    a = None
    if b is None:
        a, b = compute_ab(fc=fo, rb=rm-ym, ym=ym)
    else:
        a = fo ** 2
    rh = np.array([get_real_height(frq=fp_i, a=a, b=b, rm=rm, qp_type=qp_type) for fp_i in fp_data_test])
    return rh

def solve_for_r(f, a, b, rm, qp_type):
    f = f**2
    discriminant = (f - a) / b if qp_type == 'neg_to_pos' else -(f - a) / b
    if discriminant < 0: return -1, -1
    s = np.sqrt(discriminant)
    return rm / (1 - s), rm / (1 + s)

def get_real_height(frq, a, b, rm, qp_type):
    r1, r2 = solve_for_r(frq, a, b, rm, qp_type)
    if r1 == -1 and r2 == -1: return -1
    return max(r1, r2) if qp_type == 'neg_to_pos' else min(r1, r2)

def compute_ab(fc, rb, ym):
    """
    need f_ci, r_bi, y_mi
    """
    a = fc**2
    b = a*((rb/ym)**2)
    return a,b

def y_der_negative(r, b, rm):
    return -2 * b * rm * (1 - rm / r) * (1 / r**2)

def y_der_positive(r, b, rm):
    return 2 * b * rm * (1 - rm / r) * (1 / r**2)

def y_negative(r, a, b, rm):
    return a - b * ((1 - rm / r)**2)

def y_positive(r, a, b, rm):
    return a + b * ((1 - rm / r)**2)

def find_rm(ri, a, ai1, bi1, rmi1, tipo):
    deriv = y_der_negative(ri, bi1, rmi1) if tipo == 'neg_to_pos' else y_der_positive(ri, bi1, rmi1)
    y_type = y_negative(ri, ai1, bi1, rmi1) if tipo == 'neg_to_pos' else y_positive(ri, ai1, bi1, rmi1)

    return (ri**2 * deriv) / (2 * (y_type - a) + ri * deriv)

def find_b(ri, a, ai1, bi1, rmi1, tipo):
    deriv = y_der_negative(ri, bi1, rmi1) if tipo == 'neg_to_pos' else y_der_positive(ri, bi1, rmi1)
    y_type = y_negative(ri, ai1, bi1, rmi1) if tipo == 'neg_to_pos' else y_positive(ri, ai1, bi1, rmi1)
    cte = 1.0 if tipo == 'neg_to_pos' else -1.0

    numerator = 2 * (y_type - a) + ri * deriv
    denominator = 4 * (y_type - a)
    return cte * (numerator**2) / denominator

def L2_ERROR(cmp1, cmp2):
    assert len(cmp1) == len(cmp2), 'the lenghts of compare1 and compare2 are different'
    difference = abs(cmp1 - cmp2)
    errors = np.array([difference])
    return np.sqrt(np.sum(errors**2)) / len(errors)

def check_E_layer_existance(frq,vh,foE):
    index_fE = np.searchsorted(frq, foE, side='right')
    return index_fE >= 2

def print_centered(text):
    # Get the terminal width
    terminal_width = os.get_terminal_size().columns
    
    # Calculate the center position
    centered_text = text.center(terminal_width)
    
    # Print the text in bright green using ANSI escape codes
    print(f'\033[1;92m{centered_text}\033[0m')  # \033[1;92m is for bright green text, \033[0m resets the color

def print_yellow(text):
    # Print the text in yellow using ANSI escape codes
    print(f'\033[93m{text}\033[0m')  # \033[93m is for yellow text, \033[0m resets the color
