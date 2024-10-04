import numpy as np
import matplotlib.pyplot as plt

def generate_values(A, B, N, p):
    """
    Generates N values in the range [A, B] with higher density near B.

    Parameters:
    A (float): Lower bound of the range.
    B (float): Upper bound of the range.
    N (int): Number of values to generate.
    p (float): Power parameter (p > 1) to control the concentration.

    Returns:
    List[float]: List of generated values.
    """
    return [A + (B - A) * ((i / N) ** (1 / p)) for i in range(1, N + 1)]

import math

def generate_values_exp(A, B, N, k):
    """
    Generates N values in the range [A, B] with higher density near B using exponential mapping.

    Parameters:
    A (float): Lower bound of the range.
    B (float): Upper bound of the range.
    N (int): Number of values to generate.
    k (float): Exponential rate parameter (k > 0).

    Returns:
    List[float]: List of generated values.
    """
    return [A + (B - A) * (1 - math.exp(-k * i / N)) / (1 - math.exp(-k)) for i in range(1, N + 1)]

arr = generate_values(A=0.1, B=2.4, N=30, p=3)
arr = np.array(arr)

arr2 = generate_values_exp(A=0.1, B=2.4, N=30, k=1)
arr2 = np.array(arr2)

plt.plot(arr, [3]*len(arr),'o',markersize=2)
plt.plot(arr2, [2]*len(arr),'o',markersize=2)
plt.xlim((0,2.4))
plt.show()