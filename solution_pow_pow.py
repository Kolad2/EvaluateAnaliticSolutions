import mpmath as mp
from qmath import QBinom, QPochhammer
import numpy as np


def c_pow_pow(n, a, b, t=1):
    res = 2 ** (-n*b) * (1 - 2 ** (-n*a)) * np.exp((1-2**(-n*a))*t)
    return res