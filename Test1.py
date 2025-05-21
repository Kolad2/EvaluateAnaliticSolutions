import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from qmath import QBinom, QPochhammer



def main():
    q = 0.5
    c = get_qc(10, 2, q)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(c)
    plt.show()

def get_c(n_max, t=1):
    c = np.zeros(n_max)
    return c

def get_qc(n_max, t = 1, q = 0.5):
    r"""
    c_{n}(t) = \frac{1}{(q;q)_n}\sum_{j = 0}^{n}(-1)^{j}
    \binom{n}{j}_q q ^ {j(j - 1) / 2}
    e^{-q^{n - j}t}
    """
    t =  mp.mpf(t)
    q = mp.mpf(q)
    qbinom = QBinom(q)
    qpoch = QPochhammer(q)
    c = np.zeros(n_max)
    for n in range(n_max):
        p = qpoch(q, n)
        term = 0
        sign = 1
        for j in range(min(n+1, 10)):
            term_1 = qbinom(n, j)
            term_2 = mp.exp(-mp.power(q, n - j)*t) # \exp(q^{n-j}*t)
            term_3 = mp.power(q, j * (j - 1) // 2) # q^{j(j-1)/2}
            term = term + sign*term_1 * term_2 * term_3
            sign = -sign
        c[n] = term / p
    return c



if __name__ == "__main__":
    main()