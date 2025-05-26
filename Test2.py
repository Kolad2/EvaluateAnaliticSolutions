import numpy as np
import matplotlib.pyplot as plt
from solutions import c_tr_pow_limited
from qmath import QBinom, QPochhammer
import mpmath as mp

def main():
    q = 0.5
    t = 0
    n_max = 40
    t = np.array([5., 25, 125, 625])
    lt = np.array([5, 10, 25])
    #
    n = np.arange(-n_max, 0, 0.1)*1.
    c = np.array([c_tr_pow_limited(_n, q) for _n in n])
    #
    qpoch = QPochhammer(q)
    qpoch_inf = qpoch(q, 41)  # (q;q)_\infty
    ln_qpoch_inf = mp.log(qpoch_inf)
    hat_n = mp.log(1 - ln_qpoch_inf) / mp.log(mp.mpf(q))

    #c_2 = np.exp(-q**n) / qpoch_inf
    #c_3 = 1 - q ** (-n)

    fig = plt.figure(figsize=(10, 5))
    ax = [fig.add_subplot(1,1,i+1) for i in range(1)]
    for _c in c:
        ax[0].plot(n, c, "-o", color="black", markersize=3)
    #ax[0].plot(n, c_2)
    ax[0].plot([hat_n, hat_n], [0, 1])
    ax[0].set_ylim([0, 1])
    plt.show()





if __name__ == "__main__":
    main()