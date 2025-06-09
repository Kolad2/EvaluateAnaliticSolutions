import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from numpy.ma.core import where
from scipy.special import erf
from qmath import QBinom, QPochhammer
from solutions import c_tr_pow, c_tr_const
from scipy.optimize import minimize

def main():
    q = 0.9
    t = 0
    n_max = 60
    t = np.array([25.0*q**(-30)])
    lt = np.array([20, 50, 100])
    n = np.arange(0, n_max)
    #%%
    #c_q = [c_tr_pow(n_max, _t, q) for _t in t]

    fig = plt.figure(figsize=(10, 5))
    ax = [fig.add_subplot(1,1,i+1) for i in range(1)]

    for _t in t:
        cdf_cq = get_trpow(n_max, _t, q)

        def functional(x):
            return np.max(np.abs(cdf_cq - get_erf(n, x[0], x[1])))

        theta = minimize(
            functional,
            np.array([40, 1]),
            bounds=((1e-3, None), (1e-3, None)),
            method='nelder-mead',
            tol=1e-3
        ).x
        lognorm = get_erf(n, theta[0], theta[1])
        ax[0].plot(-n, cdf_cq, "-o", color="black", markersize=6,linewidth=3)
    ax[0].plot(-n, lognorm, "--", color="grey",linewidth=3)
    ax[0].set_ylim([-0.01, 1.01])
    n_max_10 = np.ceil(n_max * np.log(2) / np.log(10))
    xticks_10 = np.arange(0, -n_max_10 - 1, -2)
    xticks = xticks_10 * (np.log(10) / np.log(2))
    ax[0].set_xticks(xticks)
    xticklabels = [f'$10^{{{int(l)}}}$' for l in xticks_10]
    ax[0].set_xticklabels(xticklabels)
    ax[0].set_xlim([-50, -35])
    #ax[0].set_xlim([-40, -30])
    plt.show()

def get_erf(n, mu=0, sgm=1):
    return (1-erf((n-mu)/sgm*np.sqrt(2)))*(1/2)


def get_trpow(n_max, pt, q):
    n = np.arange(0, n_max)
    res = c_tr_pow(n_max, pt, q)
    res = np.where(res > 0, res, 0)
    res = np.cumsum(res * (2 ** n))
    res = 1 - res / res[-1]
    return res

if __name__ == "__main__":
    main()