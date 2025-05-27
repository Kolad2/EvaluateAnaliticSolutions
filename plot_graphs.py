import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from qmath import QBinom, QPochhammer
from solutions import c_tr_pow, c_tr_const

def main():
    q = 0.9
    t = 0
    n_max = 60
    t = np.array([25.0, 25.0*q**(-15), 25.0*q**(-30), 25.0*q**(-45)])
    lt = np.array([20, 50, 100])
    n = np.arange(0, n_max)
    #%%
    c_q = [c_tr_pow(n_max, _t, q) for _t in t]
    c = [c_tr_const(n_max, 0.5*_t) for _t in lt]

    fig = plt.figure(figsize=(10, 5))
    ax = [fig.add_subplot(1,2,i+1) for i in range(2)]
    for _c in c:
        _c = np.where(_c > 0, _c, np.nan)
        _c = np.log(_c) / np.log(2)
        _c = np.where(_c > -10, _c, np.nan) + n
        ax[0].plot(-n, _c, "-o", color="black", markersize=3)
    for _c in c_q:
        _c = np.where(_c > 0, _c, np.nan)
        _c = np.log(_c)/np.log(2)
        _c = np.where(_c > -10, _c, np.nan) + n
        ax[1].plot(-n,_c, "-o", color="black", markersize=3)

    ax[0].plot([0, -n_max], [-2, n_max - 2], "--", color="blue", markersize=3)
    ax[1].plot([0, -n_max], [-2, n_max - 2], "--", color="blue", markersize=3)

    n_max_10 = np.ceil(n_max*np.log(2)/np.log(10))
    xticks_10 = np.arange(0, -n_max_10-1, -4)
    xticks = xticks_10*(np.log(10)/np.log(2))
    ax[0].set_xticks(xticks)
    ax[1].set_xticks(xticks)
    xticklabels = [f'$10^{{{int(l)}}}$' for l in xticks_10]
    ax[0].set_xticklabels(xticklabels)
    ax[1].set_xticklabels(xticklabels)
    #
    n_max_10 = np.ceil(n_max * np.log(2) / np.log(10))
    yticks_10 = np.arange(0, n_max_10 - 1, 4)
    yticks = yticks_10 * (np.log(10) / np.log(2))
    ax[0].set_yticks(yticks)
    ax[1].set_yticks(yticks)
    yticklabels = [f'$10^{{{int(l)}}}$' for l in yticks_10]
    ax[0].set_yticklabels(yticklabels)
    ax[1].set_yticklabels(yticklabels)
    #
    ax[0].set_xlabel("Относительный размер, $x/x_0$")
    ax[1].set_xlabel("Относительный размер, $x/x_0$")
    ax[0].set_ylabel("Концентрация частиц, $c(x,t)/N_0$")
    ax[1].set_ylabel("Концентрация частиц, $c(x,t)/N_0$")
    plt.show()



if __name__ == "__main__":
    main()