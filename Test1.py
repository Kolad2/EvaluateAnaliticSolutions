import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from qmath import QBinom, QPochhammer



def main():
    q = 0.8
    t = 0
    n_max = 30
    t = np.array([5., 25, 125, 625])
    lt = np.array([5, 10, 25,])
    n = np.arange(0, n_max)
    #%%
    c_q = [get_qc(n_max, _t, q) for _t in t]
    #sigma_q = [np.sqrt(np.sum((n**2)*(_c - np.sum(n*_c)))) for _c in c_q]
    c = [get_c(n_max, 0.5*_t) for _t in lt]
    mu = [np.sum(n * _c) for _c in c]
    sigma = [np.sum((n**2)*(_c - np.sum(n*_c))) for _c in c]

    fig = plt.figure()
    ax = [fig.add_subplot(2,2,i+1) for i in range(4)]
    ax[0].set_xlabel("n")
    ax[1].set_xlabel("n")
    for _c in c:
        _c = np.where(_c >= 0.0001, _c, np.nan)
        _c = np.log(_c) + n*np.log(2)
        ax[0].plot(n, _c, "-o", color="black", markersize=3)
    #ax[0].set_ylim([np.log(0.001), 1.1])
    for _c in c_q:
        _c = np.where(_c >= 0.0001, _c, np.nan)
        _c = np.log(_c) + n*np.log(2)
        ax[1].plot(_c, "-o", color="black", markersize=3)
    #ax[1].set_ylim([-80, 1.1])

    ax[2].plot(lt, mu)
    plt.show()


def get_c(n_max, t=1):
    """
    Вычисляет c_n(pt) = (2pt)^n / n! * e^(-pt) для n от 0 до n_max-1.
    """
    t = mp.mpf(t)
    # Вычисляем общий множитель e^(-pt)
    exp_t = mp.exp(-t)
    # Инициализируем список результатов
    c = np.zeros(n_max)
    # Вычисляем c_n для каждого n
    for n in range(n_max):
        numerator = mp.power(t, n)  # (2pt)^n
        factorial_n = mp.factorial(n)  # n!
        c[n] = (numerator / factorial_n) * exp_t

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