import numpy as np
import matplotlib.pyplot as plt
from qmath import QBinom, QPochhammer
import mpmath as mp
from scipy.special import gamma

def main():
    q = 0.93
    t = 0
    t_max = 40
    #
    qpoch = QPochhammer(q)
    qpoch_inf = qpoch(q, 15)  # (q;q)_\infty
    #
    t = np.arange(0, t_max, 0.1)[1:]
    c1 = np.array([c_tr_pow_limited(_t, 0.90) for _t in t])
    c2 = np.array([c_tr_pow_limited(_t, 0.91) for _t in t])

    #
    fig = plt.figure(figsize=(10, 5))
    ax = [fig.add_subplot(1,1,i+1) for i in range(1)]

    ax[0].plot(t, np.log(c1), color="black")
    ax[0].plot(t, np.log(c2), color="black")
    d = 19
    alpha = 1.5
    c1_1 = alpha ** (d + 1) * (t ** d) * np.exp(-alpha * t) / gamma(d + 1)
    ax[0].plot(t,  np.log(c1_1))
    d = 20
    alpha = 1.9
    ax[0].plot(t, alpha**(d + 1)*(t ** d)*np.exp(-alpha * t) / gamma(d + 1))
    #ax[0].set_ylim([0, 1])
    plt.show()


def c_tr_pow_limited(t, q=0.5, max_terms=20):
    r"""
    Вычисляет N_n по формуле:
    c_n = \frac{2^n}{(q;q)_{\infty}} \sum_{j=0}^{\infty} (-1)^j \frac{q^{j(j-1)/2}}{(q;q)_j} e^{-q^{n-j}}

    Параметры:
    - n: индекс
    - q: параметр q (0 < q < 1)
    - max_terms: максимальное количество слагаемых в сумме (для аппроксимации бесконечного ряда)

    Возвращает:
    - Значение N_n
    """
    q = mp.mpf(q)
    logq = mp.log(q)
    t = mp.mpf(t)

    qpoch = QPochhammer(q)
    qpoch_inf = qpoch(q, 41) # (q;q)_\infty
    # Вычисляем сумму
    res = mp.mpf(0)
    sign = 1
    for j in range(max_terms):
        q_poch_j = qpoch(q, j)  # (q;q)_j
        exponent = mp.exp(-mp.power(q,-j)*t + j * (j - 1) / 2 * logq)
        term = sign * exponent / q_poch_j
        res = res + term
        sign = (-1) * sign

    res = res / qpoch_inf
    return float(res)



if __name__ == "__main__":
    main()