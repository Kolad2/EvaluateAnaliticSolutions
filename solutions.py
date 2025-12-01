import mpmath as mp
from qmath import QBinom, QPochhammer
import numpy as np


def c_tr_const(n_max, t=1):
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


def c_tr_pow(n_max, t = 1, q = 0.5):
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
        for j in range(min(n+1, 10000)):
            term_1 = qbinom(n, j)
            term_2 = mp.exp(-mp.power(q, n - j)*t) # \exp(q^{n-j}*t)
            term_3 = mp.power(q, j * (j - 1) // 2) # q^{j(j-1)/2}
            term = term + sign*term_1 * term_2 * term_3
            sign = -sign
        c[n] = term / p
    return c


def c_tr_pow_limited(n, q=0.5, max_terms=20):
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
    n = mp.mpf(n)

    qpoch = QPochhammer(q)
    qpoch_inf = qpoch(q, 41)  # (q;q)_\infty
    # Вычисляем сумму
    res = mp.mpf(0)
    sign = 1
    for j in range(max_terms):
        q_poch_j = qpoch(q, j)  # (q;q)_j
        exponent = mp.exp(-mp.power(q, (n - j)) + j * (j - 1) / 2 * logq)
        term = sign * exponent / q_poch_j
        res = res + term
        sign = (-1) * sign

    res = res / qpoch_inf
    return float(res)


def c_pow_pow(n, a, b, t=1):
    gamma = (b + 2)/a
    res = 2.0 ** (-n*b) * (1 - 2.0 ** (-n*a)) ** gamma * np.exp((1-2.0**(-n*a))*t)
    return res
