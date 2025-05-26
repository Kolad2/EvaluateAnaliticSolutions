import mpmath as mp
from qmath import QBinom, QPochhammer


def c_tr_pow_limited(n, q=0.5, max_terms=50):
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
    qpoch_inf = qpoch(q, 41) # (q;q)_\infty
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