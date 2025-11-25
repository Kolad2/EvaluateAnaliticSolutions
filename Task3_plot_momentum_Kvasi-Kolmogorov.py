import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from qmath import QBinom, QPochhammer
from solutions import c_tr_pow, c_tr_const

def main():



    n_max = 60
    #t = np.array([25.0, 25.0*q**(-15), 25.0*q**(-30), 25.0*q**(-45)])
    lt = np.array([20, 50, 100])

    #%%
    t = 10
    q = 0.91

    eta_2, eta_3 = evaluate_moments(t, q)
    print(eta_2, eta_3)


def evaluate_moments(t=1, q=0.9):
    n_max = 60
    n = np.arange(0, n_max)
    c_q = c_tr_pow(n_max, t, q)
    n_0 = np.sum(c_q)
    m_1 = np.sum(n * c_q / n_0)
    eta_2 = np.sqrt(np.sum((n - m_1) ** 2 * c_q / n_0))
    eta_3 = np.cbrt(np.sum((n - m_1) ** 3 * c_q / n_0))

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(n,c_q)
    plt.show()

    return eta_2, eta_3



if __name__ == "__main__":
    main()