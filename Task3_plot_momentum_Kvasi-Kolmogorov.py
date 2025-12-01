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
    q = 0.92

    eta_2, eta_3 = evaluate_moments(t, q)
    print(eta_2, eta_3)


def evaluate_moments(t=1, q=0.9):
    n_max = 60
    n = np.arange(0, n_max)
    c = c_tr_pow(n_max, t, q)
    
    if not np.all(c > 0):
        Warning("Wrong Solution, c < 0")
    
    n_0 = np.sum(c)
    pdf = c / n_0
    c = None
    m_1 = np.sum(n * pdf)
    eta_2 = np.sqrt(np.sum((n - m_1) ** 2 * pdf))
    eta_3 = np.cbrt(np.sum((n - m_1) ** 3 * pdf))

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(n, pdf)
    plt.show()

    return eta_2, eta_3



if __name__ == "__main__":
    main()