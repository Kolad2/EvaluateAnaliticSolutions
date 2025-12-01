import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from qmath import QBinom, QPochhammer
from solutions import c_tr_pow_limited, c_tr_const


def main():
    #generate_data()
    
    a_s, sigma_s = load_data()

    K = np.mean(a_s * sigma_s ** 2)
    print(K, np.sqrt(K))

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(a_s, sigma_s, color='black', linestyle="-")
    plt.savefig('graphs/2.eps', format='eps', dpi=300)
    plt.show()

def generate_data():
    # a_s = np.array([0.12, 0.20, 0.25, 0.5, 0.75])
    a_s = np.array([0.13, 0.15, 0.175, 0.20, 0.225, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75])
    sigma_s = np.empty(len(a_s))
    for i, a in enumerate(a_s):
        sigma_s[i] = get_sigma(a)
    export_mtrix = np.column_stack((a_s, sigma_s))
    np.save("temp/a_sigma.npy", export_mtrix)
    
    
def load_data():
    matrix = np.load("temp/a_sigma.npy")
    a_s = matrix[:, 0]
    sigma_s = matrix[:, 1]
    return a_s, sigma_s

def get_sigma(a=0.15):
    q = 2 ** (-a)
    print("a: ", a)
    ns = np.arange(-11, 12 + 1)
    
    n = 0
    cs = np.empty(len(ns))
    shift = np.log(a) / a - 1 / a
    #print(shift)
    for i, n in enumerate(ns):
        cs[i] = c_tr_pow_limited(n + shift, float(q), max_terms=200)
        if cs[i] < 0:
            cs[i] = 0
        cs[i] = cs[i] * (2 ** i)
    
    n_0 = np.sum(cs)
    pdf = cs / n_0
    cs = None
    mu = np.sum(ns * pdf)
    sigma = np.sqrt(np.sum((ns - mu) ** 2 * pdf))
    #
    # fig = plt.figure(figsize=(10, 5))
    # ax = fig.add_subplot(1, 1, 1)
    # ax.plot(ns, pdf)
    # plt.show()
    #
    return sigma


def evaluate_moments(t=1, q=0.9):
    
    
    h = 0.000001
    if not np.all(c > 0):
        Warning("Solution divergence, c < 0")
    c[c < h] = 0
    c = c * (2 ** n)
    
    n_0 = np.sum(c)
    pdf = c / n_0
    
    c = None
    mu = np.sum(n * pdf)
    sigma = np.sqrt(np.sum((n - mu) ** 2 * pdf))
    eta_3 = np.cbrt(np.sum((n - mu) ** 3 * pdf))
    
    
    
    return mu, sigma, eta_3


if __name__ == "__main__":
    main()