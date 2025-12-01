import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from qmath import QBinom, QPochhammer
from solutions import c_tr_pow, c_tr_const

def main():
    # n_max = 60
    #
    # t = 300
    #
    # a1 = 0.15
    # a2 = 0.25
    # a3 = 0.50
    # mu1, sigma1, eta1 = get_mu_sigma_eta(a1, 10, num=30)
    # mu2, sigma2, eta2 = get_mu_sigma_eta(a2, 12, num=30)
    # mu3, sigma3, eta3 = get_mu_sigma_eta(a3, 25, num=30)
    #
    # export_matrix = np.column_stack((mu1, sigma1, mu2, sigma2, mu3, sigma3))
    # np.save('temp/sigmas.npy', export_matrix)
    mu1, sigma1, mu2, sigma2, mu3, sigma3 = load_data()
    
    fig = plt.figure(figsize=(5, 5))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(mu1, sigma1, color='black', linestyle="-", label="0.15")
    ax1.plot(mu2, sigma2, color='black', linestyle="--", label="0.25")
    ax1.plot(mu3, sigma3, color='black', linestyle="-.", label="0.50")
    ax1.set_xlim([2, 30])
    ax1.legend(loc='lower right', ncol=3)
    plt.savefig('graphs/1.eps', format='eps', dpi=300)
    plt.show()


def load_data():
    matrix = np.load('temp/sigmas.npy')
    mu1 = matrix[:, 0]
    sigma1 = matrix[:, 1]
    mu2 = matrix[:, 2]
    sigma2 = matrix[:, 3]
    mu3 = matrix[:, 4]
    sigma3 = matrix[:, 5]
    return mu1, sigma1, mu2, sigma2, mu3, sigma3


def get_mu_sigma_eta(a, stop=10, num=10):
    q = 2 ** (-a)
    ts = np.logspace(0, stop, num, base=2)
    mu = np.empty(len(ts))
    sigma = np.empty(len(ts))
    eta = np.empty(len(ts))
    for i, t in enumerate(ts):
        mu[i], sigma[i], eta[i] = evaluate_moments(float(t), q)
    return mu, sigma, eta


def evaluate_moments(t=1, q=0.9):
    n_max = 60
    n = np.arange(0, n_max)
    c = c_tr_pow(n_max, t, q)
    
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

    # fig = plt.figure(figsize=(10, 5))
    # ax = fig.add_subplot(1, 1, 1)
    # ax.plot(n, pdf)
    # plt.show()

    return mu, sigma, eta_3


if __name__ == "__main__":
    main()