import numpy as np
import matplotlib.pyplot as plt


def main():
    mu1, sigma1, mu2, sigma2, mu3, sigma3 = load_data_mu_sigma()
    a_s, sigma_s = load_data_a_sigma()
    K = np.mean(a_s * ((sigma_s) ** 2)) * (np.log(2)) ** 2
    print(K, np.sqrt(K))
    exit()
    
    fig = plt.figure(figsize=(8.5, 4))
    
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(mu1 * np.log(2), sigma1 * np.log(2), color='black', linestyle="-", label="0.15")
    ax1.plot(mu2 * np.log(2), sigma2 * np.log(2), color='black', linestyle="--", label="0.25")
    ax1.plot(mu3 * np.log(2), sigma3 * np.log(2), color='black', linestyle="-.", label="0.50")
    ax1.set_xlim([2, 25])
    ax1.set_xlabel(r"$\mu_e$")
    ax1.set_ylabel(r"$\sigma_e$")
    ax1.legend(loc='lower right', ncol=3)
    ax1.set_title("а)")
    
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(a_s, sigma_s * np.log(2), color='black', linestyle="-")
    ax2.set_xlabel(r"$\alpha$")
    ax2.set_ylabel(r"$\sigma_e$")
    ax2.set_title("б)")
    
    ax1.set_ylim(ax2.get_ylim())
    
    plt.savefig('graphs/1.eps', format='eps', dpi=300)
    plt.show()


def load_data_a_sigma():
    matrix = np.load("temp/a_sigma.npy")
    a_s = matrix[:, 0]
    sigma_s = matrix[:, 1]
    return a_s, sigma_s


def load_data_mu_sigma():
    matrix = np.load('temp/sigmas.npy')
    mu1 = matrix[:, 0]
    sigma1 = matrix[:, 1]
    mu2 = matrix[:, 2]
    sigma2 = matrix[:, 3]
    mu3 = matrix[:, 4]
    sigma3 = matrix[:, 5]
    return mu1, sigma1, mu2, sigma2, mu3, sigma3


if __name__ == '__main__':
    main()