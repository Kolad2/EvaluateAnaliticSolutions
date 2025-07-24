import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from qmath import QBinom, QPochhammer
from solutions import c_tr_pow, c_tr_const, c_pow_pow


def main():
    n_max = 60
    n = np.arange(0, n_max)
    
    y = c_pow_pow(n, 0.5, 1, 0.5)
    
    fig = plt.figure(figsize=(10, 5))
    ax = [fig.add_subplot(1, 1, 1)]
    ax[0].plot(n, y)
    plt.show()


if __name__ == "__main__":
    main()
    