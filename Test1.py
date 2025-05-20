import numpy as np
import matplotlib.pyplot as plt

def main():
    x = np.arange(1,5)
    y = x^2
    plt.figure()
    plt.plot(x,y)
    plt.show()


if __name__ == '__main__':
    main()