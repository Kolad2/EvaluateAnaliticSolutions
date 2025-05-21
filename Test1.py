import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from qmath import QBinom

def main():
    q = 0.5
    binom = QBinom(q)
    print(binom(5, 2))  # Печатает q-биномиальный коэффициент C(5,2)_q


if __name__ == "__main__":
    main()