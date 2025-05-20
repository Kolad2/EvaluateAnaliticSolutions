import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp


class QPoch:
    def __init__(self, q):
        self.q = q


class QBinom:
    def __init__(self, q):
        self.q = mp.mpf(q)
        self.lnq = mp.log(q)
        self._res = []
        self._q_pow = {}

    def q_pow(self, m):
        """Вычисляет q^m с мемоизацией."""
        if m not in self._q_pow:
            if m == 0:
                self._q_pows[m] = mp.mpf(1)
            else:
                self._q_pow[m] = mp.exp(m * self.lnq)
        return self._q_pow[m]

    def q_binom_n_1(self, n):
        if self._res[n] is not None:
            return

    def get_row(self, m):
        """Возвращает строку q-биномиальных коэффициентов для степени m.

        Использует рекуррентную формулу:
        B(m,k)_q = B(m-1,k)_q + q^(m-k)*B(m-1,k-1)_q
        """
        # Если строка уже вычислена - возвращаем её
        if m < len(self._rows):
            return self._rows[m]

        # Вычисляем все предыдущие строки, если нужно
        for i in range(len(self._rows), m + 1):
            prev_row = self._rows[i - 1]
            current_row = [mp.mpf(1)]  # Первый элемент = 1

            for k in range(1, i):
                term = prev_row[k] + self.q_pow(i - k) * prev_row[k - 1]
                current_row.append(term)

            current_row.append(mp.mpf(1))  # Последний элемент = 1
            self._rows.append(current_row)

        return self._rows[m]

    def __call__(self, n, k):
        if k <= 0:
            return 0
        if k == 1:
            return self.q_binom_n_1(n)
        return self.lnq


def main():
    q = 0.5
    binom = QBinom(q)
    binom(5, 6)

    # x = np.arange(1,5)
    # y = x^2
    # plt.figure()
    # plt.plot(x,y)
    # plt.show()


if __name__ == '__main__':
    main()