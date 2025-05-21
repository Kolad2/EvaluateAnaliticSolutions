import mpmath as mp


class QBinom:
    def __init__(self, q):
        self.q = mp.mpf(q)
        self.lnq = mp.log(q)
        self._q_pow = {}
        self._rows = [[mp.mpf(1)]]  # Инициализируем первую строку [1]

    def q_pow(self, m):
        """Вычисляет q^m с мемоизацией."""
        if m not in self._q_pow:
            self._q_pow[m] = mp.exp(m * self.lnq)
        return self._q_pow[m]

    def evaluate_row(self, n):
        """
        Возвращает строку q-биномиальных коэффициентов для степени n.
        Использует рекуррентную формулу:
        \binom{n}{i}_q = \binom{n-1}{i}_q + q^{n-i}*\binom{n-1}{i-1}_q
        """
        upper_row = self.get_row(n - 1)
        row = [mp.mpf(0)] * (n + 1)
        row[0] = mp.mpf(1)
        row[n] = mp.mpf(1)
        for i in range(1, n):
            row[i] = upper_row[i] + self.q_pow(n - i) * upper_row[i - 1]
        return row

    def get_row(self, n):
        while n >= len(self._rows):
            self._rows.append(self.evaluate_row(len(self._rows)))
        return self._rows[n]

    def __call__(self, n, k):
        if k < 0 or k > n:
            return mp.mpf(0)
        return self.get_row(n)[k]