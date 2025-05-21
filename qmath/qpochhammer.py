import mpmath as mp


class QPochhammer:
    def __init__(self, q):
        self.q = mp.mpf(q)
        self.lnq = mp.log(q)
        self._q_pow = {}  # Для мемоизации степеней q
        self._pochhammer_cache = {}  # Кэш для значений (a; q)_n

    def q_pow(self, k):
        """Вычисляет q^k с мемоизацией."""
        if k not in self._q_pow:
            self._q_pow[k] = mp.exp(k * self.lnq)
        return self._q_pow[k]

    def evaluate_pochhammer(self, a, n):
        """
            (a; q)_n = произведение (1 - a*q^k) для k от 0 до n-1
        """
        if n == 0:
            return mp.mpf(1)

        result = mp.mpf(1)
        for k in range(n):
            result *= (1 - a * self.q_pow(k))
        return result

    def __call__(self, a, n):
        a_mp = mp.mpf(a)
        cache_key = (a_mp, n)

        if cache_key not in self._pochhammer_cache:
            self._pochhammer_cache[cache_key] = self.evaluate_pochhammer(a_mp, n)

        return self._pochhammer_cache[cache_key]