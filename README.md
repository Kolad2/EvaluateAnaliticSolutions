# EvaluateAnaliticSolutions

```
pip install numpy
pip install matplotlib
pip install mpmath
pip install scipy
```

## Колмогоровское решение

$$
c_n(pt) = \frac{(2pt)^{n}}{n!}e^{-pt}
$$

## Степенная скорость

Вычисление функции

$$
c_{n}(pt) = \frac{2^n}{(q;q)_n}\sum_{j=0}^{n}
(-1)^{j}
\binom{n}{j}_qq^{j(j-1)/2}
e^{-q^{n-j}pt}
$$

Или в более компактном виде

$$
c_n(pt)=
q^{n(n-1)/2}
\frac{(2p_0t)^n}{[n]_q!}(-1)^nD_q^n[e^{-pt}]
$$

Связь через

$$
D_q^n e^{ax} = (-1)^n\frac{q^{-n(n-1)/2}}{(1 - q)^n x^{n}} \sum_{j=0}^n (-1)^{n -j} \binom{n}{j}_q q^{(n-j)(n-j-1)/2} e^{aq^j x}
$$

Аналитическое продожение, запись в виде ряда тейлора

$$
c_n(t)=
q^{n(n-1)/2}
(2t)^n
\sum_{j=0}^{\infty} (-1)^j \binom{j+n}{j}_q\frac{t^{j}}{(j + n)!} 
$$



---

### **Рекуррентные формулы для q-биномиальных коэффициентов**

#### **1. Основная рекуррента (аналог треугольника Паскаля)**

$$
\binom{n}{k}_q = \binom{n-1}{k}_q + q^{n-k} \cdot \binom{n-1}{k-1}_q
$$