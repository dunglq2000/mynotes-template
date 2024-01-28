```python
from sage.all import *

def small_roots(f, M, X):
    d = f.degree()

    P = f.base_ring()
    B = []
    for i in range(d):
        b = [0] * (d + 1)
        b[i] = M * X**i
        B.append(b)
    B.append([j*X**i for i,j in enumerate(f.coefficients())])
    B = Matrix(ZZ, B)
    B = B.LLL()
    bf = B[0]
    g = sum((j//(X**i)) * x**i for i,j in enumerate(bf))
    roots = g.roots()
    return [root[0] for root in roots]

M = 10001
X = 10
P = PolynomialRing(ZZ, 'x')
Q = PolynomialRing(Zmod(M), 'xn')
x = P.gen()
xn = Q.gen()
f = x**3 + 10*x**2 + 5000*x - 222
g = f.change_ring(Q).subs(x=xn)
print(small_roots(f, M, X))
print(g.small_roots(X=10))
```

Coppersmith mở rộng

```python
from sage.all import *

def small_roots(f, M, X):
    d = f.degree()
    B = []
    for i in range(d):
        b = [0] * (2*d)
        b[i] = M * X**i
        B.append(b)
    for i in range(d):
        g = x**i * f
        b = [v * X**u for u, v in enumerate(g.coefficients(sparse=False))]
        b = b + [0] * (2*d - len(b))
        B.append(b)
    B = Matrix(ZZ, B)
    print(B)
    B = B.LLL()
    bf = B[0]
    g = sum((j // (X**i)) * x**i for i, j in enumerate(bf))
    roots = g.roots()
    return [root[0] for root in roots]

M = 10001
X = 10
P = PolynomialRing(ZZ, 'x')
x = P.gen()
f = x**3 + 10*x**2 + 5000*x - 222
print(small_roots(f, M, X))
```

Coppersmith nâng cao phần 2.

```python
from sage.all import *

def small_roots(f, M, h = None, epsilon = None, X = None):
    d = f.degree()
    if not h:
        h = d
    if not epsilon:
        epsilon = 1/(d*h)
    if not X:
        X = round(0.5*M**(1/d-epsilon))
    B = []
    for j in range(h):
        g = M**(h-1-j) * f**j
        for i in range(d):
            k = g * x**i
            b = [v * X**u for u, v in enumerate(k.coefficients(sparse=False))]
            b = b + [0] * (d*h - len(b))
            B.append(b)

    B = Matrix(ZZ, B)
    B = B.LLL()
    bf = B[0]
    g = sum((j // (X**i)) * x**i for i, j in enumerate(bf))
    roots = g.roots()
    return [root[0] for root in roots]

# Test theo sách giáo khoa

M = (2**30 + 3)*(2**32 + 15)
P = PolynomialRing(ZZ, 'x')
x = P.gen()
f = 1942528644709637042 + 1234567890123456789*x + 987654321987654321*x**2 + x**3
print(small_roots(f, M))

# Tự test

M = (2**20 + 7)*(2**21 + 17)
P = PolynomialRing(ZZ, 'x')
x = P.gen()
f = x**3 + (2**25 - 2883584)*x**2 + 46976195*x + 227
print(small_roots(f, M, X = 2**9))
```