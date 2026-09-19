import random
import math
from typing import List, Tuple, Any, Optional

# =========================
# Utilidades generales
# =========================

def _is_prime(p: int) -> bool:
    if p <= 1:
        return False
    if p <= 3:
        return True
    if p % 2 == 0 or p % 3 == 0:
        return False
    i = 5
    while i * i <= p:
        if p % i == 0 or p % (i + 2) == 0:
            return False
        i += 6
    return True


def _strip_trailing_zeros(coefs: List[int]) -> Tuple[int, ...]:
    i = len(coefs) - 1
    while i >= 0 and coefs[i] == 0:
        i -= 1
    if i < 0:
        return tuple()
    return tuple(coefs[:i+1])

def _strip_trailing_objs(coefs: List[Any], is_zero_fn) -> Tuple[Any, ...]:
    i = len(coefs) - 1
    while i >= 0 and is_zero_fn(coefs[i]):
        i -= 1
    if i < 0:
        return tuple()
    return tuple(coefs[:i+1])

def _int_to_base_digits(n: int, base: int) -> List[int]:
    if n < 0:
        raise ValueError("El entero debe ser no negativo en esta conversión.")
    if n == 0:
        return [0]
    ds = []
    while n:
        ds.append(n % base)
        n //= base
    return ds  # coef 0,1,2,...

def _base_digits_to_int(digs: List[int], base: int) -> int:
    v = 0
    powb = 1
    for d in digs:
        v += d * powb
        powb *= base
    return v

def _parse_poly_over_integers(s: str, var: str) -> List[int]:
    """
    Parser simple de polinomios con coeficientes enteros:
    ejemplos: "0", "1", "x", "3x^2+2x+1", "-x^3+2", "x^5-10x+7"
    Soporta espacios y signos.
    Devuelve lista de coeficientes enteros (sin mod).
    """
    s = s.replace(' ', '')
    if s == '' or s == '+':
        return [0]
    # normaliza signos para split
    s = s.replace('-', '+-')
    if s.startswith('+-'):
        s = s[1:]
    parts = s.split('+')
    coefs = {}
    for term in parts:
        if term == '' or term == '+':
            continue
        if var in term:
            # algo como "3x^2", "-x^3", "x", "5x"
            before, _, after = term.partition(var)
            # coeficiente
            if before in ('', '+'):
                c = 1
            elif before == '-':
                c = -1
            elif before.endswith('*'):
                before = before[:-1]
                c = int(before) if before not in ('', '+', '-') else (1 if before in ('', '+') else -1)
            else:
                c = int(before)
            # exponente
            if after.startswith('^'):
                k = int(after[1:])
            elif after == '':
                k = 1
            else:
                raise ValueError(f"Término mal formado: {term}")
        else:
            # constante
            c = int(term)
            k = 0
        coefs[k] = coefs.get(k, 0) + c
    if not coefs:
        return [0]
    deg = max(coefs.keys())
    res = [0] * (deg + 1)
    for k, c in coefs.items():
        res[k] = c
    return res

# =========================
# 𝔽_p
# =========================

class cuerpo_fp:
    class _Elem:
        __slots__ = ('_v',)
        def __init__(self, v: int):
            self._v = v
        def __repr__(self):
            return f"FpElem({self._v})"

    def __init__(self, p): # construye el cuerpo de p elementos Fp = Z/pZ
        if not _is_prime(p):
            raise ValueError("p debe ser primo para que Z/pZ sea un cuerpo.")
        self.p = int(p)
        self._rng = random.Random()

    def _mk(self, v: int):
        return cuerpo_fp._Elem(v % self.p)

    # básicos
    def cero(self):
        return self._mk(0)

    def uno(self):
        return self._mk(1)

    def elem_de_int(self, n):
        return self._mk(int(n))

    def elem_de_str(self, s):
        s = s.strip()
        if s == '':
            s = '0'
        return self._mk(int(s))

    def conv_a_int(self, a):
        return a._v

    def conv_a_str(self, a):
        return str(a._v)

    # aritmética
    def suma(self, a, b):
        return self._mk(a._v + b._v)

    def inv_adit(self, a):
        return self._mk(-a._v)

    def mult(self, a, b):
        return self._mk(a._v * b._v)

    def pot(self, a, k):
        k = int(k)
        if k < 0:
            return self.pot(self.inv_mult(a), -k)
        return self._mk(pow(a._v, k, self.p))

    def inv_mult(self, a):
        if a._v == 0:
            raise ZeroDivisionError("0 no tiene inverso multiplicativo en 𝔽_p.")
        return self._mk(pow(a._v, self.p - 2, self.p))  # Fermat

    # predicados
    def es_cero(self, a):
        return a._v == 0

    def es_uno(self, a):
        return a._v == 1

    def es_igual(self, a, b):
        return a._v == b._v

    # aleatorio y tablas
    def aleatorio(self):
        return self._mk(self._rng.randrange(self.p))

    def tabla_suma(self):
        M = []
        for i in range(self.p):
            row = []
            for j in range(self.p):
                row.append((i + j) % self.p)
            M.append(row)
        return M

    def tabla_mult(self):
        M = []
        for i in range(self.p):
            row = []
            for j in range(self.p):
                row.append((i * j) % self.p)
            M.append(row)
        return M

    def tabla_inv_adit(self):
        return [(-i) % self.p for i in range(self.p)]

    def tabla_inv_mult(self):
        res = ['*'] * self.p  # índice 0 '*'
        for i in range(1, self.p):
            res[i] = pow(i, self.p - 2, self.p)
        return res

    def cuadrado_latino(self, a):
        if self.es_cero(a):
            raise ValueError("a debe ser no nulo.")
        A = self.conv_a_int(a)
        M = []
        for i in range(self.p):
            row = []
            for j in range(self.p):
                val = (A * i + j) % self.p
                row.append(val)
            M.append(row)
        return M


# =========================
# 𝔽_p[x]
# =========================

class anillo_fp_x:
    class _Poly:
        __slots__ = ('_c',)  # tupla de coeficientes en 0..p-1 (sin ceros de cola)
        def __init__(self, coefs: Tuple[int, ...]):
            self._c = coefs
        def __repr__(self):
            return f"FpPoly{self._c}"

    def __init__(self, fp: cuerpo_fp, var='x'):
        self.fp = fp
        self.p = fp.p
        self.var = var
        self._rng = random.Random()

    def _mk(self, coefs: List[int]):
        coefs = [c % self.p for c in coefs]
        return anillo_fp_x._Poly(_strip_trailing_zeros(coefs))

    def _zero_tuple(self):
        return tuple()

    # básicos
    def cero(self):
        return self._mk([])

    def uno(self):
        return self._mk([1])

    def elem_de_tuple(self, a):
        # admite ints o elementos de 𝔽_p
        lst = []
        for c in a:
            if isinstance(c, cuerpo_fp._Elem):
                lst.append(self.fp.conv_a_int(c))
            else:
                lst.append(int(c))
        return self._mk(lst)

    def elem_de_int(self, a):
        digs = _int_to_base_digits(int(a), self.p)
        return self._mk(digs)

    def elem_de_str(self, s):
        coefs = _parse_poly_over_integers(s, self.var)
        return self._mk(coefs)

    def conv_a_tuple(self, a):
        return a._c

    def conv_a_int(self, a):
        return _base_digits_to_int(list(a._c), self.p)

    def conv_a_str(self, a):
        c = a._c
        if not c:
            return "0"
        terms = []
        for i, ai in enumerate(c):
            if ai == 0:
                continue
            coef = ai
            if i == 0:
                terms.append(str(coef))
            elif i == 1:
                if coef == 1:
                    terms.append(self.var)
                else:
                    terms.append(f"{coef}{self.var}")
            else:
                if coef == 1:
                    terms.append(f"{self.var}^{i}")
                else:
                    terms.append(f"{coef}{self.var}^{i}")
        if not terms:
            return "0"
        # combinar con +, pero los coef ya están en [0,p-1]
        return " + ".join(terms)

    # aritmética
    def suma(self, a, b):
        n = max(len(a._c), len(b._c))
        res = [0] * n
        for i in range(n):
            ai = a._c[i] if i < len(a._c) else 0
            bi = b._c[i] if i < len(b._c) else 0
            res[i] = (ai + bi) % self.p
        return self._mk(res)

    def inv_adit(self, a):
        return self._mk([(-ai) % self.p for ai in a._c])

    def mult(self, a, b):
        if not a._c or not b._c:
            return self.cero()
        n = len(a._c) + len(b._c) - 1
        res = [0] * n
        for i, ai in enumerate(a._c):
            if ai == 0:
                continue
            for j, bj in enumerate(b._c):
                if bj == 0:
                    continue
                res[i + j] = (res[i + j] + ai * bj) % self.p
        return self._mk(res)

    def mult_por_escalar(self, a, e):
        if isinstance(e, cuerpo_fp._Elem):
            e = self.fp.conv_a_int(e)
        e = int(e) % self.p
        return self._mk([(ai * e) % self.p for ai in a._c])

    def grado(self, a):
        return len(a._c) - 1

    # división
    def divmod(self, a, b):
        if not b._c:
            raise ZeroDivisionError("división por 0 en 𝔽_p[x].")
        a_co = list(a._c)
        b_co = b._c
        if not a_co:
            return (self.cero(), self.cero())
        m = len(a_co) - 1
        n = len(b_co) - 1
        if m < n:
            return (self.cero(), self._mk(a_co))
        inv_lc_b = pow(b_co[-1], self.p - 2, self.p)
        q = [0] * (m - n + 1)
        r = a_co[:]
        for k in range(m - n, -1, -1):
            coef = (r[n + k] * inv_lc_b) % self.p
            q[k] = coef
            if coef != 0:
                for j in range(n + 1):
                    r[j + k] = (r[j + k] - coef * b_co[j]) % self.p
        return (self._mk(q), self._mk(r))

    def div(self, a, b):
        q, _ = self.divmod(a, b)
        return q

    def mod(self, a, b):
        _, r = self.divmod(a, b)
        return r

    # gcd
    def gcd(self, a, b):
        A = a
        B = b
        while B._c:
            A, B = B, self.mod(A, B)
        # mónico
        if not A._c:
            return A
        inv_lc = pow(A._c[-1], self.p - 2, self.p)
        return self.mult_por_escalar(A, inv_lc)

    def gcd_ext(self, a, b):
        # EEA polinómico
        x0, x1 = self.uno(), self.cero()
        y0, y1 = self.cero(), self.uno()
        A, B = a, b
        while B._c:
            q, r = self.divmod(A, B)
            A, B = B, r
            x0, x1 = x1, self.suma(x0, self.inv_adit(self.mult(q, x1)))
            y0, y1 = y1, self.suma(y0, self.inv_adit(self.mult(q, y1)))
        # A = gcd(a,b). Hacerlo mónico
        if not A._c:
            return (A, self.cero(), self.cero())
        inv_lc = pow(A._c[-1], self.p - 2, self.p)
        g = self.mult_por_escalar(A, inv_lc)
        x = self.mult_por_escalar(x0, inv_lc)
        y = self.mult_por_escalar(y0, inv_lc)
        return (g, x, y)

    def inv_mod(self, a, b):
        g, x, _ = self.gcd_ext(a, b)
        if self.grado(g) != 0 or (g._c and g._c[0] != 1):
            raise ZeroDivisionError("a no es invertible módulo b.")
        return self.mod(x, b)

    def pot_mod(self, a, k, b):
        k = int(k)
        if k < 0:
            ai = self.inv_mod(a, b)
            return self.pot_mod(ai, -k, b)
        res = self.uno()
        base = self.mod(a, b)
        while k > 0:
            if k & 1:
                res = self.mod(self.mult(res, base), b)
            base = self.mod(self.mult(base, base), b)
            k >>= 1
        return res

    # tests & factorización (Rabin + Cantor–Zassenhaus)
    def _deriv(self, f):
        if not f._c:
            return self.cero()
        res = [0] * max(0, len(f._c) - 1)
        for i in range(1, len(f._c)):
            res[i - 1] = (f._c[i] * i) % self.p
        return self._mk(res)

    def es_cero(self, a):
        return len(a._c) == 0

    def es_uno(self, a):
        return a._c == (1,)

    def es_igual(self, a, b):
        return a._c == b._c

    def es_irreducible(self, f):
        # Rabin test sobre 𝔽_p, versión estándar y robusta
        n = self.grado(f)
        if n <= 0:
            return False
        if n == 1:
            return True

        # hacer f mónico
        if f._c[-1] != 1:
            inv_lc = pow(f._c[-1], self.p - 2, self.p)
            f = self.mult_por_escalar(f, inv_lc)

        # squarefree: gcd(f, f') = 1
        if not self.es_uno(self.gcd(f, self._deriv(f))):
            return False

    
        # helper: x en F_p[x]/(f)
        x = self.elem_de_tuple((0, 1))

        # factoriza n en primos
        def _prime_factors(m):
            pf = set()
            d = 2
            while d * d <= m:
                while m % d == 0:
                    pf.add(d)
                    m //= d
                d += 1 if d == 2 else 2  # 2,3,5,7,...
            if m > 1:
                pf.add(m)
            return sorted(pf)

        # Frobenius^k: eleva a p^k (iterando k veces la potencia p)
        def _frobenius_pow(a, k, modf):
            h = a
            for _ in range(k):
                h = self.pot_mod(h, self.p, modf)
            return h

        # Condición de Rabin: para cada primo r|n, gcd(x^{p^{n/r}} - x, f) = 1
        for r in _prime_factors(n):
            h = _frobenius_pow(x, n // r, f)
            g = self.gcd(self.suma(h, self.inv_adit(x)), f)  # h - x
            if self.grado(g) > 0:
                return False

        # Y además x^{p^n} ≡ x (mod f)
        h = _frobenius_pow(x, n, f)
        return self.es_igual(h, x)

    # --- Factorización Cantor–Zassenhaus ---
    def _squarefree_decomp(self, f):
        """
        Devuelve lista de (fi, mi) tal que f = ∏ fi^{mi}, cada fi libre de cuadrados.
        Maneja correctamente:
          - f' = 0  (f = g(x^p))
          - factores con multiplicidad múltiplo de p (recurre en la raíz p-ésima)
        """
        if self.es_cero(f):
            return []

        res = []
        df = self._deriv(f)
        p = self.p

        if self.es_cero(df):
            # f = g(x^p) con g en F_p[x]. Tomar raíz p-ésima y descomponer recursivamente.
            g_co = []
            for i in range(0, len(f._c), p):
                g_co.append(f._c[i])
            g = self._mk(g_co)
            sub = self._squarefree_decomp(g)
            for (h, m) in sub:
                res.append((h, m * p))
            return res

        g = self.gcd(f, df)
        w = self.div(f, g)
        i = 1
        # Parte coprima con f': produce factores con multiplicidad i=1,2,3,...
        while not self.es_uno(w):
            y = self.gcd(w, g)
            z = self.div(w, y)
            if not self.es_uno(z):
                res.append((z, i))
            w = y
            g = self.div(g, y)
            i += 1

        # Resto: multiplicidades múltiplo de p  -> tomar raíz p-ésima y multiplicidad m*p
        if not self.es_uno(g):
            g_co = []
            for k in range(0, len(g._c), p):
                g_co.append(g._c[k])
            g_root = self._mk(g_co)
            sub = self._squarefree_decomp(g_root)
            for (h, m) in sub:
                res.append((h, m * p))

        return res


    def _distinct_degree_decomp(self, f):
        """
        Para f squarefree, devuelve lista [(f1,1),...,(fd,1)] donde fi es el producto
        de irreducibles de grado i.
        """
        x = self.elem_de_tuple((0, 1))
        res = []
        R = f
        q = self.p
        i = 1
        h = x
        while 2 * i <= self.grado(R):
            h = self.pot_mod(h, q, R)
            g = self.gcd(self.suma(h, self.inv_adit(x)), R)  # h - x
            if not self.es_uno(g):
                res.append((g, i))
                R = self.div(R, g)
                h = self.mod(h, R) if R._c else self.cero()
            i += 1
        if R._c:
            res.append((R, self.grado(R)))
        return res

    def _random_poly(self, deg_max):
        # polinomio aleatorio de grado < deg_max
        if deg_max <= 0:
            return self.cero()
        co = [self._rng.randrange(self.p) for _ in range(deg_max)]
        # evitar cero
        if all(c == 0 for c in co):
            co[0] = 1
        return self._mk(co)

    def _edf_split(self, f, d):
        """
        Equal Degree Factorization: asume f squarefree y todos los irreducibles de grado d.
        Devuelve la lista de irreducibles (cada uno con multiplicidad 1).

        Robusto:
        - d == 1: determinista por gcd con (x + c)
        - d >= 2: Cantor–Zassenhaus con límite de reintentos y prueba con b±1
        """
        # --- Caso especial: d == 1 ---
        if d == 1:
            x = self.elem_de_tuple((0, 1))
            R = f
            out = []
            for c in range(self.p):
                if self.es_cero(R):
                    break
                Lc = self.suma(x, self.elem_de_tuple((c % self.p,)))  # x + c
                t = self.gcd(R, Lc)
                if not self.es_uno(t) and not self.es_cero(t):
                    out.append(t)
                    R = self.div(R, t)
            if R._c and self.grado(R) == 1:
                out.append(R)
            return out

        # --- Caso general: d >= 2 ---
        factors = [f]
        done = []
        qd = pow(self.p, d)
        target = (qd - 1) // 2

        while factors:
            g = factors.pop()
            if self.grado(g) == d:
                done.append(g)
                continue
            MAX_TRIES = 64
            tries = 0
            split_ok = False
            while tries < MAX_TRIES:
                tries += 1
                a = self._random_poly(self.grado(g))
                if not a._c:
                    continue
                b = self.pot_mod(a, target, g)  # a^{(q^d-1)/2} mod g
                # intentar con b - 1
                t = self.gcd(self.suma(b, self.inv_adit(self.uno())), g)
                if not self.es_uno(t) and not self.es_igual(t, g):
                    factors.append(t)
                    factors.append(self.div(g, t))
                    split_ok = True
                    break
                # intentar con b + 1
                t = self.gcd(self.suma(b, self.uno()), g)
                if not self.es_uno(t) and not self.es_igual(t, g):
                    factors.append(t)
                    factors.append(self.div(g, t))
                    split_ok = True
                    break
            if not split_ok:
                # Como fallback: lo damos por irreducible (suele ocurrir muy raramente)
                done.append(g)
        return done


# =========================
# 𝔽_q = 𝔽_p[a]/<g(a)>
# =========================

class cuerpo_fq:
    class _Elem:
        __slots__ = ('_c',)  # tupla de coeficientes en 0..p-1
        def __init__(self, coefs: Tuple[int, ...]):
            self._c = coefs
        def __repr__(self):
            return f"FqElem{self._c}"

    def __init__(self, fp: cuerpo_fp, g, var='a'):
        """
        Construye 𝔽_q = 𝔽_p[a]/<g(a)>.

        g DEBE ser una tupla/lista de coeficientes en 𝔽_p:
           g = (c0, c1, ..., cn)   con c0 término independiente.
        Cada ci puede ser int o cuerpo_fp._Elem.
        """
        self.fp = fp
        self.var = var
        self.Px = anillo_fp_x(fp, var=var)  # para manejar g y reducción

        # --- Solo se admite list/tuple de coeficientes en 𝔽_p ---
        if not isinstance(g, (list, tuple)):
            raise TypeError("g debe ser una tupla/lista de coeficientes en 𝔽_p (c0, c1, ..., cn).")

        # Normalizar coeficientes (admite ints o elementos de 𝔽_p); reduce mod p y quita ceros de cola
        coefs = []
        for c in g:
            if isinstance(c, cuerpo_fp._Elem):
                coefs.append(self.fp.conv_a_int(c))
            else:
                coefs.append(int(c))
        poly = self.Px._mk(coefs)

        # Validaciones
        if not poly._c:
            raise ValueError("g no puede ser el polinomio cero.")
        if self.Px.grado(poly) == 0:
            raise ValueError("g debe tener grado ≥ 1.")
        if not self.Px.es_irreducible(poly):
            raise ValueError("g debe ser irreducible para construir 𝔽_q.")

        self.g = poly
        self.m = self.Px.grado(poly)
        self.q = fp.p ** self.m
        self._rng = random.Random()

    def _reduce(self, coefs: List[int]):
        # reducir modulo g
        poly = self.Px._mk(coefs)
        _, r = self.Px.divmod(poly, self.g)
        # devolver tupla (long < m)
        return cuerpo_fq._Elem(_strip_trailing_zeros(list(r._c)))

    def cero(self):
        return self._reduce([])

    def uno(self):
        return self._reduce([1])

    def elem_de_tuple(self, a):
        lst = []
        for c in a:
            if isinstance(c, cuerpo_fp._Elem):
                lst.append(self.fp.conv_a_int(c))
            else:
                lst.append(int(c))
        return self._reduce(lst)

    def elem_de_int(self, a):
        digs = _int_to_base_digits(int(a), self.fp.p)
        return self._reduce(digs)

    def elem_de_str(self, s):
        # interpretar s como polinomio en 'a' sobre 𝔽_p
        poly = self.Px.elem_de_str(s)
        return self._reduce(list(poly._c))

    def conv_a_tuple(self, a):
        return a._c

    def conv_a_int(self, a):
        return _base_digits_to_int(list(a._c), self.fp.p)

    def conv_a_str(self, a):
        c = a._c
        if not c:
            return "0"
        terms = []
        for i, ai in enumerate(c):
            if ai == 0:
                continue
            if i == 0:
                terms.append(str(ai))
            elif i == 1:
                if ai == 1:
                    terms.append(self.var)
                else:
                    terms.append(f"{ai}{self.var}")
            else:
                if ai == 1:
                    terms.append(f"{self.var}^{i}")
                else:
                    terms.append(f"{ai}{self.var}^{i}")
        return " + ".join(terms) if terms else "0"

    # aritmética
    def suma(self, a, b):
        n = max(len(a._c), len(b._c))
        res = [0] * n
        for i in range(n):
            ai = a._c[i] if i < len(a._c) else 0
            bi = b._c[i] if i < len(b._c) else 0
            res[i] = (ai + bi) % self.fp.p
        return self._reduce(res)

    def inv_adit(self, a):
        return self._reduce([(-ai) % self.fp.p for ai in a._c])

    def mult(self, a, b):
        # multiplicación como polinomios y reducción
        if not a._c or not b._c:
            return self.cero()
        tmp = [0] * (len(a._c) + len(b._c) - 1)
        for i, ai in enumerate(a._c):
            if ai == 0:
                continue
            for j, bj in enumerate(b._c):
                if bj == 0:
                    continue
                tmp[i + j] = (tmp[i + j] + ai * bj) % self.fp.p
        return self._reduce(tmp)

    def pot(self, a, k):
        k = int(k)
        if k < 0:
            return self.pot(self.inv_mult(a), -k)
        res = self.uno()
        base = a
        while k > 0:
            if k & 1:
                res = self.mult(res, base)
            base = self.mult(base, base)
            k >>= 1
        return res

    def inv_mult(self, a):
        if not a._c:
            raise ZeroDivisionError("0 no es invertible en 𝔽_q.")
        # EEA sobre 𝔽_p[x] entre a(x) y g(x)
        ax = self.Px._mk(list(a._c))
        g = self.g
        G, X, _Y = self.Px.gcd_ext(ax, g)
        # G debe ser 1
        if not self.Px.es_uno(G):
            raise ZeroDivisionError("El elemento no es invertible (gcd != 1).")
        X = self.Px.mod(X, g)
        return self._reduce(list(X._c))

    # predicados
    def es_cero(self, a):
        return len(a._c) == 0

    def es_uno(self, a):
        return a._c == (1,)

    def es_igual(self, a, b):
        return a._c == b._c

    # aleatorio y tablas
    def aleatorio(self):
        co = [self._rng.randrange(self.fp.p) for _ in range(self.m)]
        return self._reduce(co)

    def tabla_suma(self):
        M = []
        for i in range(self.q):
            ai = self.elem_de_int(i)
            row = []
            for j in range(self.q):
                bj = self.elem_de_int(j)
                row.append(self.conv_a_int(self.suma(ai, bj)))
            M.append(row)
        return M

    def tabla_mult(self):
        M = []
        for i in range(self.q):
            ai = self.elem_de_int(i)
            row = []
            for j in range(self.q):
                bj = self.elem_de_int(j)
                row.append(self.conv_a_int(self.mult(ai, bj)))
            M.append(row)
        return M

    def tabla_inv_adit(self):
        return [self.conv_a_int(self.inv_adit(self.elem_de_int(i))) for i in range(self.q)]

    def tabla_inv_mult(self):
        res = ['*'] * self.q
        for i in range(1, self.q):
            a = self.elem_de_int(i)
            res[i] = self.conv_a_int(self.inv_mult(a))
        return res

    def cuadrado_latino(self, a):
        if self.es_cero(a):
            raise ValueError("a debe ser no nulo.")
        M = []
        for i in range(self.q):
            ei = self.elem_de_int(i)
            row = []
            for j in range(self.q):
                ej = self.elem_de_int(j)
                val = self.suma(self.mult(a, ei), ej)
                row.append(self.conv_a_int(val))
            M.append(row)
        return M


# =========================
# 𝔽_q[x]
# =========================

class anillo_fq_x:
    class _Poly:
        __slots__ = ('_c',)  # tupla de coeficientes (objetos de 𝔽_q), sin ceros de cola
        def __init__(self, coefs: Tuple[Any, ...]):
            self._c = coefs
        def __repr__(self):
            return f"FqPoly{self._c}"

    def __init__(self, fq: cuerpo_fq, var='x'):
        if var == fq.var:
            raise ValueError("La variable del anillo y la del cuerpo no deben coincidir.")
        self.fq = fq
        self.var = var
        self.q = fq.q
        self._rng = random.Random()

    def _mk(self, coefs: List[Any]):
        # coefs son elementos de 𝔽_q (objetos). Normaliza ceros de cola.
        return anillo_fq_x._Poly(_strip_trailing_objs(coefs, self.fq.es_cero))

    # básicos
    def cero(self):
        return self._mk([])

    def uno(self):
        return self._mk([self.fq.uno()])

    def elem_de_tuple(self, a):
        lst = []
        for c in a:
            if isinstance(c, cuerpo_fq._Elem):
                lst.append(c)
            else:
                # intenta como entero en 𝔽_q
                lst.append(self.fq.elem_de_int(int(c)))
        return self._mk(lst)

    def elem_de_int(self, a):
        # descompone en base q (coeficientes en 0..q-1) y los mapea a 𝔽_q
        digs = _int_to_base_digits(int(a), self.q)
        coefs = [self.fq.elem_de_int(d) for d in digs]
        return self._mk(coefs)

    def _parse_list_coeffs(self, s: str):
        # parsing del tipo "[coef0, coef1, ...]" con coefi strings de 𝔽_q
        t = s.strip()
        if not (t.startswith('[') and t.endswith(']')):
            raise ValueError("Se esperaba una lista de coeficientes entre corchetes.")
        inner = t[1:-1].strip()
        if inner == '':
            return []
        parts = [u.strip() for u in inner.split(',')]
        coefs = [self.fq.elem_de_str(u) for u in parts]
        return coefs

    def elem_de_str(self, s):
        t = s.strip()
        if t.startswith('['):
            return self._mk(self._parse_list_coeffs(t))
        # fallback muy simple: polinomio con coeficientes enteros
        # (si quieres coeficientes generales, usa la forma de lista)
        coefs_int = _parse_poly_over_integers(t, self.var)
        coefs = [self.fq.elem_de_int(ci % self.fq.fp.p) for ci in coefs_int]
        return self._mk(coefs)

    def conv_a_tuple(self, a):
        # devuelve la tupla de coeficientes (objetos de 𝔽_q)
        return a._c

    def conv_a_int(self, a):
        # convierte cada coef a int base p y después evalúa en base q
        digs = [self.fq.conv_a_int(ci) for ci in a._c]
        return _base_digits_to_int(digs, self.q)

    def conv_a_str(self, a):
        c = a._c
        if not c:
            return "0"
        terms = []
        for i, ai in enumerate(c):
            if self.fq.es_cero(ai):
                continue
            cs = self.fq.conv_a_str(ai)
            if i == 0:
                terms.append(cs)
            elif i == 1:
                if self.fq.es_uno(ai):
                    terms.append(self.var)
                else:
                    terms.append(f"({cs}){self.var}")
            else:
                if self.fq.es_uno(ai):
                    terms.append(f"{self.var}^{i}")
                else:
                    terms.append(f"({cs}){self.var}^{i}")
        return " + ".join(terms) if terms else "0"

    # aritmética
    def suma(self, a, b):
        n = max(len(a._c), len(b._c))
        res = []
        for i in range(n):
            ai = a._c[i] if i < len(a._c) else self.fq.cero()
            bi = b._c[i] if i < len(b._c) else self.fq.cero()
            res.append(self.fq.suma(ai, bi))
        return self._mk(res)

    def inv_adit(self, a):
        return self._mk([self.fq.inv_adit(ai) for ai in a._c])

    def mult(self, a, b):
        if not a._c or not b._c:
            return self.cero()
        n = len(a._c) + len(b._c) - 1
        res = [self.fq.cero() for _ in range(n)]
        for i, ai in enumerate(a._c):
            if self.fq.es_cero(ai):
                continue
            for j, bj in enumerate(b._c):
                if self.fq.es_cero(bj):
                    continue
                res[i + j] = self.fq.suma(res[i + j], self.fq.mult(ai, bj))
        return self._mk(res)

    def mult_por_escalar(self, a, e):
        if not isinstance(e, cuerpo_fq._Elem):
            e = self.fq.elem_de_int(int(e))
        return self._mk([self.fq.mult(ai, e) for ai in a._c])

    def grado(self, a):
        return len(a._c) - 1

    # división
    def divmod(self, a, b):
        if not b._c:
            raise ZeroDivisionError("división por 0 en 𝔽_q[x].")
        a_co = list(a._c)
        b_co = list(b._c)
        if not a_co:
            return (self.cero(), self.cero())
        m = len(a_co) - 1
        n = len(b_co) - 1
        if m < n:
            return (self.cero(), self._mk(a_co))
        inv_lc_b = self.fq.inv_mult(b_co[-1])
        q = [self.fq.cero() for _ in range(m - n + 1)]
        r = a_co[:]
        for k in range(m - n, -1, -1):
            coef = self.fq.mult(r[n + k], inv_lc_b)
            q[k] = coef
            if not self.fq.es_cero(coef):
                for j in range(n + 1):
                    r[j + k] = self.fq.suma(r[j + k], self.fq.inv_adit(self.fq.mult(coef, b_co[j])))
        return (self._mk(q), self._mk(r))

    def div(self, a, b):
        q, _ = self.divmod(a, b)
        return q

    def mod(self, a, b):
        _, r = self.divmod(a, b)
        return r

    # gcd
    def gcd(self, a, b):
        A, B = a, b
        while B._c:
            A, B = B, self.mod(A, B)
        if not A._c:
            return A
        inv_lc = self.fq.inv_mult(A._c[-1])
        return self.mult_por_escalar(A, inv_lc)

    def gcd_ext(self, a, b):
        x0, x1 = self.uno(), self.cero()
        y0, y1 = self.cero(), self.uno()
        A, B = a, b
        while B._c:
            q, r = self.divmod(A, B)
            A, B = B, r
            x0, x1 = x1, self.suma(x0, self.inv_adit(self.mult(q, x1)))
            y0, y1 = y1, self.suma(y0, self.inv_adit(self.mult(q, y1)))
        if not A._c:
            return (A, self.cero(), self.cero())
        inv_lc = self.fq.inv_mult(A._c[-1])
        g = self.mult_por_escalar(A, inv_lc)
        x = self.mult_por_escalar(x0, inv_lc)
        y = self.mult_por_escalar(y0, inv_lc)
        return (g, x, y)

    def inv_mod(self, a, b):
        g, x, _ = self.gcd_ext(a, b)
        if self.grado(g) != 0 or not self.fq.es_uno(g._c[0]):
            raise ZeroDivisionError("a no es invertible módulo b.")
        return self.mod(x, b)

    def pot_mod(self, a, k, b):
        k = int(k)
        if k < 0:
            ai = self.inv_mod(a, b)
            return self.pot_mod(ai, -k, b)
        res = self.uno()
        base = self.mod(a, b)
        while k > 0:
            if k & 1:
                res = self.mod(self.mult(res, base), b)
            base = self.mod(self.mult(base, base), b)
            k >>= 1
        return res

    # predicados
    def es_cero(self, a):
        return len(a._c) == 0

    def es_uno(self, a):
        return len(a._c) == 1 and self.fq.es_uno(a._c[0])

    def es_igual(self, a, b):
    # Igualdad coeficiente a coeficiente en Fq
      if len(a._c) != len(b._c):
          return False
      for ai, bi in zip(a._c, b._c):
          if not self.fq.es_igual(ai, bi):
              return False
      return True

    # irreducibilidad y factorización sobre 𝔽_q
    def _deriv(self, f):
        if not f._c:
            return self.cero()
        res = [self.fq.cero() for _ in range(max(0, len(f._c) - 1))]
        for i in range(1, len(f._c)):
            # derivada: i * a_i, pero i está en Z -> mapeamos i mod p dentro de F_q (entero en F_q)
            # Como la característica es p, i puede anularse cuando p|i.
            ai = f._c[i]
            if i % self.fq.fp.p == 0:
                coef = self.fq.cero()
            else:
                coef_int = (i % self.fq.fp.p)
                coef = self.fq.mult(ai, self.fq.elem_de_int(coef_int))
            res[i - 1] = coef
        return self._mk(res)

    def es_irreducible(self, f):
        # Rabin en 𝔽_q, con q = p^m
        n = self.grado(f)
        if n <= 0:
            return False
        if n == 1:
            return True

        # mónico
        if not self.fq.es_uno(f._c[-1]):
            inv_lc = self.fq.inv_mult(f._c[-1])
            f = self.mult_por_escalar(f, inv_lc)

        # squarefree: gcd(f, f') = 1
        if not self.es_uno(self.gcd(f, self._deriv(f))):
            return False

        # x en 𝔽_q[x]/(f)
        x = self.elem_de_tuple((self.fq.cero(), self.fq.uno()))

        p = self.fq.fp.p
        q = self.q  # = p^m

        def _prime_factors(m):
            pf = set()
            d = 2
            while d * d <= m:
                while m % d == 0:
                    pf.add(d)
                    m //= d
                d += 1 if d == 2 else 2
            if m > 1:
                pf.add(m)
            return sorted(pf)

        # Frobenius_q^k: eleva a q^k (iterando k veces la potencia q)
        def _frobenius_q_pow(a, k, modf):
            h = a
            for _ in range(k):
                h = self.pot_mod(h, q, modf)
            return h

        # Condición de Rabin sobre 𝔽_q: para r|n primo, gcd(x^{q^{n/r}} - x, f) = 1
        for r in _prime_factors(n):
            h = _frobenius_q_pow(x, n // r, f)
            g = self.gcd(self.suma(h, self.inv_adit(x)), f)  # h - x
            if self.grado(g) > 0:
                return False

        # y x^{q^n} ≡ x (mod f)
        h = _frobenius_q_pow(x, n, f)
        return self.es_igual(h, x)

    def _random_poly(self, deg_max):
        if deg_max <= 0:
            return self.cero()
        co = [self.fq.aleatorio() for _ in range(deg_max)]
        # evitar cero
        if all(self.fq.es_cero(c) for c in co):
            co[0] = self.fq.uno()
        return self._mk(co)

    def _squarefree_decomp(self, f):
        """
        Devuelve lista de (fi, mi) tal que f = ∏ fi^{mi}, cada fi libre de cuadrados.
        Maneja correctamente el caso f' = 0 (f = g(x^p)).
        """
        if self.es_cero(f):
            return []
        res = []
        df = self._deriv(f)
        p = self.fq.fp.p

        if self.es_cero(df):
            # f = g(x^p). Tomar raíz p-ésima (en el índice) y repetir.
            step = p
            g_co = []
            for i in range(0, len(f._c), step):
                g_co.append(f._c[i])
            g = self._mk(g_co)
            for (h, m) in self._squarefree_decomp(g):
                res.append((h, m * step))
            return res

        g = self.gcd(f, df)
        w = self.div(f, g)
        i = 1
        while not self.es_uno(w):
            y = self.gcd(w, g)
            z = self.div(w, y)
            if not self.es_uno(z):
                res.append((z, i))
            w = y
            g = self.div(g, y)
            i += 1

        if not self.es_uno(g):
            # Los factores que quedan tienen multiplicidad múltiplo de p
            # quitamos raíz p-ésima repetidamente acumulando multiplicidades
            step = p
            while not self.es_uno(g):
                g_co = []
                for k in range(0, len(g._c), step):
                    g_co.append(g._c[k])
                g = self._mk(g_co)
                res.append((g, i * step))
        return res

    def _distinct_degree_decomp(self, f):
        x = self.elem_de_tuple((self.fq.cero(), self.fq.uno()))
        res = []
        R = f
        q = self.q
        i = 1
        h = x
        while 2 * i <= self.grado(R):
            h = self.pot_mod(h, q, R)
            g = self.gcd(self.suma(h, self.inv_adit(x)), R)
            if not self.es_uno(g):
                res.append((g, i))
                R = self.div(R, g)
                h = self.mod(h, R) if R._c else self.cero()
            i += 1
        if R._c:
            res.append((R, self.grado(R)))
        return res

    def _edf_split(self, f, d):
        """
        Equal Degree Factorization: asume que f es squarefree y que TODOS sus
        factores irreducibles tienen grado d. Devuelve la lista de dichos
        irreducibles (cada uno con multiplicidad 1).

        Robusto:
        - d == 1: determinista por búsqueda de raíces/gcd con (x + c)
        - d >= 2: Cantor–Zassenhaus con límite de reintentos
        """
        # --- Caso especial: d == 1 (lineales) ---
        if d == 1:
            x = self.elem_de_tuple((self.fq.cero(), self.fq.uno()))
            R = f
            out = []
            # Recorremos todos los c in F_q: lineales (x + c)
            for i in range(self.q):
                if self.es_cero(R):
                    break
                c = self.fq.elem_de_int(i)
                Lc = self.suma(x, self.elem_de_tuple((c,)))  # x + c
                # Extraer todas las copias (en principio squarefree ⇒ a lo sumo 1)
                t = self.gcd(R, Lc)
                if not self.es_uno(t) and not self.es_cero(t):
                    out.append(t)
                    R = self.div(R, t)
            # Si quedó un lineal suelto (por colisiones o construcción), añádelo
            if R._c and self.grado(R) == 1:
                out.append(R)
            return out

        # --- Caso general: d >= 2 ---
        factors = [f]
        done = []
        qd = pow(self.q, d)
        target = (qd - 1) // 2

        while factors:
            g = factors.pop()
            if self.grado(g) == d:
                done.append(g)
                continue
            # Limitar reintentos con distintos 'a'
            MAX_TRIES = 64
            tries = 0
            split_ok = False
            while tries < MAX_TRIES:
                tries += 1
                a = self._random_poly(self.grado(g))
                if self.es_cero(a):
                    continue
                b = self.pot_mod(a, target, g)         # a^{(q^d-1)/2} mod g
                # Primero con b - 1
                t = self.gcd(self.suma(b, self.inv_adit(self.uno())), g)
                if not self.es_uno(t) and not self.es_igual(t, g):
                    g1 = t
                    g2 = self.div(g, t)
                    factors.append(g1)
                    factors.append(g2)
                    split_ok = True
                    break
                # Luego probamos también con b + 1 (técnica clásica)
                t = self.gcd(self.suma(b, self.uno()), g)
                if not self.es_uno(t) and not self.es_igual(t, g):
                    g1 = t
                    g2 = self.div(g, t)
                    factors.append(g1)
                    factors.append(g2)
                    split_ok = True
                    break
            if not split_ok:
                # Como último recurso: devuelve g tal cual (debería ser ya irreducible)
                # o reintenta con otra semilla aleatoria si prefieres.
                done.append(g)
        return done
