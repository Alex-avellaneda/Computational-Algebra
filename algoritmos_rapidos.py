# ============================================================
# Funciones auxiliares
# ============================================================

import cuerpos_finitos as cf

# --- Karatsuba para 𝔽_p[x] (coeficientes enteros mod p) ---

def _karatsuba_fp_coefs(p, A, B, THRESH=32):
    """
    Multiplica dos listas de coeficientes A, B (enteros) mod p
    usando Karatsuba y devuelve una lista de enteros (ya mod p).
    """
    if not A or not B:
        return []
    if min(len(A), len(B)) <= THRESH:
        n = len(A) + len(B) - 1
        res = [0] * n
        for i, ai in enumerate(A):
            if ai == 0:
                continue
            for j, bj in enumerate(B):
                if bj == 0:
                    continue
                res[i + j] = (res[i + j] + ai * bj) % p
        return res

    m = min(len(A), len(B)) // 2
    A0 = A[:m]
    A1 = A[m:]
    B0 = B[:m]
    B1 = B[m:]

    Z0 = _karatsuba_fp_coefs(p, A0, B0, THRESH)
    Z2 = _karatsuba_fp_coefs(p, A1, B1, THRESH)

    # (A0 + A1), (B0 + B1)
    maxA = max(len(A0), len(A1))
    Ap = [(A0[i] if i < len(A0) else 0) + (A1[i] if i < len(A1) else 0)
          for i in range(maxA)]
    maxB = max(len(B0), len(B1))
    Bp = [(B0[i] if i < len(B0) else 0) + (B1[i] if i < len(B1) else 0)
          for i in range(maxB)]
    Ap = [x % p for x in Ap]
    Bp = [x % p for x in Bp]

    Z1 = _karatsuba_fp_coefs(p, Ap, Bp, THRESH)

    # Z1 = Z1 - Z0 - Z2
    maxZ = max(len(Z1), len(Z0), len(Z2))
    ZZ1 = [0] * maxZ
    for i in range(len(Z1)):
        ZZ1[i] = (ZZ1[i] + Z1[i]) % p
    for i in range(len(Z0)):
        ZZ1[i] = (ZZ1[i] - Z0[i]) % p
    for i in range(len(Z2)):
        ZZ1[i] = (ZZ1[i] - Z2[i]) % p
    Z1 = [x % p for x in ZZ1]

    n = len(A) + len(B) - 1
    res = [0] * n
    for i, v in enumerate(Z0):
        res[i] = (res[i] + v) % p
    for i, v in enumerate(Z1):
        if i + m < n:
            res[i + m] = (res[i + m] + v) % p
    for i, v in enumerate(Z2):
        if i + 2 * m < n:
            res[i + 2 * m] = (res[i + 2 * m] + v) % p
    return res


# --- Karatsuba para 𝔽_q[x] (coeficientes objetos de 𝔽_q) ---

def _karatsuba_fq_coefs(fq, A, B, THRESH=16):
    """
    Multiplica dos listas de coeficientes A, B (elementos de 𝔽_q)
    usando Karatsuba y devuelve una lista de elementos de 𝔽_q.
    """
    if not A or not B:
        return []
    if min(len(A), len(B)) <= THRESH:
        n = len(A) + len(B) - 1
        res = [fq.cero() for _ in range(n)]
        for i, ai in enumerate(A):
            if fq.es_cero(ai):
                continue
            for j, bj in enumerate(B):
                if fq.es_cero(bj):
                    continue
                res[i + j] = fq.suma(res[i + j], fq.mult(ai, bj))
        return res

    m = min(len(A), len(B)) // 2
    A0 = A[:m]
    A1 = A[m:]
    B0 = B[:m]
    B1 = B[m:]

    Z0 = _karatsuba_fq_coefs(fq, A0, B0, THRESH)
    Z2 = _karatsuba_fq_coefs(fq, A1, B1, THRESH)

    def _add_lists(X, Y):
        n = max(len(X), len(Y))
        out = []
        for i in range(n):
            xi = X[i] if i < len(X) else fq.cero()
            yi = Y[i] if i < len(Y) else fq.cero()
            out.append(fq.suma(xi, yi))
        return out

    Ap = _add_lists(A0, A1)
    Bp = _add_lists(B0, B1)

    Z1 = _karatsuba_fq_coefs(fq, Ap, Bp, THRESH)

    # Z1 = Z1 - Z0 - Z2
    nZ = max(len(Z1), len(Z0), len(Z2))
    ZZ1 = [fq.cero() for _ in range(nZ)]
    for i in range(len(Z1)):
        ZZ1[i] = fq.suma(ZZ1[i], Z1[i])
    for i in range(len(Z0)):
        ZZ1[i] = fq.suma(ZZ1[i], fq.inv_adit(Z0[i]))
    for i in range(len(Z2)):
        ZZ1[i] = fq.suma(ZZ1[i], fq.inv_adit(Z2[i]))
    Z1 = ZZ1

    n = len(A) + len(B) - 1
    res = [fq.cero() for _ in range(n)]
    for i, v in enumerate(Z0):
        res[i] = fq.suma(res[i], v)
    for i, v in enumerate(Z1):
        if i + m < n:
            res[i + m] = fq.suma(res[i + m], v)
    for i, v in enumerate(Z2):
        if i + 2 * m < n:
            res[i + 2 * m] = fq.suma(res[i + 2 * m], v)
    return res


# ============================================================
# 𝔽_p[x]: Karatsuba, Toeplitz, división, FFT
# ============================================================

# input: fpx -> anillo_fp_x
# input: f -> polinomio (objeto opaco creado por fpx)
# input: g -> polinomio (objeto opaco creado por fpx)
# output: f*g calculado usando el método de Karatsuba
def fp_x_mult_karatsuba(fpx, f, g):
    A = list(fpx.conv_a_tuple(f))
    B = list(fpx.conv_a_tuple(g))
    if not A or not B:
        return fpx.cero()
    coefs = _karatsuba_fp_coefs(fpx.p, A, B)
    return fpx.elem_de_tuple(coefs)

# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fp_x.mult_fast = fp_x_mult_karatsuba


# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera columna de una
#    matriz de Toeplitz inferior T de nxn)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_inf_vec(fp, n, a, b):
    """
    T Toeplitz inferior nxn con primera columna a = (a0,...,a_{n-1}),
    b = vector de longitud n.
    Devuelve T*b.
    Se implementa como convolución A(x)*B(x), donde
        A(x) = sum a_i x^i, B(x) = sum b_i x^i,
    usando _karatsuba_fp_coefs sobre enteros mod p.
    """
    p = fp.p

    # pasamos de elementos de Fp a enteros mod p
    A_int = [fp.conv_a_int(x) % p for x in a]
    B_int = [fp.conv_a_int(x) % p for x in b]

    # convolución
    C_int = _karatsuba_fp_coefs(p, A_int, B_int)

    # truncamos/rellenamos a longitud n
    if len(C_int) < n:
        C_int += [0] * (n - len(C_int))
    else:
        C_int = C_int[:n]

    # devolvemos como elementos de Fp
    return tuple(fp.elem_de_int(ci) for ci in C_int)


# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera fila de una
#    matriz de Toeplitz superior T de nxn)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_sup_vec(fp, n, a, b):
    """
    Multiplica una Toeplitz superior T por un vector b.
    a[k] = T[0,k] = coef de la primera fila.
    Formula:
        (Tb)_i = sum_{j=i..n-1} a[j-i] * b[j]
    """
    res = []
    for i in range(n):
        acc = fp.cero()
        for j in range(i, n):
            acc = fp.suma(acc, fp.mult(a[j - i], b[j]))
        res.append(acc)
    return tuple(res)


# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud 2*n-1 de elementos de fp (primera fila de una
#    matriz de Toeplitz completa T de nxn seguida de la primera columna
#    excepto el elemento de la esquina)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_vec(fp, n, a, b):
    """
    Multiplica una Toeplitz completa T (n×n) por un vector b.

    La Toeplitz se define a partir del vector a de longitud 2n-1:

        a[0..n-1]     = primera fila  (T[0,j])
        a[n..2n-2]    = primera columna sin la esquina (T[i,0], i>=1)

    Reconstruimos:
        row[j] = T[0,j]                  para j = 0..n-1
        col[i] = T[i,0]                  para i = 0..n-1   (añadimos row[0] al inicio)

    Fórmula de Toeplitz completa:
        T[i,j] = row[j-i]   si j >= i
               = col[i-j]   si j <  i

    Devuelve un vector de Fp de longitud n.
    """

    # reconstruimos primera fila y primera columna completas
    row = a[:n]                # (row[0], row[1], ..., row[n-1])
    col = (row[0],) + a[n:]    # (row[0], col[1], ..., col[n-1])

    res = []
    for i in range(n):
        acc = fp.cero()
        for j in range(n):
            if j >= i:
                # parte superior / diagonal
                acc = fp.suma(acc, fp.mult(row[j - i], b[j]))
            else:
                # parte inferior
                acc = fp.suma(acc, fp.mult(col[i - j], b[j]))
        res.append(acc)

    return tuple(res)


# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera columna de una
#    matriz de Toeplitz inferior T de nxn)... suponemos a[0] != 0
# output: primera columna de T^(-1) -> tupla de longitud n de elementos de
#    fp (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz inferior
def fp_toep_inf_inv(fp, n, a):
    """
    Inversa Toeplitz inferior en Fp usando inverso de serie formal (Newton).
    a = (a0,...)
    Devuelve la primera columna de T^{-1}.
    """
    p = fp.p

    # Convertir a enteros mod p
    a_int = [(fp.conv_a_int(x) % p) for x in a]

    a0 = a_int[0] % p
    if a0 == 0:
        raise ZeroDivisionError("a[0] = 0: Toeplitz no invertible")

    # b = [a0^{-1}] (entero mod p)
    b = [pow(a0, p - 2, p)]
    m = 1

    while m < n:
        # AB = A * B (convolución truncada)
        AB = _karatsuba_fp_coefs(p, a_int, b)

        L = 2 * m
        if len(AB) < L:
            AB += [0] * (L - len(AB))
        else:
            AB = AB[:L]

        # C = 2 - AB
        C = []
        for i, c in enumerate(AB):
            if i == 0:
                C.append((2 - c) % p)
            else:
                C.append((-c) % p)

        # B_new = B*C (truncado)
        B_new = _karatsuba_fp_coefs(p, b, C)
        if len(B_new) < L:
            B_new += [0] * (L - len(B_new))
        else:
            B_new = B_new[:L]

        b = B_new
        m *= 2

    b = b[:n]

    # devolver como elementos Fp
    return tuple(fp.elem_de_int(x % p) for x in b)


# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera fila de una
#    matriz de Toeplitz superior T de nxn)... suponemos a[0] != 0
# output: primera fila de T^(-1) -> tupla de longitud n de elementos de
#    fp (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz superior
def fp_toep_sup_inv(fp, n, a):
    """
    Inversa de Toeplitz superior. Usamos la identidad
    T_sup(a) = J T_inf(a) J, de modo que T_sup^{-1} tiene como
    primera fila la misma tupla que la primera columna de T_inf^{-1}.
    """
    col = fp_toep_inf_inv(fp, n, a)
    return col


# input: fpx -> anillo_fp_x
# input: f -> polinomio (objeto opaco creado por fpx)
# input: g -> polinomio no nulo (objeto opaco creado por fpx)
# output: q -> cociente
# output: r -> resto
# se cumple que f = g*q+r, r=0 o deg(r)<deg(g)
# reformular el problema en términos de matrices de Toeplitz y luego usar
# las funciones de arriba para obtener q y r
def fp_x_divmod(fpx, f, g):
    fp = fpx.fp

    # Aceptar también tuplas/listas como atajo
    if isinstance(f, (tuple, list)):
        f = fpx.elem_de_tuple(tuple(f))
    if isinstance(g, (tuple, list)):
        g = fpx.elem_de_tuple(tuple(g))

    # Obtener coeficientes ascendentes (enteros mod p)
    f_co = list(fpx.conv_a_tuple(f))
    g_co = list(fpx.conv_a_tuple(g))

    # Limpiar ceros finales
    while len(f_co) > 1 and f_co[-1] == 0:
        f_co.pop()
    while len(g_co) > 1 and g_co[-1] == 0:
        g_co.pop()

    if not g_co:
        raise ZeroDivisionError("división por 0 en 𝔽_p[x]")
    if not f_co:
        return (fpx.cero(), fpx.cero())

    m = len(f_co) - 1   # grado de f
    n = len(g_co) - 1   # grado de g

    if m < n:
        # cociente 0, resto = f
        return (fpx.cero(), f)

    k = m - n + 1       # nº coeficientes del cociente

    # 1) Vector h = (f_m, f_{m-1}, ..., f_n)
    F_desc = [f_co[m - i] for i in range(k)]
    h = tuple(fp.elem_de_int(x) for x in F_desc)

    # 2) Columna Toeplitz: col[i] = g_{n-i} para i <= n, 0 si i>n
    G_desc = [g_co[n - i] for i in range(n + 1)]  # (g_n, ..., g_0)

    col = []
    for i in range(k):
        if i <= n:
            col.append(fp.elem_de_int(G_desc[i]))
        else:
            col.append(fp.cero())
    col = tuple(col)

    # 3) Resolver Toeplitz: T * Q_desc = h usando columna de T^{-1}
    Tinv_col = fp_toep_inf_inv(fp, k, col)
    Q_desc = fp_toep_inf_vec(fp, k, Tinv_col, h)

    # 4) Convertir Q_desc (descendente) → coef. ascendentes del cociente
    q_coefs = [
        int(Q_desc[k - 1 - i]._v) % fp.p
        for i in range(k)
    ]
    q = fpx.elem_de_tuple(tuple(q_coefs))

    # 5) Resto = f - g*q
    prod = fpx.mult(q, g)
    r = fpx.suma(f, fpx.inv_adit(prod))

    r_co = list(fpx.conv_a_tuple(r))
    while len(r_co) > 1 and r_co[-1] == 0:
        r_co.pop()
    r = fpx.elem_de_tuple(tuple(r_co))

    return (q, r)

# mantener el "método rápido" en la clase
cf.anillo_fp_x.divmod_fast = fp_x_divmod


# input: fp -> cuerpo_fp
# input: g -> elemento del grupo multiplicativo fp* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fp
# output: DFT_{n,g}(a) -> tupla de longitud n de elementos de fp
# utilizar el algoritmo de Cooley-Tuckey
def fp_fft(fp, g, k, a):
    n = 1 << k
    if n == 1:
        return (a[0],)
    if len(a) != n:
        raise ValueError("La longitud de a debe ser 2**k.")

    m = n // 2
    a_even = tuple(a[0::2])
    a_odd = tuple(a[1::2])

    g2 = fp.pot(g, 2)
    y_even = fp_fft(fp, g2, k - 1, a_even)
    y_odd = fp_fft(fp, g2, k - 1, a_odd)

    res = [None] * n
    w = fp.uno()
    for j in range(m):
        w_y_odd = fp.mult(w, y_odd[j])
        res[j] = fp.suma(y_even[j], w_y_odd)
        res[j + m] = fp.suma(y_even[j], fp.inv_adit(w_y_odd))
        w = fp.mult(w, g)
    return tuple(res)


# input: fp -> cuerpo_fp
# input: g -> elemento del grupo multiplicativo fp* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fp
# output: IDFT_{n,g}(a) -> tupla de longitud n de elementos de fp
# recordar que IDFT_{n,g} = n^(-1) * DFT_{n,g^(-1)}
def fp_ifft(fp, g, k, a):
    n = 1 << k
    if len(a) != n:
        raise ValueError("La longitud de a debe ser 2**k.")
    g_inv = fp.inv_mult(g)
    y = fp_fft(fp, g_inv, k, a)
    n_inv = fp.inv_mult(fp.elem_de_int(n))
    return tuple(fp.mult(n_inv, yi) for yi in y)


# ============================================================
# 𝔽_q[x]: Karatsuba, Toeplitz, división, FFT
# ============================================================

# input: fqx -> anillo_fq_x
# input: f -> polinomio (objeto opaco creado por fqx)
# input: g -> polinomio (objeto opaco creado por fqx)
# output: f*g calculado usando el método de Karatsuba
def fq_x_mult_karatsuba(fqx, f, g):
    fq = fqx.fq
    A = list(fqx.conv_a_tuple(f))
    B = list(fqx.conv_a_tuple(g))
    if not A or not B:
        return fqx.cero()
    coefs = _karatsuba_fq_coefs(fq, A, B)
    return fqx.elem_de_tuple(coefs)

# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fq_x.mult_fast = fq_x_mult_karatsuba


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera columna de una
#    matriz de Toeplitz inferior T de nxn)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_inf_vec(fq, n, a, b):
    # Interpretamos a y b como coeficientes polinómicos
    fqx = cf.anillo_fq_x(fq)  # usar variable distinta de fq.var

    A = fqx.elem_de_tuple(a)   # A(x) = sum a[i] x^i
    B = fqx.elem_de_tuple(b)   # B(x) = sum b[i] x^i

    # Convolución usando Karatsuba
    C = fq_x_mult_karatsuba(fqx, A, B)
    co = fqx.conv_a_tuple(C)

    # Resultado: (Tb)_i = c_i
    return tuple(co[i] if i < len(co) else fq.cero() for i in range(n))


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera fila de una
#    matriz de Toeplitz superior T de nxn)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_sup_vec(fq, n, a, b):
    # Implementación directa análoga a la versión sobre Fp
    res = []
    for i in range(n):
        acc = fq.cero()
        for j in range(i, n):
            acc = fq.suma(acc, fq.mult(a[j - i], b[j]))
        res.append(acc)
    return tuple(res)


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud 2*n-1 de elementos de fq (primera fila de una
#    matriz de Toeplitz completa T de nxn seguida de la primera columna
#    excepto el elemento de la esquina)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_vec(fq, n, a, b):
    # Igual que la versión sobre Fp, pero con operaciones en Fq
    row = a[:n]
    col = (row[0],) + a[n:]  # reconstruimos primera columna completa

    res = []
    for i in range(n):
        acc = fq.cero()
        for j in range(n):
            if j >= i:
                acc = fq.suma(acc, fq.mult(row[j - i], b[j]))
            else:
                acc = fq.suma(acc, fq.mult(col[i - j], b[j]))
        res.append(acc)
    return tuple(res)


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera columna de una
#    matriz de Toeplitz inferior T de nxn)... suponemos a[0] != 0
# output: primera columna de T^(-1) -> tupla de longitud n de elementos de
#    fq (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz inferior
def fq_toep_inf_inv(fq, n, a):
    """
    Devuelve la primera columna de T^{-1}, donde T es Toeplitz
    inferior n×n sobre Fq con primera columna a = (a0,...,a_{n-1}),
    y a0 != 0.

    Se usa el inverso de serie formal vía Newton:
        A(x) = sum_{i>=0} a_i x^i
        B(x) ≈ 1/A(x)  (mod x^n)

    Devuelve (b0, b1, ..., b_{n-1}) con B(x)=∑ b_i x^i.
    """

    if n <= 0:
        return tuple()

    # a0 debe ser invertible en Fq
    if fq.es_cero(a[0]):
        raise ZeroDivisionError("a[0] = 0: la Toeplitz no es invertible")

    # Convertimos a lista (los elementos son ya elementos de Fq)
    A = list(a)

    # B_0 = a0^{-1}
    B = [fq.inv_mult(a[0])]
    m = 1

    # IMPORTANTE: 2 como 1+1 en F_q, no elem_de_int(2), que en tu implementación es un polinomio en 'a'
    two = fq.suma(fq.uno(), fq.uno())

    # Newton: duplicar precisión hasta llegar a n
    while m < n:

        # Convolución AB (Karatsuba sobre Fq)
        AB = _karatsuba_fq_coefs(fq, A, B)

        # Truncar AB a longitud 2m
        L = 2*m
        if len(AB) < L:
            AB = AB + [fq.cero()] * (L - len(AB))
        else:
            AB = AB[:L]

        # C(x) = 2 - A(x)B(x)
        C = []
        for i, c in enumerate(AB):
            if i == 0:
                # 2 - c0
                C.append(fq.suma(two, fq.inv_adit(c)))
            else:
                # -ci
                C.append(fq.inv_adit(c))

        # B_new = B * C (mod x^{2m})
        B_new = _karatsuba_fq_coefs(fq, B, C)
        if len(B_new) < L:
            B_new = B_new + [fq.cero()] * (L - len(B_new))
        else:
            B_new = B_new[:L]

        B = B_new
        m *= 2

    # Truncar a longitud n exacta
    if len(B) < n:
        B = B + [fq.cero()] * (n - len(B))

    return tuple(B[:n])


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera fila de una
#    matriz de Toeplitz superior T de nxn)... suponemos a[0] != 0
# output: primera fila de T^(-1) -> tupla de longitud n de elementos de
#    fq (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz superior
def fq_toep_sup_inv(fq, n, a):
    """
    Análogo a fp_toep_sup_inv: la primera fila de T_sup^{-1}(a)
    coincide con la primera columna de T_inf^{-1}(a).
    """
    col = fq_toep_inf_inv(fq, n, a)
    return col


# input: fqx -> anillo_fq_x
# input: f -> polinomio (objeto opaco creado por fqx)
# input: g -> polinomio no nulo (objeto opaco creado por fqx)
# output: q -> cociente
# output: r -> resto
# se cumple que f = g*q+r, r=0 o deg(r)<deg(g)
# reformular el problema en términos de matrices de Toeplitz y luego usar
# las funciones de arriba para obtener q y r
def fq_x_divmod(fqx, f, g):
    fq = fqx.fq

    # Aceptar tuplas/listas como atajo
    if isinstance(f, (tuple, list)):
        f = fqx.elem_de_tuple(tuple(f))
    if isinstance(g, (tuple, list)):
        g = fqx.elem_de_tuple(tuple(g))

    # Coeficientes ascendentes (elementos de Fq)
    f_co = list(fqx.conv_a_tuple(f))
    g_co = list(fqx.conv_a_tuple(g))

    # Limpiar ceros finales
    while len(f_co) > 1 and fq.es_cero(f_co[-1]):
        f_co.pop()
    while len(g_co) > 1 and fq.es_cero(g_co[-1]):
        g_co.pop()

    if not g_co:
        raise ZeroDivisionError("división por 0 en 𝔽_q[x]")
    if not f_co:
        return (fqx.cero(), fqx.cero())

    m = len(f_co) - 1  # grado f
    n = len(g_co) - 1  # grado g

    if m < n:
        # cociente 0, resto = f
        return (fqx.cero(), f)

    k = m - n + 1  # nº coef. del cociente

    # 1) h = (f_m, f_{m-1}, ..., f_n)
    F_desc = [f_co[m - i] for i in range(k)]
    h = tuple(F_desc)

    # 2) Primera columna T: col[i] = g_{n-i}
    G_desc = [g_co[n - i] for i in range(n + 1)]  # (g_n,...,g_0)

    col = []
    for i in range(k):
        if i <= n:
            col.append(G_desc[i])
        else:
            col.append(fq.cero())
    col = tuple(col)

    # 3) Resolver Toeplitz inferior: T * Q_desc = h
    Tinv_col = fq_toep_inf_inv(fq, k, col)
    Q_desc = fq_toep_inf_vec(fq, k, Tinv_col, h)

    # 4) Convertir descendente → ascendente
    q_coefs = [Q_desc[k - 1 - i] for i in range(k)]
    q = fqx.elem_de_tuple(tuple(q_coefs))

    # 5) Resto r = f - g*q
    prod = fqx.mult(q, g)
    r = fqx.suma(f, fqx.inv_adit(prod))

    r_co = list(fqx.conv_a_tuple(r))
    while len(r_co) > 1 and fq.es_cero(r_co[-1]):
        r_co.pop()

    return (q, fqx.elem_de_tuple(tuple(r_co)))


cf.anillo_fq_x.divmod_fast = fq_x_divmod


# input: fq -> cuerpo_fq
# input: g -> elemento del grupo multiplicativo fq* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a q-1
# input: a -> tupla de longitud n de elementos de fq
# output: DFT_{n,g}(a) -> tupla de longitud n de elementos de fq
# utilizar el algoritmo de Cooley-Tuckey
def fq_fft(fq, g, k, a):
    n = 1 << k
    if n == 1:
        return (a[0],)
    if len(a) != n:
        raise ValueError("La longitud de a debe ser 2**k.")

    m = n // 2
    a_even = tuple(a[0::2])
    a_odd = tuple(a[1::2])

    g2 = fq.pot(g, 2)
    y_even = fq_fft(fq, g2, k - 1, a_even)
    y_odd = fq_fft(fq, g2, k - 1, a_odd)

    res = [None] * n
    w = fq.uno()
    for j in range(m):
        w_y_odd = fq.mult(w, y_odd[j])
        res[j] = fq.suma(y_even[j], w_y_odd)
        res[j + m] = fq.suma(y_even[j], fq.inv_adit(w_y_odd))
        w = fq.mult(w, g)
    return tuple(res)


# input: fq -> cuerpo_fq
# input: g -> elemento del grupo multiplicativo fq* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fq
# output: IDFT_{n,g}(a) -> tupla de longitud n de elementos de fq
# recordar que IDFT_{n,g} = n^(-1) * DFT_{n,g^(-1)}
def fq_ifft(fq, g, k, a):
    n = 1 << k
    if len(a) != n:
        raise ValueError("La longitud de a debe ser 2**k.")
    g_inv = fq.inv_mult(g)
    y = fq_fft(fq, g_inv, k, a)

    # Queremos n como elemento del subcuerpo primo F_p dentro de F_q:
    p = fq.fp.p
    n_mod_p = n % p                # n visto en Z/pZ
    scalar_n = fq.elem_de_int(n_mod_p)   # esto sí es constante en F_q
    n_inv = fq.inv_mult(scalar_n)

    return tuple(fq.mult(n_inv, yi) for yi in y)