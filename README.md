# Computational Algebra: Finite Fields, Polynomial Factorization & Fast Arithmetic

A from-scratch implementation of finite field arithmetic, polynomial
factorization over finite fields, and fast (sub-quadratic) polynomial
arithmetic — coursework for *Computational Algebra* (Mathematics-Physics
double degree, Universidad Complutense de Madrid).

Everything is built from first principles in pure Python: no external
algebra libraries (SymPy, SageMath, etc.) are used — every field operation,
factorization routine and fast-arithmetic algorithm implements the
underlying mathematics directly.

## Contents

### `cuerpos_finitos.py` — finite fields and factorization

- **`cuerpo_fp`** — the prime field 𝔽ₚ = ℤ/pℤ: arithmetic, binary
  exponentiation, multiplicative inverses via Fermat's little theorem,
  addition/multiplication tables, Latin squares.
- **`anillo_fp_x`** — the polynomial ring 𝔽ₚ[x]: division with remainder,
  (extended) GCD, and full **polynomial factorization via
  Cantor–Zassenhaus**:
  - *SFD* (square-free decomposition)
  - *DDF* (distinct-degree factorization)
  - *EDF* (equal-degree factorization, via random splitting)
  - **Rabin's irreducibility test**, used as a fast primality-style check
    before attempting factorization
- **`cuerpo_fq`** — extension fields 𝔽_q = 𝔽ₚ[x]/⟨g(x)⟩ for irreducible g,
  with the same operation set plus element enumeration and Latin squares.
- **`anillo_fq_x`** — polynomials over 𝔽_q: the same SFD/DDF/EDF pipeline
  adapted to extension-field coefficients.

### `algoritmos_rapidos.py` — fast arithmetic

Sub-quadratic algorithms layered on top of the field arithmetic above,
over both 𝔽ₚ[x] and 𝔽_q[x]:

- **Karatsuba multiplication** — O(n^log₂3) polynomial multiplication,
  replacing the naive O(n²) convolution.
- **Toeplitz matrix–vector products** and **Toeplitz matrix inversion**
  via Newton iteration on formal power series (doubling precision each
  step) — used to express polynomial division as a Toeplitz system and
  solve it without repeated long division.
- **Fast Fourier Transform over finite fields** (Cooley–Tukey radix-2 NTT)
  and its inverse, using an element of order n = 2^k in 𝔽ₚ* (or 𝔽_q*) as
  the transform's root of unity.

## Why this matters

Fields and their extensions are the algebraic backbone of coding theory
and cryptography (Reed–Solomon codes, AES's 𝔽₂⁸, elliptic-curve
cryptography over 𝔽_p). Polynomial factorization over finite fields
(Cantor–Zassenhaus) is a core primitive there and in computer algebra
systems. The fast-arithmetic half is the more general point: the same
Karatsuba/Toeplitz/NTT machinery that speeds up polynomial arithmetic
here is exactly what underlies fast large-integer multiplication,
polynomial-time algorithms in computational number theory, and FFT-based
convolution more broadly — the number-theoretic transform is the
finite-field analogue of the FFT used throughout scientific computing.

## Usage

```python
import cuerpos_finitos as cf

Fp = cf.cuerpo_fp(7)                    # F_7
Fpx = cf.anillo_fp_x(Fp)                # F_7[x]

f = Fpx.elem_de_str("x^4 + 1")
print(Fpx.es_irreducible(f))            # Rabin's test
print(Fpx.factorizar(f))                # Cantor-Zassenhaus factorization
```

```python
import algoritmos_rapidos as fast

# Fast multiplication in F_p[x] via Karatsuba
h = fast.fp_x_mult_karatsuba(Fpx, f, f)

# Number-theoretic transform (Cooley-Tukey) over F_p
# g: element of order n = 2^k in F_p*
y = fast.fp_fft(Fp, g, k, a)
```

## Requirements

Pure Python (standard library only — `random`, `re`, `itertools`, `math`,
`typing`).
