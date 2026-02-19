"""
McEliece Cryptographic QUBO 생성기

Mandrà et al. (arXiv:2308.09704, FGCS 2025) 기반:
McEliece 암호 프로토콜의 공개키를 Ising 스핀 시스템으로 캐스팅하여,
암호학적으로 어려운 planted QUBO 인스턴스를 생성.

핵심 아이디어:
  1. Goppa 코드 C[N, k, d] 생성 (GF(2^m) 위)
  2. McEliece 키 생성: G' = S·G·P (공개키)
  3. 인코딩: q' = q·G' + ε (메시지 q에 에러 ε 추가)
  4. Ising 변환: 각 컬럼의 비영 행 인덱스로 p-local 스핀 상호작용 생성
  5. 차수 축소: p-body → 2-body (Rosenberg 보조변수)
  6. QUBO 변환: Ising → QUBO (s = 2x - 1)

장점:
  - 암호학적 보안에 기반한 computational hardness
  - ground state가 McEliece 복호화와 동치
  - 파라미터 (m, t)로 난이도 조절 가능
"""

import random
import itertools
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from qubo_utils import calculate_energy, save_qubo_edgelist, print_q_matrix, print_qubo_formula


# ============================================================
# 1. GF(2^m) 유한체 클래스
# ============================================================

# 기약다항식 테이블 (정수 표현): m → irreducible polynomial over GF(2)
# 예: m=4 → x^4 + x + 1 = 0b10011 = 19
IRREDUCIBLE_POLYS = {
    2: 0b111,        # x^2 + x + 1
    3: 0b1011,       # x^3 + x + 1
    4: 0b10011,      # x^4 + x + 1
    5: 0b100101,     # x^5 + x^2 + 1
    6: 0b1000011,    # x^6 + x + 1
    7: 0b10000011,   # x^7 + x + 1
    8: 0b100011101,  # x^8 + x^4 + x^3 + x^2 + 1
    9: 0b1000010001, # x^9 + x^4 + 1
    10: 0b10000001001, # x^10 + x^3 + 1
}


class GF2m:
    """GF(2^m) 유한체 산술 (log/antilog 룩업 테이블)."""

    def __init__(self, m):
        if m not in IRREDUCIBLE_POLYS:
            raise ValueError(f"m={m} 지원 안 함. 가능: {sorted(IRREDUCIBLE_POLYS.keys())}")
        self.m = m
        self.q = 1 << m          # 체의 원소 수 = 2^m
        self.mod_poly = IRREDUCIBLE_POLYS[m]
        self._build_tables()

    def _build_tables(self):
        """log/antilog 룩업 테이블 구축."""
        q = self.q
        self.exp_table = [0] * (2 * q)   # antilog: exp_table[i] = α^i
        self.log_table = [0] * q          # log: log_table[a] = i where α^i = a

        val = 1
        for i in range(q - 1):
            self.exp_table[i] = val
            self.log_table[val] = i
            val <<= 1
            if val >= q:
                val ^= self.mod_poly
        # α^(q-1) = 1 (순환)
        for i in range(q - 1, 2 * q):
            self.exp_table[i] = self.exp_table[i - (q - 1)]

    def add(self, a, b):
        """GF(2^m) 덧셈 = XOR."""
        return a ^ b

    def mul(self, a, b):
        """GF(2^m) 곱셈 (log 테이블 사용)."""
        if a == 0 or b == 0:
            return 0
        return self.exp_table[self.log_table[a] + self.log_table[b]]

    def inv(self, a):
        """GF(2^m) 역원: a^(-1) = a^(q-2)."""
        if a == 0:
            raise ZeroDivisionError("0의 역원은 없음")
        # α^i의 역원 = α^(q-1-i)
        return self.exp_table[(self.q - 1) - self.log_table[a]]

    def div(self, a, b):
        """GF(2^m) 나눗셈."""
        if b == 0:
            raise ZeroDivisionError("0으로 나눌 수 없음")
        if a == 0:
            return 0
        return self.exp_table[(self.log_table[a] - self.log_table[b]) % (self.q - 1)]

    def poly_eval(self, coeffs, x):
        """
        다항식 평가: coeffs[0] + coeffs[1]*x + ... + coeffs[d]*x^d.

        Horner's method 사용.
        """
        result = 0
        for c in reversed(coeffs):
            result = self.add(self.mul(result, x), c)
        return result

    def poly_mul(self, p1, p2):
        """두 다항식의 곱 (GF(2^m) 위)."""
        if not p1 or not p2:
            return [0]
        result = [0] * (len(p1) + len(p2) - 1)
        for i, a in enumerate(p1):
            for j, b in enumerate(p2):
                result[i + j] = self.add(result[i + j], self.mul(a, b))
        return result

    def poly_is_irreducible_gf2(self, poly_int, degree):
        """
        GF(2) 위의 다항식이 기약인지 검사.

        Ben-Or 알고리즘: x^(2^i) mod f(x) for i=1..degree//2.
        각 단계에서 gcd(x^(2^i) - x, f(x)) == 1이어야 함.
        """
        # GF(2) 위의 다항식 산술
        def gf2_poly_mod(a, m_poly, deg_m):
            """a mod m_poly (GF(2) 위)."""
            a = int(a)
            m_poly = int(m_poly)
            while a.bit_length() - 1 >= deg_m:
                a ^= m_poly << (a.bit_length() - 1 - deg_m)
            return a

        def gf2_poly_mulmod(a, b, m_poly, deg_m):
            """(a * b) mod m_poly (GF(2) 위)."""
            result = 0
            a = int(a)
            b = int(b)
            while b:
                if b & 1:
                    result ^= a
                a <<= 1
                b >>= 1
            return gf2_poly_mod(result, m_poly, deg_m)

        def gf2_poly_powmod(base, exp, m_poly, deg_m):
            """base^exp mod m_poly (GF(2) 위)."""
            result = 1  # = 1 (상수 다항식)
            base = gf2_poly_mod(base, m_poly, deg_m)
            while exp > 0:
                if exp & 1:
                    result = gf2_poly_mulmod(result, base, m_poly, deg_m)
                exp >>= 1
                base = gf2_poly_mulmod(base, base, m_poly, deg_m)
            return result

        def gf2_poly_gcd(a, b):
            """gcd(a, b) (GF(2) 위)."""
            while b:
                # a mod b
                while a.bit_length() >= b.bit_length():
                    a ^= b << (a.bit_length() - b.bit_length())
                a, b = b, a
            return a

        # x^(2^i) mod f for i=1..degree//2
        for i in range(1, degree // 2 + 1):
            # x^(2^i) mod poly_int
            x_pow = gf2_poly_powmod(0b10, 1 << i, poly_int, degree)
            # x^(2^i) + x (= x^(2^i) XOR x in GF(2))
            x_pow_minus_x = x_pow ^ 0b10
            g = gf2_poly_gcd(x_pow_minus_x, poly_int)
            if g != 1:
                return False
        return True


# ============================================================
# 2. Goppa 코드 생성
# ============================================================

def create_goppa_code(m, t, seed=None):
    """
    GF(2^m) 위의 이진 Goppa 코드 생성.

    N = 2^m (코드워드 길이), k ≥ N - m*t (메시지 길이)

    Args:
        m: GF(2^m) 확장 차수
        t: 에러 정정 능력 (d ≥ 2t+1)
        seed: 난수 시드

    Returns:
        G: k×N 이진 생성 행렬 (numpy array)
        goppa_poly: Goppa 다항식 계수 (GF(2^m) 원소)
        support: 서포트 L (GF(2^m) 원소 리스트, 길이 N)
    """
    rng = random.Random(seed)
    gf = GF2m(m)
    N = gf.q  # 2^m

    # 1. 랜덤 기약 Goppa 다항식 g(z) of degree t over GF(2^m)
    goppa_poly = _find_irreducible_goppa_poly(gf, t, rng)

    # 2. 서포트 L: g(α_i) ≠ 0인 α_i들만 포함
    support = [alpha for alpha in range(N) if gf.poly_eval(goppa_poly, alpha) != 0]
    N = len(support)  # 실제 코드워드 길이

    # 3. 패리티 검사 행렬 H 구성
    # H는 t×N 행렬 (GF(2^m) 원소), 이진 전개 후 mt×N
    H_bin = _build_parity_check_matrix(gf, goppa_poly, support, t)

    # 4. H_bin의 영공간으로 G 계산 (systematic form)
    G = _null_space_gf2(H_bin, N)

    return G, goppa_poly, support


def _find_irreducible_goppa_poly(gf, t, rng):
    """
    GF(2^m) 위에서 degree t의 기약다항식 찾기.

    t=1이면 모든 1차 다항식이 기약.
    t≥2이면 랜덤 생성 후 기약성 검사.
    """
    if t == 1:
        # g(z) = z + a (a ≠ 0)
        a = rng.randint(1, gf.q - 1)
        return [a, 1]  # a + z

    # t ≥ 2: 랜덤 monic 다항식 생성, 기약성 검사
    max_attempts = 10000
    for _ in range(max_attempts):
        # 랜덤 monic degree-t 다항식: coeffs[0] + ... + coeffs[t-1]*z^(t-1) + z^t
        coeffs = [rng.randint(0, gf.q - 1) for _ in range(t)] + [1]
        # 상수항이 0이면 z가 인수 → skip
        if coeffs[0] == 0:
            continue
        if _is_irreducible_over_gf2m(gf, coeffs, t):
            return coeffs

    raise RuntimeError(f"기약다항식을 {max_attempts}회 내에 찾지 못함 (m={gf.m}, t={t})")


def _is_irreducible_over_gf2m(gf, poly_coeffs, degree):
    """
    GF(2^m) 위의 다항식이 기약인지 검사.

    GF(2^m)의 모든 원소에서 evaluate → 근이 있으면 가약.
    degree=2일 때는 이것으로 충분. degree≥3이면 추가 검사 필요.
    """
    # 모든 GF(2^m) 원소에서 근 검사
    for alpha in range(gf.q):
        if gf.poly_eval(poly_coeffs, alpha) == 0:
            return False

    if degree <= 2:
        return True

    # degree ≥ 3: 2차 이상의 인수가 있는지도 검사해야 함
    # 간단한 방법: 모든 가능한 1차~(degree//2)차 인수로 나눗셈 시도
    # 여기서는 1차 인수만 없으면 OK인 경우 (degree ≤ 3)
    if degree == 3:
        # 3차에서 1차 인수가 없으면 기약
        return True

    # degree ≥ 4: 2차 인수까지 검사
    # 모든 monic 2차 다항식으로 나눗셈 시도
    for a0 in range(gf.q):
        for a1 in range(gf.q):
            divisor = [a0, a1, 1]  # a0 + a1*z + z^2
            quotient, remainder = _poly_divmod_gf2m(gf, poly_coeffs, divisor)
            if all(r == 0 for r in remainder):
                return False

    return True


def _poly_divmod_gf2m(gf, dividend, divisor):
    """GF(2^m) 위의 다항식 나눗셈: dividend = quotient * divisor + remainder."""
    dividend = list(dividend)
    deg_d = len(divisor) - 1
    deg_n = len(dividend) - 1

    if deg_n < deg_d:
        return [0], dividend

    quotient = [0] * (deg_n - deg_d + 1)
    for i in range(deg_n - deg_d, -1, -1):
        if len(dividend) - 1 >= i + deg_d:
            coeff = gf.div(dividend[i + deg_d], divisor[deg_d])
            quotient[i] = coeff
            for j in range(deg_d + 1):
                dividend[i + j] = gf.add(dividend[i + j], gf.mul(coeff, divisor[j]))

    # remainder = dividend[:deg_d]
    remainder = dividend[:deg_d] if deg_d > 0 else [0]
    return quotient, remainder


def _build_parity_check_matrix(gf, goppa_poly, support, t):
    """
    Goppa 코드의 이진 패리티 검사 행렬 구성.

    H[l, i] = α_i^l / g(α_i), l=0,...,t-1, i=0,...,N-1 (GF(2^m) 원소)
    → 이진 전개: 각 GF(2^m) 원소를 m비트로 → mt × N 이진 행렬

    Returns:
        H_bin: mt × N numpy array (GF(2))
    """
    N = len(support)
    m = gf.m

    # g(α_i)의 역원 미리 계산
    g_inv = []
    for alpha in support:
        g_alpha = gf.poly_eval(goppa_poly, alpha)
        g_inv.append(gf.inv(g_alpha))

    # H 행렬 (GF(2^m) 원소): t × N
    H_gf = np.zeros((t, N), dtype=int)
    for i, alpha in enumerate(support):
        alpha_power = 1  # α_i^0 = 1
        for l in range(t):
            H_gf[l, i] = gf.mul(alpha_power, g_inv[i])
            alpha_power = gf.mul(alpha_power, alpha)

    # 이진 전개: t × N → mt × N
    H_bin = np.zeros((m * t, N), dtype=np.uint8)
    for l in range(t):
        for i in range(N):
            val = H_gf[l, i]
            for b in range(m):
                H_bin[l * m + b, i] = (val >> b) & 1

    return H_bin


def _null_space_gf2(H, N):
    """
    GF(2) 위의 행렬 H의 영공간(kernel) 계산 → 생성 행렬 G.

    Gaussian elimination으로 H를 행 사다리꼴로 변환 후,
    자유변수에 대응하는 기저 벡터로 G 구성.

    Returns:
        G: k × N numpy array (GF(2)), k = N - rank(H)
    """
    mt, n = H.shape
    assert n == N

    # 행 사다리꼴 변환 (GF(2))
    H_rref = H.copy()
    pivot_cols = []
    row = 0
    for col in range(n):
        # 피봇 찾기
        found = False
        for r in range(row, mt):
            if H_rref[r, col] == 1:
                # 행 교환
                H_rref[[row, r]] = H_rref[[r, row]]
                found = True
                break
        if not found:
            continue

        pivot_cols.append(col)
        # 소거
        for r in range(mt):
            if r != row and H_rref[r, col] == 1:
                H_rref[r] = H_rref[r] ^ H_rref[row]
        row += 1

    rank = len(pivot_cols)
    k = N - rank

    if k <= 0:
        raise ValueError(f"Goppa 코드의 차원이 0 이하 (rank={rank}, N={N}). "
                         f"m 또는 t를 조정하세요.")

    # 자유변수 열
    free_cols = [c for c in range(N) if c not in pivot_cols]

    # 생성 행렬: 각 자유변수에 대해 기저 벡터 하나
    G = np.zeros((k, N), dtype=np.uint8)
    for idx, fc in enumerate(free_cols):
        G[idx, fc] = 1
        # 피봇 열에 대한 값 채우기
        for prow, pc in enumerate(pivot_cols):
            if prow < mt and H_rref[prow, fc] == 1:
                G[idx, pc] = 1

    return G


# ============================================================
# 3. McEliece 키 생성 & 인코딩
# ============================================================

def _random_invertible_binary_matrix(k, rng):
    """랜덤 가역 k×k 이진 행렬 생성 (GF(2) 위)."""
    max_attempts = 100
    for _ in range(max_attempts):
        S = np.array([[rng.randint(0, 1) for _ in range(k)] for _ in range(k)],
                     dtype=np.uint8)
        if _gf2_rank(S) == k:
            return S
    raise RuntimeError(f"{max_attempts}회 내에 가역 행렬을 찾지 못함 (k={k})")


def _gf2_rank(M):
    """GF(2) 위의 행렬 rank 계산."""
    mat = M.copy()
    rows, cols = mat.shape
    rank = 0
    for col in range(cols):
        found = False
        for r in range(rank, rows):
            if mat[r, col] == 1:
                mat[[rank, r]] = mat[[r, rank]]
                found = True
                break
        if not found:
            continue
        for r in range(rows):
            if r != rank and mat[r, col] == 1:
                mat[r] = mat[r] ^ mat[rank]
        rank += 1
    return rank


def _gf2_matmul(A, B):
    """GF(2) 위의 행렬 곱셈."""
    # numpy 정수 곱 후 mod 2
    return (A.astype(int) @ B.astype(int)) % 2


def mceliece_keygen(G, seed=None):
    """
    McEliece 키 생성: G_pub = S · G · P

    Args:
        G: k×N 생성 행렬
        seed: 난수 시드

    Returns:
        G_pub: k×N 공개키 행렬
        S: k×k 스크램블 행렬
        P: N×N 치환 행렬
        P_perm: 치환 인덱스 배열
    """
    rng = random.Random(seed)
    k, N = G.shape

    # S: 랜덤 가역 k×k 이진 행렬
    S = _random_invertible_binary_matrix(k, rng)

    # P: 랜덤 N×N 치환 행렬
    perm = list(range(N))
    rng.shuffle(perm)
    P = np.zeros((N, N), dtype=np.uint8)
    for i, p in enumerate(perm):
        P[i, p] = 1

    # G_pub = S · G · P
    SG = _gf2_matmul(S, G)
    G_pub = _gf2_matmul(SG, P)

    return G_pub, S, P, perm


def _to_systematic_form(G_pub):
    """
    공개키를 체계적 형태 [I_k | R]로 변환.

    열 치환을 통해 처음 k개 열이 단위 행렬이 되도록 변환.

    Returns:
        G_sys: 체계적 형태의 k×N 행렬
        col_perm: 열 치환 인덱스 (원래 인덱스)
        inv_col_perm: 역치환 인덱스
    """
    k, N = G_pub.shape
    mat = G_pub.copy()
    col_perm = list(range(N))

    # Gaussian elimination with column pivoting
    for row in range(k):
        # 피봇 찾기
        found = False
        for col in range(row, N):
            if mat[row, col] == 1:
                # 열 교환
                if col != row:
                    mat[:, [row, col]] = mat[:, [col, row]]
                    col_perm[row], col_perm[col] = col_perm[col], col_perm[row]
                found = True
                break
        if not found:
            # 아래 행에서 피봇 찾기
            for r in range(row + 1, k):
                for col in range(row, N):
                    if mat[r, col] == 1:
                        mat[[row, r]] = mat[[r, row]]
                        if col != row:
                            mat[:, [row, col]] = mat[:, [col, row]]
                            col_perm[row], col_perm[col] = col_perm[col], col_perm[row]
                        found = True
                        break
                if found:
                    break
            if not found:
                raise ValueError("G_pub를 체계적 형태로 변환 불가 (rank 부족)")

        # 소거
        for r in range(k):
            if r != row and mat[r, row] == 1:
                mat[r] = mat[r] ^ mat[row]

    # 역치환 계산
    inv_col_perm = [0] * N
    for i, c in enumerate(col_perm):
        inv_col_perm[c] = i

    return mat, col_perm, inv_col_perm


def mceliece_encrypt(target, G_pub, t, seed=None):
    """
    McEliece 인코딩: q' = target · G_pub + ε

    Args:
        target: 이진 문자열 (길이 k)
        G_pub: k×N 공개키 행렬
        t: 에러 가중치
        seed: 난수 시드

    Returns:
        ciphertext: 길이 N의 numpy 배열 (GF(2))
        error: 길이 N의 에러 벡터 (가중치 t)
    """
    rng = random.Random(seed)
    k, N = G_pub.shape
    msg = np.array([int(b) for b in target], dtype=np.uint8)

    # 인코딩: msg · G_pub
    codeword = _gf2_matmul(msg.reshape(1, k), G_pub).flatten()

    # 랜덤 에러 벡터 (가중치 t)
    error = np.zeros(N, dtype=np.uint8)
    error_positions = rng.sample(range(N), t)
    for pos in error_positions:
        error[pos] = 1

    ciphertext = (codeword + error) % 2

    return ciphertext, error


# ============================================================
# 4. p-local → QUBO 변환
# ============================================================

def _ciphertext_to_ising_terms(G_sys, ciphertext, k):
    """
    체계적 형태의 공개키와 ciphertext로부터 Ising 상호작용 항 생성.

    각 컬럼 j의 기여: H_j = [1 - (-1)^{q'_j} · Π_{i∈I_j} σ_i] / 2
    여기서 I_j는 컬럼 j에서 비영(1)인 행 인덱스, σ_i = 2x_i - 1.

    체계적 형태 [I_k | R]:
    - 처음 k개 컬럼: 가중치 1 (각 컬럼에 정확히 1개 비영) → 1-body 항
    - 나머지 N-k개 컬럼: 가중치 w (≥ 2) → w-body 항 → 차수 축소 필요

    Returns:
        terms: [(sign, var_indices), ...]
            sign: +1 또는 -1 ((-1)^{q'_j})
            var_indices: 컬럼에 관여하는 메시지 변수 인덱스 리스트
    """
    _, N = G_sys.shape
    terms = []

    for j in range(N):
        col = G_sys[:, j]
        active_rows = [i for i in range(k) if col[i] == 1]
        sign = (-1) ** int(ciphertext[j])
        terms.append((sign, active_rows))

    return terms


def _reduce_multibody_to_qubo(terms, k):
    """
    p-body Ising 항들을 반복적 Rosenberg 차수 축소로 QUBO로 변환.

    각 term: H_j = [1 - sign_j · Π_{i∈I_j} (1-2x_i)] / 2

    Π(1-2x_i)를 전개하면 x_i의 다항식. 3차 이상이면:
    - 두 변수의 곱 x_a·x_b를 보조변수 y로 대체
    - 패널티: M·(x_a·x_b - 2·x_a·y - 2·x_b·y + 3·y) 추가
    - 곱 결과를 (2y-1)로 대체하여 나머지 전개 계속

    Returns:
        Q: QUBO dict {(i,j): weight} (i <= j)
        constant: 상수항
        num_aux: 사용된 보조변수 수
    """
    Q = {}
    constant = 0.0
    next_aux = k  # 보조변수 시작 인덱스

    def add_to_Q(i, j, val):
        key = (min(i, j), max(i, j))
        Q[key] = Q.get(key, 0) + val

    for sign, var_indices in terms:
        w = len(var_indices)

        if w == 0:
            # 상수항: H_j = (1 - sign) / 2
            constant += (1 - sign) / 2.0
            continue

        if w == 1:
            # 1-body: H_j = [1 - sign·(1-2x_i)] / 2
            i = var_indices[0]
            if sign == 1:
                # [1 - (1-2x)] / 2 = x
                add_to_Q(i, i, 1.0)
            else:
                # [1 - (-(1-2x))] / 2 = [1 + 1 - 2x] / 2 = 1 - x
                constant += 1.0
                add_to_Q(i, i, -1.0)
            continue

        # w ≥ 2: 다체 상호작용 → 전개 후 2차까지 축소
        # Π_{i∈I} (1 - 2x_i)를 직접 전개
        # (1-2x_a)(1-2x_b) = 1 - 2x_a - 2x_b + 4x_a·x_b
        #
        # 반복적 축소: 한 쌍씩 처리
        # 현재 "multi-linear polynomial"을 계수로 추적
        #
        # 대신 더 직접적인 접근: 전체 곱을 전개
        #   Π_{i=1}^{w} (1-2x_i) = Σ_{S ⊆ I} (-2)^|S| Π_{i∈S} x_i
        #
        # 이를 QUBO로 변환하려면 3차 이상 항에 보조변수 도입

        # 접근: 반복적 쌍 축소
        # current_factors = list of (coeff, variables) representing the product so far
        # 초기: (1, []) (상수 1)
        # 각 변수 x_i를 곱: factor * (1 - 2*x_i) 적용

        # 다항식을 {frozenset(vars): coeff} 딕셔너리로 표현
        poly = {frozenset(): 1.0}  # 시작: 상수 1

        for xi in var_indices:
            # poly * (1 - 2*x_i)
            new_poly = {}
            for vars_set, coeff in poly.items():
                # coeff * 1
                new_poly[vars_set] = new_poly.get(vars_set, 0) + coeff
                # coeff * (-2 * x_i)
                new_vars = vars_set | frozenset([xi])
                new_poly[new_vars] = new_poly.get(new_vars, 0) - 2 * coeff
            poly = new_poly

            # 3차 이상 항이 있으면 즉시 축소
            max_degree = max(len(s) for s in poly)
            while max_degree > 2:
                poly, next_aux = _reduce_degree_once(poly, Q, next_aux)
                max_degree = max(len(s) for s in poly) if poly else 0

        # H_j = [1 - sign · poly] / 2 = 1/2 - sign/2 · poly
        constant += 0.5
        for vars_set, coeff in poly.items():
            scaled = -sign * coeff / 2.0
            if len(vars_set) == 0:
                constant += scaled
            elif len(vars_set) == 1:
                (i,) = vars_set
                add_to_Q(i, i, scaled)
            elif len(vars_set) == 2:
                i, j = sorted(vars_set)
                add_to_Q(i, j, scaled)
            else:
                # 이 시점에서 3차 이상이 남아있으면 안 됨
                raise RuntimeError(f"차수 축소 후에도 {len(vars_set)}차 항이 남음: {vars_set}")

    num_aux = next_aux - k

    # 0에 가까운 항 제거
    Q = {key: val for key, val in Q.items() if abs(val) > 1e-15}

    return Q, constant, num_aux


def _reduce_degree_once(poly, Q_penalties, next_aux):
    """
    다항식에서 가장 빈번한 변수쌍을 보조변수로 대체하여 차수를 1 줄임.

    x_a·x_b → y (보조변수), 패널티: M·(x_a·x_b - 2·x_a·y - 2·x_b·y + 3·y)

    Args:
        poly: {frozenset(var_indices): coeff}
        Q_penalties: 전체 Q dict (패널티 추가용)
        next_aux: 다음 보조변수 인덱스

    Returns:
        new_poly: 축소된 다항식
        next_aux: 업데이트된 보조변수 인덱스
    """
    # 3차 이상 항에서 가장 빈번한 변수 쌍 찾기
    pair_count = {}
    for vars_set, coeff in poly.items():
        if len(vars_set) >= 3:
            sorted_vars = sorted(vars_set)
            for i in range(len(sorted_vars)):
                for j in range(i + 1, len(sorted_vars)):
                    pair = (sorted_vars[i], sorted_vars[j])
                    pair_count[pair] = pair_count.get(pair, 0) + 1

    if not pair_count:
        return poly, next_aux

    # 가장 빈번한 쌍 선택
    best_pair = max(pair_count, key=pair_count.get)
    xa, xb = best_pair
    y = next_aux
    next_aux += 1

    # 패널티 강도 M 계산: 목적함수 항들의 최대 절댓값의 10배
    max_coeff = max(abs(c) for c in poly.values()) if poly else 1.0
    M = max(10.0 * max_coeff, 10.0)

    def add_penalty(i, j, val):
        key = (min(i, j), max(i, j))
        Q_penalties[key] = Q_penalties.get(key, 0) + val

    # 패널티: M·(x_a·x_b - 2·x_a·y - 2·x_b·y + 3·y)
    add_penalty(xa, xb, M)
    add_penalty(xa, y, -2.0 * M)
    add_penalty(xb, y, -2.0 * M)
    add_penalty(y, y, 3.0 * M)

    # 다항식에서 x_a·x_b를 y로 대체
    new_poly = {}
    for vars_set, coeff in poly.items():
        if xa in vars_set and xb in vars_set:
            # x_a·x_b → y
            new_vars = (vars_set - frozenset([xa, xb])) | frozenset([y])
            new_poly[new_vars] = new_poly.get(new_vars, 0) + coeff
        else:
            new_poly[vars_set] = new_poly.get(vars_set, 0) + coeff

    return new_poly, next_aux


# ============================================================
# 5. 메인 진입점
# ============================================================

def create_qubo_mceliece(target, m=None, t=2, seed=None):
    """
    McEliece Cryptographic QUBO 생성 (메인 진입점).

    Args:
        target: 이진 문자열 또는 길이 (정수)
        m: GF(2^m) 확장 차수 (None이면 자동 선택)
        t: 에러 가중치 (기본 2)
        seed: 난수 시드

    Returns:
        Q: QUBO dict {(i,j): weight}
        info: 메타정보 딕셔너리
    """
    # target이 정수면 길이로 해석
    if isinstance(target, int):
        rng = random.Random(seed if seed is not None else 42)
        target = ''.join(str(rng.randint(0, 1)) for _ in range(target))

    k = len(target)

    # m 자동 선택
    if m is None:
        m = _auto_select_m(k, t)

    # 시드 분리
    if seed is not None:
        goppa_seed = seed
        key_seed = seed + 100
        enc_seed = seed + 200
    else:
        goppa_seed = None
        key_seed = None
        enc_seed = None

    # 1. Goppa 코드 생성
    G, goppa_poly, support = create_goppa_code(m, t, seed=goppa_seed)
    actual_k = G.shape[0]
    N = G.shape[1]  # 실제 코드워드 길이 (support 크기)

    # target 길이 조정 (실제 k가 다를 수 있음)
    if k > actual_k:
        raise ValueError(
            f"target 길이 k={k}가 실제 코드 차원 {actual_k}를 초과 "
            f"(m={m}, t={t}, N={N}). m을 늘리거나 t를 줄이세요."
        )
    if k < actual_k:
        # 패딩: 앞에 0 추가
        target = '0' * (actual_k - k) + target
        k = actual_k

    # 2. McEliece 키 생성 → 공개키 G_pub = S·G·P
    G_pub, S_mat, P_mat, P_perm = mceliece_keygen(G, seed=key_seed)

    # 3. 체계적 형태로 변환: G_sys = [I_k | R]
    G_sys, col_perm, inv_col_perm = _to_systematic_form(G_pub)

    # 4. 체계적 형태에서 직접 인코딩
    #    c = target · G_sys + ε (체계적 프레임에서)
    #    → 첫 k비트가 target이므로 ground state = target
    rng_enc = random.Random(enc_seed)
    msg = np.array([int(b) for b in target], dtype=np.uint8)
    codeword_sys = _gf2_matmul(msg.reshape(1, k), G_sys).flatten()

    # 랜덤 에러 벡터 (가중치 t)
    error = np.zeros(N, dtype=np.uint8)
    error_positions = rng_enc.sample(range(N), t)
    for pos in error_positions:
        error[pos] = 1

    ciphertext = (codeword_sys + error) % 2

    # 5. Ising 항 생성
    terms = _ciphertext_to_ising_terms(G_sys, ciphertext, k)

    # 6. QUBO 변환 (차수 축소 포함)
    Q, constant, num_aux = _reduce_multibody_to_qubo(terms, k)

    # target에서의 에너지 계산
    # 보조변수의 최적값 결정
    aux_str = _compute_aux_values(Q, target, k, num_aux)
    full_target = target + aux_str
    target_energy = calculate_energy(full_target, Q)

    info = {
        'n': k,
        'N': N,
        'm': m,
        't': t,
        'num_aux': num_aux,
        'total_vars': k + num_aux,
        'target_energy': target_energy,
        'constant_offset': constant,
        'error_weight': int(np.sum(error)),
        'goppa_poly': goppa_poly,
        'full_target': full_target,
    }

    return Q, info


def _auto_select_m(k, t):
    """k ≤ N - m*t를 만족하는 최소 m 찾기. N은 support 크기 (≤ 2^m)."""
    for m in sorted(IRREDUCIBLE_POLYS.keys()):
        # 보수적 추정: t=1이면 N=2^m-1, t≥2면 N=2^m (기약이면 근 없음)
        N = (1 << m) - (1 if t == 1 else 0)
        expected_k = N - m * t
        if expected_k >= k:
            return m
    raise ValueError(f"k={k}, t={t}에 맞는 m을 찾지 못함 (최대 m={max(IRREDUCIBLE_POLYS.keys())})")


def _compute_aux_values(Q, target, k, num_aux):
    """
    주어진 target에 대해 보조변수의 최적값을 탐욕적으로 결정.

    각 보조변수를 순서대로, 0/1 중 에너지가 낮은 값 선택.
    """
    if num_aux == 0:
        return ''

    total_vars = k + num_aux
    current = list(target) + ['0'] * num_aux

    # 각 보조변수에 대해 최적값 결정
    for aux_idx in range(k, total_vars):
        # try 0
        current[aux_idx] = '0'
        e0 = calculate_energy(''.join(current), Q)
        # try 1
        current[aux_idx] = '1'
        e1 = calculate_energy(''.join(current), Q)

        current[aux_idx] = '0' if e0 <= e1 else '1'

    # 반복 최적화 (여러 라운드)
    for _ in range(5):
        changed = False
        for aux_idx in range(k, total_vars):
            old_val = current[aux_idx]
            current[aux_idx] = '0'
            e0 = calculate_energy(''.join(current), Q)
            current[aux_idx] = '1'
            e1 = calculate_energy(''.join(current), Q)
            best = '0' if e0 <= e1 else '1'
            if best != old_val:
                changed = True
            current[aux_idx] = best
        if not changed:
            break

    return ''.join(current[k:])


# ============================================================
# 검증 함수
# ============================================================

def extract_original_solution(sample, k):
    """SA 결과에서 원래 k개 변수만 추출."""
    return ''.join(str(sample.get(i, 0)) for i in range(k))


def verify_ground_state(Q, target, info, num_random_samples=10000):
    """
    target + 최적 aux가 ground state인지 검증.

    Returns:
        (is_local_min, is_likely_global, stats_dict)
    """
    k = info['n']
    num_aux = info['num_aux']
    total_n = k + num_aux
    full_target = info['full_target']
    target_energy = info['target_energy']

    # 1. Single-flip 이웃 검사
    min_flip_delta = float('inf')
    for i in range(total_n):
        flipped = list(full_target)
        flipped[i] = '0' if flipped[i] == '1' else '1'
        flipped_str = ''.join(flipped)
        flipped_energy = calculate_energy(flipped_str, Q)
        delta = flipped_energy - target_energy
        min_flip_delta = min(min_flip_delta, delta)

    is_local_min = min_flip_delta > -1e-10

    # 2. 랜덤 샘플 검사
    lower_count = 0
    min_random_energy = float('inf')

    for _ in range(num_random_samples):
        rand_orig = ''.join(str(random.randint(0, 1)) for _ in range(k))
        rand_aux = _compute_aux_values(Q, rand_orig, k, num_aux)
        rand_state = rand_orig + rand_aux
        rand_energy = calculate_energy(rand_state, Q)

        if rand_energy < target_energy - 1e-10:
            lower_count += 1
        min_random_energy = min(min_random_energy, rand_energy)

    is_likely_global = (lower_count == 0)

    stats = {
        'target_energy': target_energy,
        'min_flip_delta': min_flip_delta,
        'is_local_min': is_local_min,
        'lower_count': lower_count,
        'num_random_samples': num_random_samples,
        'min_random_energy': min_random_energy,
        'is_likely_global': is_likely_global,
        'full_target': full_target,
        'total_vars': total_n,
    }

    return is_local_min, is_likely_global, stats


def verify_brute_force(Q, target, info):
    """n + num_aux <= 20일 때 brute force로 ground state 완전 검증."""
    k = info['n']
    num_aux = info['num_aux']
    total_n = k + num_aux

    if total_n > 20:
        return None

    full_target = info['full_target']
    target_energy = info['target_energy']

    best_energy = float('inf')
    best_state = None
    num_degenerate = 0

    for bits in itertools.product([0, 1], repeat=total_n):
        state = ''.join(map(str, bits))
        energy = calculate_energy(state, Q)
        if energy < best_energy - 1e-10:
            best_energy = energy
            best_state = state
            num_degenerate = 1
        elif abs(energy - best_energy) < 1e-10:
            num_degenerate += 1

    return {
        'target_energy': target_energy,
        'best_energy': best_energy,
        'best_state': best_state,
        'best_original': best_state[:k] if best_state else None,
        'is_ground_state': abs(best_energy - target_energy) < 1e-10,
        'energy_gap': target_energy - best_energy,
        'num_degenerate': num_degenerate,
        'total_vars': total_n,
    }


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    target = "10110"
    m = None
    t = 2
    seed = None

    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if set(arg) <= {'0', '1'} and len(arg) >= 2:
            target = arg
        elif arg.isdigit():
            length = int(arg)
            if length < 2:
                length = 2
            random.seed(42)
            target = ''.join(str(random.randint(0, 1)) for _ in range(length))
            print(f"[설정] 길이 {length}의 랜덤 목표 해 생성")
        else:
            target = arg

    if len(sys.argv) > 2:
        m = int(sys.argv[2])

    if len(sys.argv) > 3:
        t = int(sys.argv[3])

    if len(sys.argv) > 4:
        seed = int(sys.argv[4])

    k = len(target)
    print("=" * 60)
    print("McEliece Cryptographic QUBO 생성기")
    print("=" * 60)
    print(f"Target: {target} (k={k})")
    if m is not None:
        print(f"m={m}, t={t}, N={1 << m}")
    else:
        print(f"t={t} (m은 자동 선택)")
    if seed is not None:
        print(f"Seed: {seed}")

    # QUBO 생성
    Q, info = create_qubo_mceliece(target, m=m, t=t, seed=seed)

    print(f"\n[코드 정보]")
    print(f"  GF(2^{info['m']}), N={info['N']}, k={info['n']}")
    print(f"  에러 가중치 t={info['t']}")
    print(f"  보조변수 수: {info['num_aux']}")
    print(f"  전체 QUBO 변수: {info['total_vars']}")
    print(f"  QUBO 비영 항 수: {len(Q)}")

    # 행렬 출력 (작을 때만)
    if info['total_vars'] <= 12:
        print_q_matrix(Q, info['total_vars'])
        print_qubo_formula(Q)

    # 에너지 검증
    print(f"\n[에너지 검증]")
    print(f"  Target 에너지: {info['target_energy']:.6f}")
    print(f"  상수 오프셋: {info['constant_offset']:.6f}")

    # Brute force 검증
    total_n = info['total_vars']
    if total_n <= 20:
        print(f"\n[Brute Force 검증] (전체 변수 {total_n}개)")
        bf = verify_brute_force(Q, target, info)
        print(f"  Target 에너지: {bf['target_energy']:.6f}")
        print(f"  최소 에너지:   {bf['best_energy']:.6f}")
        print(f"  최소 상태:     {bf['best_state']}")
        print(f"  원래 변수:     {bf['best_original']}")
        print(f"  축퇴도:        {bf['num_degenerate']}")
        if bf['is_ground_state']:
            print("  Ground State 검증 성공!")
        else:
            print(f"  Ground State 검증 실패 (gap: {bf['energy_gap']:.6f})")
    else:
        print(f"\n[통계적 검증] (전체 변수 {total_n}개)")
        is_local, is_global, stats = verify_ground_state(Q, target, info,
                                                          num_random_samples=5000)
        print(f"  Target 에너지: {stats['target_energy']:.6f}")
        print(f"  최소 flip delta: {stats['min_flip_delta']:.6f}")
        print(f"  Local minimum: {'OK' if is_local else 'FAIL'}")
        print(f"  랜덤 {stats['num_random_samples']}개 중 더 낮은 에너지: {stats['lower_count']}개")
        print(f"  Global minimum (추정): {'OK' if is_global else 'FAIL'}")

    # 결과 저장
    output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(output_dir, exist_ok=True)
    filename_target = f"{target[:5]}_{k}"
    output_file = os.path.join(output_dir, f"mceliece_{filename_target}_m{info['m']}_t{t}.txt")
    save_qubo_edgelist(Q, output_file, info['full_target'])
    print(f"\n[저장 완료] {output_file}")
