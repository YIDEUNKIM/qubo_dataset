"""
QUBO 공유 유틸리티 함수

모든 방법론에서 공통으로 사용하는 함수들.
"""

import csv


def calculate_energy(x, Q):
    """에너지 계산"""
    energy = 0
    x_vec = [int(bit) for bit in x]

    for (i, j), weight in Q.items():
        if i == j:
            energy += weight * x_vec[i]
        else:
            energy += weight * x_vec[i] * x_vec[j]

    return energy


def print_q_matrix(Q, n, max_display=12):
    """Q 행렬을 n x n 그리드 형태로 출력"""
    if n > max_display:
        print(f"\nQ 행렬이 너무 큼 ({n}x{n}). 처음 {max_display}x{max_display}만 표시:")
        display_n = max_display
    else:
        print(f"\nQ 행렬 ({n}x{n}):")
        display_n = n

    # 헤더
    print("      ", end="")
    for i in range(display_n):
        print(f"x_{i:<6}", end="")
    print("\n" + "      " + "--------" * display_n)

    # 행렬 내용
    for i in range(display_n):
        print(f"x_{i:<3} |", end="")
        for j in range(display_n):
            if i <= j:
                val = Q.get((i, j), 0)
                print(f"{val:>7.2f} ", end="")
            else:
                print(f"{'':>8}", end="")
        print()

    if n > max_display:
        print(f"      ... ({n - max_display}개 행/열 생략)")


def print_qubo_formula(Q, max_terms=20):
    """QUBO 목적 함수를 수식 형태로 출력"""
    terms = []
    sorted_keys = sorted(Q.keys())

    for (i, j) in sorted_keys:
        weight = Q[(i, j)]
        if abs(weight) < 0.001:
            continue

        sign = "+" if weight > 0 else "-"
        abs_weight = abs(weight)

        if i == j:
            term = f"{sign} {abs_weight:.2f}x_{i}"
        else:
            term = f"{sign} {abs_weight:.2f}x_{i}x_{j}"

        terms.append(term)

    if not terms:
        print("E(x) = 0")
        return

    if len(terms) > max_terms:
        formula = " ".join(terms[:max_terms]) + f" ... (+{len(terms) - max_terms}개 항)"
    else:
        formula = " ".join(terms)

    if formula.startswith("+ "):
        formula = formula[2:]

    print(f"\nQUBO 수식:\nE(x) = {formula}")


def save_qubo_edgelist(Q, filepath, target=None):
    """Edge List 형식으로 저장"""
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        if target:
            writer.writerow(['# target', target])
        writer.writerow(['i', 'j', 'weight'])

        for (i, j), weight in sorted(Q.items()):
            if abs(weight) > 1e-10:
                writer.writerow([i, j, f'{weight:.6f}'])
