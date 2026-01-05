import itertools
import random


def create_dense_qubo_for_target(target_binary_string, density=0.5, balance_expectation=True):
    """
    PDF에 명시된 방법대로 QUBO를 생성합니다.
    
    Step 1: 대각 항 - b_i=0이면 +a_i*x_i, b_i=1이면 -a_i*x_i
    Step 2: 비대각 항 - optimal이 아닌 조합에 각각 독립적인 페널티 부여
    
    Args:
        target_binary_string: 목표 해 (예: "10110")
        density: 상호작용 추가 확률
        balance_expectation: True면 기댓값 0을 위한 비율 조정 적용
    """
    n = len(target_binary_string)
    Q = {}
    
    print(f"목표 해: {target_binary_string}")
    print(f"기댓값 균형 조정: {'ON' if balance_expectation else 'OFF'}")
    print("-" * 50)
    
    # Step 1: 대각 항 생성
    # b_i = 0이면 a_i * x_i, b_i = 1이면 -a_i * x_i
    for i in range(n):
        bit = int(target_binary_string[i])
        a_i = random.uniform(1, 3)  # 임의의 양수
        if bit == 0:
            Q[(i, i)] = Q.get((i, i), 0) + a_i
        else:
            Q[(i, i)] = Q.get((i, i), 0) - a_i
    
    # Step 2: 비대각 항 생성
    # optimal이 아닌 3가지 조합에 각각 독립적인 페널티 부여
    interaction_count = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() > density:
                continue
            
            bit_i = int(target_binary_string[i])
            bit_j = int(target_binary_string[j])
            
            # 기댓값 0을 위한 샘플링 함수
            # PDF: E[r0], E[r0']를 다른 값의 3배로 샘플링
            def sample_r(is_special=False):
                if balance_expectation and is_special:
                    return random.uniform(3, 9)  # 3배 범위
                else:
                    return random.uniform(1, 3)
            
            if bit_i == 0 and bit_j == 0:
                # optimal = 00, 페널티 대상: 01, 10, 11
                # PDF 기준: 01이 r0 (special), 10이 r1, 11이 r2
                r0 = sample_r(is_special=True)   # 01 페널티 (special)
                r1 = sample_r(is_special=False)  # 10 페널티
                r2 = sample_r(is_special=False)  # 11 페널티
                
                # 01: r0 * (1-x_i) * x_j = r0 * (x_j - x_i*x_j)
                Q[(j, j)] = Q.get((j, j), 0) + r0
                Q[(i, j)] = Q.get((i, j), 0) - r0
                
                # 10: r1 * x_i * (1-x_j) = r1 * (x_i - x_i*x_j)
                Q[(i, i)] = Q.get((i, i), 0) + r1
                Q[(i, j)] = Q.get((i, j), 0) - r1
                
                # 11: r2 * x_i * x_j
                Q[(i, j)] = Q.get((i, j), 0) + r2
                
                print(f"({i},{j}) b=00: 01페널티={r0:.2f}, 10페널티={r1:.2f}, 11페널티={r2:.2f}")
                
            elif bit_i == 0 and bit_j == 1:
                # optimal = 01, 페널티 대상: 00, 10, 11
                # PDF 기준: 00이 r0' (special)
                r0 = sample_r(is_special=True)   # 00 페널티 (special)
                r1 = sample_r(is_special=False)  # 10 페널티
                r2 = sample_r(is_special=False)  # 11 페널티
                
                # 00: r0 * (1-x_i) * (1-x_j) = r0 * (1 - x_i - x_j + x_i*x_j)
                # 상수항 무시
                Q[(i, i)] = Q.get((i, i), 0) - r0
                Q[(j, j)] = Q.get((j, j), 0) - r0
                Q[(i, j)] = Q.get((i, j), 0) + r0
                
                # 10: r1 * x_i * (1-x_j) = r1 * (x_i - x_i*x_j)
                Q[(i, i)] = Q.get((i, i), 0) + r1
                Q[(i, j)] = Q.get((i, j), 0) - r1
                
                # 11: r2 * x_i * x_j
                Q[(i, j)] = Q.get((i, j), 0) + r2
                
                print(f"({i},{j}) b=01: 00페널티={r0:.2f}, 10페널티={r1:.2f}, 11페널티={r2:.2f}")
                
            elif bit_i == 1 and bit_j == 0:
                # optimal = 10, 페널티 대상: 00, 01, 11
                r0 = sample_r(is_special=True)   # 00 페널티 (special)
                r1 = sample_r(is_special=False)  # 01 페널티
                r2 = sample_r(is_special=False)  # 11 페널티
                
                # 00: r0 * (1-x_i) * (1-x_j)
                Q[(i, i)] = Q.get((i, i), 0) - r0
                Q[(j, j)] = Q.get((j, j), 0) - r0
                Q[(i, j)] = Q.get((i, j), 0) + r0
                
                # 01: r1 * (1-x_i) * x_j = r1 * (x_j - x_i*x_j)
                Q[(j, j)] = Q.get((j, j), 0) + r1
                Q[(i, j)] = Q.get((i, j), 0) - r1
                
                # 11: r2 * x_i * x_j
                Q[(i, j)] = Q.get((i, j), 0) + r2
                
                print(f"({i},{j}) b=10: 00페널티={r0:.2f}, 01페널티={r1:.2f}, 11페널티={r2:.2f}")
                
            else:  # bit_i == 1 and bit_j == 1
                # optimal = 11, 페널티 대상: 00, 01, 10
                r0 = sample_r(is_special=True)   # 00 페널티 (special)
                r1 = sample_r(is_special=False)  # 01 페널티
                r2 = sample_r(is_special=False)  # 10 페널티
                
                # 00: r0 * (1-x_i) * (1-x_j)
                Q[(i, i)] = Q.get((i, i), 0) - r0
                Q[(j, j)] = Q.get((j, j), 0) - r0
                Q[(i, j)] = Q.get((i, j), 0) + r0
                
                # 01: r1 * (1-x_i) * x_j = r1 * (x_j - x_i*x_j)
                Q[(j, j)] = Q.get((j, j), 0) + r1
                Q[(i, j)] = Q.get((i, j), 0) - r1
                
                # 10: r2 * x_i * (1-x_j) = r2 * (x_i - x_i*x_j)
                Q[(i, i)] = Q.get((i, i), 0) + r2
                Q[(i, j)] = Q.get((i, j), 0) - r2
                
                print(f"({i},{j}) b=11: 00페널티={r0:.2f}, 01페널티={r1:.2f}, 10페널티={r2:.2f}")
            
            interaction_count += 1
    
    print(f"\n-> 총 {interaction_count}개의 상호작용 쌍 추가")
    return Q


def calculate_energy(x, Q):
    """주어진 이진 문자열 x에 대해 에너지 E = x^T Q x를 계산"""
    energy = 0
    x_vec = [int(bit) for bit in x]
    
    for (i, j), weight in Q.items():
        if i == j:
            energy += weight * x_vec[i]
        else:
            energy += weight * x_vec[i] * x_vec[j]
    
    return energy


def solve_brute_force(Q, n):
    """모든 2^n가지 가능성을 확인하여 최소 에너지 상태를 찾음"""
    if n > 20:
        print(f"\n[경고] n={n}은 브루트 포스로 검증하기에 너무 큽니다.")
        return None, None

    best_energy = float('inf')
    best_solution = None
    all_results = []
    
    for bits in itertools.product([0, 1], repeat=n):
        current_state = "".join(map(str, bits))
        energy = calculate_energy(current_state, Q)
        all_results.append((current_state, energy))
        
        if energy < best_energy:
            best_energy = energy
            best_solution = current_state
    
    return best_solution, best_energy, all_results


def print_q_matrix(Q, n):
    """Q 행렬을 n x n 그리드 형태로 출력"""
    print(f"\nQ 행렬 ({n}x{n}):")
    
    print("      ", end="")
    for i in range(n):
        print(f"x_{i:<6}", end="")
    print("\n" + "      " + "--------" * n)
    
    for i in range(n):
        print(f"x_{i:<3} |", end="")
        for j in range(n):
            if i <= j:
                val = Q.get((i, j), 0)
                print(f"{val:>7.2f} ", end="")
            else:
                print(f"{'':>8}", end="")
        print()


def print_qubo_formula(Q):
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

    formula = " ".join(terms)
    if formula.startswith("+ "):
        formula = formula[2:]
    
    print(f"\nQUBO 수식:\nE(x) = {formula}")


def verify_correctness(target, Q, n):
    """목표 해가 실제로 최소 에너지를 갖는지 검증"""
    print("\n" + "=" * 50)
    print("검증 결과")
    print("=" * 50)
    
    found_solution, min_energy, all_results = solve_brute_force(Q, n)
    
    if found_solution is None:
        return
    
    # 목표 해의 에너지
    target_energy = calculate_energy(target, Q)
    
    print(f"목표 해 '{target}'의 에너지: {target_energy:.4f}")
    print(f"브루트포스 최소 해 '{found_solution}'의 에너지: {min_energy:.4f}")
    
    if found_solution == target:
        print("\n✓ 성공: 목표 해가 유일한 최소값입니다!")
    elif abs(target_energy - min_energy) < 0.0001:
        print("\n△ 부분 성공: 목표 해가 최소값이지만 다른 해도 같은 에너지를 가짐")
        # 동일 에너지를 가진 다른 해 찾기
        same_energy = [s for s, e in all_results if abs(e - min_energy) < 0.0001 and s != target]
        if same_energy:
            print(f"   동일 에너지 해: {same_energy[:5]}")
    else:
        print(f"\n✗ 실패: 목표 해가 최소값이 아닙니다!")
        print(f"   에너지 차이: {target_energy - min_energy:.4f}")
    
    # 상위 5개 해 출력
    print("\n에너지 순위 (상위 5개):")
    sorted_results = sorted(all_results, key=lambda x: x[1])
    for rank, (state, energy) in enumerate(sorted_results[:5], 1):
        marker = " <- 목표" if state == target else ""
        print(f"  {rank}. {state}: {energy:.4f}{marker}")


if __name__ == "__main__":
    random.seed(42)  # 재현성을 위한 시드 고정
    
    target = "10110"
    
    print("=" * 50)
    print("PDF 방식 QUBO 생성 테스트")
    print("=" * 50)
    
    # Q 행렬 생성
    Q = create_dense_qubo_for_target(target, density=0.7, balance_expectation=True)
    
    # Q 행렬 출력
    print_q_matrix(Q, len(target))
    
    # 수식 출력
    print_qubo_formula(Q)
    
    # 검증
    verify_correctness(target, Q, len(target))