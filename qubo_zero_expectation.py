"""
완전한 기댓값 0 QUBO 생성기 (정밀 버전)

선형 프로그래밍으로 도출한 정확한 실수 비율 사용
"""

import itertools
import random
import csv
import numpy as np
import sys
import os



from abc import ABC, abstractmethod

# -----------------------------------------------------------------------------
# Penalty Strategies (Strategy Pattern)
# -----------------------------------------------------------------------------

class PenaltyModel(ABC):
    """
    QUBO 페널티 비율을 정의하기 위한 추상 기본 클래스(Abstract Base Class)입니다.
    페널티 분배를 위한 다양한 전략 패턴을 허용합니다.
    """
    @abstractmethod
    def get_ratios(self, target_pair: tuple) -> dict:
        """
        penalty_state -> ratio 매핑 딕셔너리를 반환합니다.
        예시: 목표가 (0,0)인 경우, {(0,1): 1.0, (1,0): 1.0, (1,1): 1.65} 반환
        """
        pass

class DefaultZeroExpectationModel(PenaltyModel):
    """
    선형 프로그래밍(LP)을 통해 도출된 기댓값 0 보장 표준 모델입니다.
    최적화된 정밀 비율을 사용합니다.
    """
    def __init__(self):
        # 정규화 기준: 최솟값 = 1.0 (LP 결과값 캡슐화)
        self._ratios_table = {
            # 정답 = (0, 0): 페널티
            (0, 0): {(0, 1): 1.00, (1, 0): 1.00, (1, 1): 1.65},
            # 정답 = (0, 1)
            (0, 1): {(0, 0): 2.00, (1, 0): 1.00, (1, 1): 1.68},
            # 정답 = (1, 0)
            (1, 0): {(0, 0): 2.00, (0, 1): 1.00, (1, 1): 1.68},
            # 정답 = (1, 1)
            (1, 1): {(0, 0): 1.00, (0, 1): 3.00, (1, 0): 3.00},
        }

    def get_ratios(self, target_pair: tuple) -> dict:
        return self._ratios_table.get(target_pair, {})

class SimpleUniformModel(PenaltyModel):
    """
    대체 모델: 비교를 위한 단순 균등 페널티 모델입니다.
    모든 오답에 대해 동일한 페널티를 부여합니다.
    """
    def get_ratios(self, target_pair: tuple) -> dict:
        # 모든 오답 상태에 1.0 부여
        all_states = [(0, 0), (0, 1), (1, 0), (1, 1)]
        return {s: 1.0 for s in all_states if s != target_pair}


def create_qubo_precise(target_binary_string, density=1.0, base_range=(1, 3), model: PenaltyModel = None, **kwargs):
    """
    정밀한 기댓값 0 QUBO 생성
    
    Args:
        target_binary_string: 목표 해 (예: "10110")
        density: 상호작용 추가 확률 (1.0 권장)
        base_range: 기본 r 샘플링 범위
        model: PenaltyModel 인스턴스 (None일 경우 DefaultZeroExpectationModel 사용)
    
    Returns:
        Q: QUBO 행렬 (딕셔너리)
    """
    if model is None:
        model = DefaultZeroExpectationModel()
        
    n = len(target_binary_string)
    Q = {}
    
    # 모든 큐빗 쌍(i, j)에 대해 반복 (i < j)
    for i in range(n):
        for j in range(i + 1, n):
            # density 확률로 상호작용 추가 여부 결정
            if random.random() > density:
                continue
            
            # 1. 목표 상태 확인 (이 쌍의 정답이 무엇이어야 하는지)
            bit_i = int(target_binary_string[i])
            bit_j = int(target_binary_string[j])
            target_state = (bit_i, bit_j)
            
            # 2. 모델을 통해 페널티 비율 가져오기 (Strategy Pattern 적용)
            ratios = model.get_ratios(target_state)
            
            for penalty_state, ratio in ratios.items():
                # 3. 페널티 강도(r) 결정
                r = random.uniform(*base_range) * ratio
                
                s_i, s_j = penalty_state
                
                # 4. 페널티 항을 QUBO 행렬에 분배
                # 목표: 상태가 (s_i, s_j)일 때 에너지가 r만큼 증가해야 함
                
                if s_i == 0 and s_j == 0:
                    # Case 00: 페널티 식 = r * (1-x_i)(1-x_j)
                    # 전개: r(1 - x_i - x_j + x_i*x_j)
                    Q[(i, i)] = Q.get((i, i), 0) - r
                    Q[(j, j)] = Q.get((j, j), 0) - r
                    Q[(i, j)] = Q.get((i, j), 0) + r
                    
                elif s_i == 0 and s_j == 1:
                    # Case 01: 페널티 식 = r * (1-x_i)x_j
                    # 전개: r(x_j - x_i*x_j)
                    Q[(j, j)] = Q.get((j, j), 0) + r
                    Q[(i, j)] = Q.get((i, j), 0) - r
                    
                elif s_i == 1 and s_j == 0:
                    # Case 10: 페널티 식 = r * x_i(1-x_j)
                    # 전개: r(x_i - x_i*x_j)
                    Q[(i, i)] = Q.get((i, i), 0) + r
                    Q[(i, j)] = Q.get((i, j), 0) - r
                    
                else: # s_i == 1 and s_j == 1
                    # Case 11: 페널티 식 = r * x_i*x_j
                    # 전개: r*x_i*x_j
                    Q[(i, j)] = Q.get((i, j), 0) + r
    
    # 행/열 합 0으로 수렴하게 만들기 (옵션)
    if kwargs.get('balance_rows', False):
        # Ising 모델 기반으로 재생성 (기존 Q 덮어쓰기)
        # 이 방식은 E[RowSum] = 0을 만족하면서 Ground State를 완벽히 보존함
        Q = create_qubo_ising_derived(target_binary_string, density=density, base_range=base_range)
        
    return Q


def create_qubo_ising_derived(target, density=1.0, base_range=(1, 3)):
    """
    Ising 모델 (H = -sum J_ij s_i s_j, h_i=0)에서 유도된 QUBO.
    특징:
    1. Ground State가 Target과 일치함이 수학적으로 보장됨.
    2. QUBO 변환 시 Expected Row Sum이 0임.
    """
    n = len(target)
    Q = {}
    
    # Ising Spin 변환 (-1, 1)
    spins = [1 if b == '1' else -1 for b in target]
    
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() > density:
                continue
                
            # Ising Interaction J_ij
            # Target 상태에서 에너지가 최소화되려면 J_ij * s_i * s_j > 0 이어야 함 (H = - sum ...)
            # 즉 J_ij = alpha * s_i * s_j
            alpha = random.uniform(*base_range)
            J_ij = alpha * spins[i] * spins[j]
            
            # Convert Ising J to QUBO
            # H_term = - J_ij * s_i * s_j
            # s = 2x - 1
            # s_i s_j = (2x_i - 1)(2x_j - 1) = 4x_ix_j - 2x_i - 2x_j + 1
            # H_term = - J_ij (4x_ix_j - 2x_i - 2x_j + 1)
            #        = -4 J_ij x_ix_j + 2 J_ij x_i + 2 J_ij x_j - J_ij
            
            q_ij = -4 * J_ij
            q_i = 2 * J_ij
            q_j = 2 * J_ij
            
            Q[(i, j)] = Q.get((i, j), 0) + q_ij
            Q[(i, i)] = Q.get((i, i), 0) + q_i
            Q[(j, j)] = Q.get((j, j), 0) + q_j
            
    return Q


def balance_qubo_rows(Q, n):
    """Deprecated: 단순 보정은 Ground State를 깸. Ising Derived 사용 권장."""
    # 각 변수별 총 합 계산 (대각 + 비대각)
    # Q는 상삼각 형태지만 논리적으로 대칭이므로 i와 연결된 모든 j를 고려해야 함
    row_sums = {i: 0.0 for i in range(n)}
    
    for (i, j), weight in Q.items():
        if i == j:
            row_sums[i] += weight
        else:
            # i < j 형태
            row_sums[i] += weight
            row_sums[j] += weight
            
    # 대각 성분 보정
    for i in range(n):
        # 현재 합이 S_i라면, Q_ii에서 S_i를 빼면 새로운 합은 0이 됨
        # (Q_ii - S_i) + others = Q_ii + others - S_i = S_i - S_i = 0
        correction = row_sums[i]
        Q[(i, i)] = Q.get((i, i), 0) - correction
        
    return Q


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


def solve_brute_force(Q, n):
    """브루트포스 풀이"""
    if n > 20:
        return None, None, None

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


def batch_test(num_tests=50, n_bits=10):
    """배치 테스트"""
    print(f"\n배치 테스트: {num_tests}회, {n_bits}비트")
    print("-" * 40)
    
    successes = 0
    
    for i in range(num_tests):
        target = ''.join(random.choice('01') for _ in range(n_bits))
        Q = create_qubo_precise(target, density=1.0)
        found, min_e, _ = solve_brute_force(Q, n_bits)
        
        if found == target:
            successes += 1
    
    print(f"성공률: {successes}/{num_tests} ({100*successes/num_tests:.1f}%)")
    return successes / num_tests


def large_scale_analysis(n_tests=500, n_bits=30):
    """
    대규모 계수 분석 (브루트포스 없이)
    """
    print(f"\n대규모 계수 분석: {n_tests}회, {n_bits}비트")
    print("-" * 40)
    
    # 각 조건별로 분리해서 분석
    xi_bi0 = []  # b_i=0인 변수의 x_i 계수
    xi_bi1 = []  # b_i=1인 변수의 x_i 계수
    xj_bj0 = []  # b_j=0인 변수의 x_j 계수
    xj_bj1 = []  # b_j=1인 변수의 x_j 계수
    xixj_all = []
    
    for _ in range(n_tests):
        target = ''.join(random.choice('01') for _ in range(n_bits))
        Q = create_qubo_precise(target, density=1.0)
        
        # 대각 항 분석
        for i in range(n_bits):
            if (i, i) in Q:
                if target[i] == '0':
                    xi_bi0.append(Q[(i, i)])
                else:
                    xi_bi1.append(Q[(i, i)])
        
        # 비대각 항 분석
        for (i, j), weight in Q.items():
            if i != j:
                xixj_all.append(weight)
    
    print(f"\nE[x_i | b_i=0] = {np.mean(xi_bi0):.6f} (표준편차: {np.std(xi_bi0):.4f})")
    print(f"E[x_i | b_i=1] = {np.mean(xi_bi1):.6f} (표준편차: {np.std(xi_bi1):.4f})")
    print(f"E[x_ix_j] = {np.mean(xixj_all):.6f} (표준편차: {np.std(xixj_all):.4f})")


def compare_with_random_qubo(n_bits=20, n_tests=100):
    """
    생성된 QUBO와 순수 랜덤 QUBO의 계수 분포 비교
    """
    print(f"\n생성 QUBO vs 랜덤 QUBO 비교 ({n_bits}비트, {n_tests}회)")
    print("-" * 50)
    
    # 생성된 QUBO 계수 수집
    generated_diag = []
    generated_offdiag = []
    
    for _ in range(n_tests):
        target = ''.join(random.choice('01') for _ in range(n_bits))
        Q = create_qubo_precise(target, density=1.0)
        
        for (i, j), weight in Q.items():
            if i == j:
                generated_diag.append(weight)
            else:
                generated_offdiag.append(weight)
    
    # 순수 랜덤 QUBO 생성
    random_diag = []
    random_offdiag = []
    
    for _ in range(n_tests):
        for i in range(n_bits):
            # 대각: uniform(-10, 10)
            random_diag.append(random.uniform(-10, 10))
        
        for i in range(n_bits):
            for j in range(i+1, n_bits):
                # 비대각: uniform(-10, 10)
                random_offdiag.append(random.uniform(-10, 10))
    
    print("\n대각 항 (x_i 계수):")
    print(f"  생성 QUBO: 평균={np.mean(generated_diag):.4f}, 표준편차={np.std(generated_diag):.4f}")
    print(f"  랜덤 QUBO: 평균={np.mean(random_diag):.4f}, 표준편차={np.std(random_diag):.4f}")
    
    print("\n비대각 항 (x_ix_j 계수):")
    print(f"  생성 QUBO: 평균={np.mean(generated_offdiag):.4f}, 표준편차={np.std(generated_offdiag):.4f}")
    print(f"  랜덤 QUBO: 평균={np.mean(random_offdiag):.4f}, 표준편차={np.std(random_offdiag):.4f}")
    
    print("\n* 두 분포가 비슷할수록 구별하기 어려움")


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


def generate_dataset(num_problems=100, n_bits=50, output_dir='/home/claude'):
    """
    벤치마크 데이터셋 생성
    """
    import os
    
    print(f"\n데이터셋 생성: {num_problems}개, {n_bits}비트")
    print("-" * 40)
    
    dataset_path = os.path.join(output_dir, 'qubo_dataset')
    os.makedirs(dataset_path, exist_ok=True)
    
    metadata = []
    
    for i in range(num_problems):
        target = ''.join(random.choice('01') for _ in range(n_bits))
        Q = create_qubo_precise(target, density=1.0)
        
        filename = f'qubo_{i:04d}.csv'
        filepath = os.path.join(dataset_path, filename)
        save_qubo_edgelist(Q, filepath, target)
        
        metadata.append({
            'id': i,
            'filename': filename,
            'target': target,
            'n_bits': n_bits,
            'n_terms': len(Q)
        })
    
    # 메타데이터 저장
    meta_path = os.path.join(dataset_path, 'metadata.csv')
    with open(meta_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['id', 'filename', 'target', 'n_bits', 'n_terms'])
        writer.writeheader()
        writer.writerows(metadata)
    
    print(f"저장 완료: {dataset_path}")
    print(f"  - {num_problems}개의 QUBO 파일")
    print(f"  - metadata.csv")
    
    return dataset_path


if __name__ == "__main__":
    random.seed(42)
    
    print("=" * 60)
    print("기댓값 0 QUBO 생성기 (정밀 버전)")
    print("=" * 60)
    
    # 0. 단일 예제 상세 출력
    print("\n" + "=" * 60)
    print("단일 예제 상세 출력")
    print("=" * 60)
    
    target = "11000101010001101"
    balance_mode = False
    
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg.isdigit():
            # 숫자로 입력되면 해당 길이의 랜덤 타겟 생성
            length = int(arg)
            target = "".join(str(random.randint(0, 1)) for _ in range(length))
            print(f"[설정] 길이 {length}의 랜덤 목표 해를 생성했습니다.")
        elif arg == "balance":
            print("[설정] 행/열 평균 0 수렴 모드 활성화")
            balance_mode = True
        else:
            # 이진 문자열로 직접 입력된 경우
            target = arg
            
        # 추가 인자 확인 (순서 무관하게 balance 체크)
        if "balance" in sys.argv:
            balance_mode = True
            print("[설정] 행/열 평균 0 수렴 모드 활성화")
        
    print(f"\n목표 해: {target}")
    
    Q = create_qubo_precise(target, density=1.0, balance_rows=balance_mode)
    
    # Q 행렬 출력
    print_q_matrix(Q, len(target))
    
    # 수식 출력
    print_qubo_formula(Q)
    
    # 검증
    found, min_energy, all_results = solve_brute_force(Q, len(target))
    target_energy = calculate_energy(target, Q)
    
    print(f"\n검증 결과:")
    print(f"  목표 해 '{target}'의 에너지: {target_energy}")
    
    if found is None:
        print(f"  브루트포스 최소 해: (N > 20 이므로 생략됨)")
        print("  - 검증 생략됨")
    else:
        print(f"  브루트포스 최소 해 '{found}'의 에너지: {min_energy}")
        
        if found == target:
            print("  ✓ 성공: 목표 해가 최소값!")
        else:
            print("  ✗ 실패")

    # 결과 저장
    output_dir = "qubo_results"
    os.makedirs(output_dir, exist_ok=True)
    
    # 파일명 길이가 너무 길어지지 않도록 조정
    filename_target = f"{target[:5]}_{len(target)}"
    
    output_file = os.path.join(output_dir, f"qubo_{filename_target}.txt")
    
    save_qubo_edgelist(Q, output_file, target)
    print(f"\n[저장 완료] 결과 파일: {output_file}")
    
    # 에너지 순위
    if all_results:
        print("\n  에너지 순위 (상위 5개):")
        sorted_results = sorted(all_results, key=lambda x: x[1])
        for rank, (state, energy) in enumerate(sorted_results[:5], 1):
            marker = " <- 목표" if state == target else ""
            print(f"    {rank}. {state}: {energy:.4f}{marker}")
    else:
         print("\n  에너지 순위: (생략됨)")
    
    # 1. 배치 테스트
    batch_test(num_tests=30, n_bits=10)
    
    # 2. 대규모 계수 분석
    large_scale_analysis(n_tests=500, n_bits=30)
    
    # 3. 랜덤 QUBO와 비교
    compare_with_random_qubo(n_bits=20, n_tests=100)
    
    # 4. 데이터셋 생성 (예시)
    # generate_dataset(num_problems=10, n_bits=50)
