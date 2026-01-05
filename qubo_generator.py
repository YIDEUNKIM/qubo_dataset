import itertools

def create_dense_qubo_for_target(target_binary_string, density=0.5):
    """
    목표 해를 유지하면서 변수 간 상호작용(비대각 성분)이 있는 복잡한 QUBO를 생성합니다.
    
    원리:
    기본 식 sum(x_i - b_i)^2 에 다음 항들을 추가해도 정답(최소값 0)은 변하지 않습니다.
    1. b_i == b_j 인 경우: (x_i - x_j)^2 추가 (둘이 같으면 0, 다르면 벌점)
    2. b_i != b_j 인 경우: (x_i + x_j - 1)^2 추가 (둘이 다르면 0, 같으면 벌점)
    """
    n = len(target_binary_string)
    Q = {}
    
    print(f"목표 해: {target_binary_string} (상호작용 포함)")
    print("-" * 30)
    
    # 1. 기본 대각 항 생성 (Base)
    # E = sum(x_i - b_i)^2
    # (x_i - b_i)^2 = x_i^2 - 2x_i b_i + b_i^2 = x_i(1 - 2b_i) + const
    for i in range(n):
        bit = int(target_binary_string[i])
        Q[(i, i)] = Q.get((i, i), 0) + (1 - 2 * bit)

    # 2. 상호작용 항 추가 (Interactions)
    import random
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            # density 확률로 상호작용 추가
            if random.random() > density:
                continue
            
            bit_i = int(target_binary_string[i])
            bit_j = int(target_binary_string[j])
            
            # 계수 (강도)를 랜덤하게 설정하여 문제의 다양성 확보
            strength = random.randint(1, 3)
            
            if bit_i == bit_j:
                # 둘의 값이 같아야 함: strength * (x_i - x_j)^2
                # = strength * (x_i^2 - 2x_i x_j + x_j^2)
                # = strength * (x_i - 2x_i x_j + x_j)  (x^2=x)
                Q[(i, i)] += strength
                Q[(j, j)] += strength
                Q[(i, j)] = Q.get((i, j), 0) - 2 * strength
                
                print(f"상호작용 추가 ({i}, {j}): 같은 값 선호 (강도 {strength})")
            else:
                # 둘의 값이 달라야 함: strength * (x_i + x_j - 1)^2
                # = strength * (x_i^2 + x_j^2 + 1 + 2x_i x_j - 2x_i - 2x_j)
                # = strength * (2x_i x_j - x_i - x_j + 1)
                Q[(i, i)] -= strength
                Q[(j, j)] -= strength
                Q[(i, j)] = Q.get((i, j), 0) + 2 * strength
                
                print(f"상호작용 추가 ({i}, {j}): 다른 값 선호 (강도 {strength})")
            
            count += 1
            
    print(f"-> 총 {count}개의 상호작용 항이 추가되었습니다.")
    return Q

def calculate_energy(x, Q):
    """
    주어진 이진 문자열 x에 대해 에너지 E = x^T Q x를 계산합니다.
    """
    energy = 0
    # x는 "101"과 같은 문자열이므로 정수 리스트로 변환
    x_vec = [int(bit) for bit in x]
    
    for (i, j), weight in Q.items():
        if i == j:
            energy += weight * x_vec[i]
        else:
            energy += weight * x_vec[i] * x_vec[j]
            
    return energy

    return best_solution, best_energy

def solve_brute_force(Q, n):
    """
    모든 2^n가지 가능성을 확인하여 최소 에너지 상태를 찾습니다.
    주의: n이 20을 넘어가면 시간이 매우 오래 걸릴 수 있습니다.
    """
    if n > 20:
        print(f"\n[경고] n={n} (2^{n} = {2**n:,} 상태)는 브루트 포스로 검증하기에 너무 큽니다.")
        print("검증을 건너뜁니다.")
        return None, None

    best_energy = float('inf')
    best_solution = None
    
    print(f"\n모든 가능한 상태 검증 (총 {2**n}개):")
    # n이 5보다 크면 모든 출력을 보여주지 않고 요약만 함
    show_all = n <= 5
    
    if show_all:
        print(f"{'상태':<10} | {'에너지':<10}")
        print("-" * 25)
    
    # 모든 2^n 조합에 대해 반복
    for i, bits in enumerate(itertools.product([0, 1], repeat=n)):
        current_state = "".join(map(str, bits))
        energy = calculate_energy(current_state, Q)
        
        if show_all:
            print(f"{current_state:<10} | {energy:<10}")
        
        if energy < best_energy:
            best_energy = energy
            best_solution = current_state
            
    return best_solution, best_energy

def print_qubo_formula(Q):
    """
    QUBO 목적 함수를 수식 형태로 출력합니다.
    E(x) = sum(Q_ii * x_i) + sum(Q_ij * x_i * x_j)
    """
    terms = []
    
    # 키를 정렬하여 순서대로 출력 (i, j)
    sorted_keys = sorted(Q.keys())
    
    for (i, j) in sorted_keys:
        weight = Q[(i, j)]
        if weight == 0:
            continue
            
        sign = "+" if weight > 0 else "-"
        abs_weight = abs(weight)
        
        # 계수가 1이면 생략 (단, 상수항이 아니므로 변수명은 항상 붙음)
        weight_str = f"{abs_weight}" if abs_weight != 1 else ""
        
        if i == j:
            # 이진 변수이므로 x_i^2 = x_i 이지만, QUBO(2차)의 의미를 살리기 위해 제곱으로 표기
            term = f"{sign} {weight_str}x_{i}^2"
        else:
            term = f"{sign} {weight_str}x_{i}x_{j}"
            
        terms.append(term)
        
    if not terms:
        print("E(x) = 0")
        return

    # 첫 번째 항의 '+' 부호 처리를 위해 조작
    formula = " ".join(terms)
    if formula.startswith("+ "):
        formula = formula[2:] # 맨 앞 "+ " 제거
    
    print(f"\nQUBO 수식:\nE(x) = {formula}")

def print_q_matrix(Q, n):
    """
    Q 행렬을 n x n 그리드 형태로 출력합니다.
    """
    print(f"\nQ 행렬 ({n}x{n}):")
    
    # 열 헤더
    print("      ", end="")
    for i in range(n):
        print(f"x_{i:<3}", end="")
    print("\n" + "      " + "-----" * n)
    
    for i in range(n):
        print(f"x_{i:<3} |", end="")
        for j in range(n):
            if i <= j:
                val = Q.get((i, j), 0)
                print(f"{val:^5}", end="")
            else:
                print(f"{' ':^5}", end="") # 하삼각은 공백 처리
        print()

if __name__ == "__main__":
    # 예제: 목표 "10110"에 대한 QUBO 생성
    target = "10110"
    
    # 1. Q 행렬 생성 (상호작용 포함 버전 사용)
    # Q = create_qubo_for_target(target) # 기존 단순 버전
    Q = create_dense_qubo_for_target(target, density=0.5)
    
    # Q 행렬 출력
    print_q_matrix(Q, len(target))
    
    # 수식 출력
    print_qubo_formula(Q)
    
    # 2. 브루트 포스로 검증 (n이 작을 때만 수행)
    found_solution, min_energy = solve_brute_force(Q, len(target))
    
    if found_solution is not None:
        print("-" * 25)
        print(f"찾은 최소 상태: {found_solution}")
        print(f"최소 에너지: {min_energy}")
        
        if found_solution == target:
            print("\n성공: QUBO 문제가 올바르게 구성되었습니다!")
        else:
            print("\n실패: 찾은 해가 목표와 일치하지 않습니다.")
