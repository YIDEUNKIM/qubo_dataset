import numpy as np

def check_zero_expectation(file_path):
    """
    QUBO 행렬의 비대각 요소와 대각 요소에 대해
    각 계수를 포함하는 행과 열의 평균이 0으로 수렴하는지 검사
    """
    # 파일에서 Q 행렬 읽기
    data = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines[2:]:  # 헤더 2줄 건너뛰기
            line = line.strip()
            if line:
                parts = line.split(',')
                i, j, weight = int(parts[0]), int(parts[1]), float(parts[2])
                data.append((i, j, weight))

    # 행렬 크기 결정
    max_idx = max(max(i, j) for i, j, _ in data)
    n = max_idx + 1

    # Q 행렬 생성 (대칭 행렬)
    Q = np.zeros((n, n))
    for i, j, weight in data:
        Q[i, j] = weight
        if i != j:
            Q[j, i] = weight  # 대칭

    print(f"Q 행렬 크기: {n} x {n}")
    print(f"총 요소 수: {len(data)}\n")

    # 대각 요소와 비대각 요소 분리
    diagonal_elements = []
    off_diagonal_elements = []

    for i, j, weight in data:
        if i == j:
            diagonal_elements.append((i, j, weight))
        else:
            off_diagonal_elements.append((i, j, weight))

    print(f"대각 요소 수: {len(diagonal_elements)}")
    print(f"비대각 요소 수: {len(off_diagonal_elements)}\n")

    # === 대각 요소 검사 ===
    print("=" * 60)
    print("대각 요소 검사 (Diagonal Elements)")
    print("=" * 60)

    diagonal_row_means = []
    for i, j, weight in diagonal_elements:
        # i번째 행의 모든 요소 평균
        row_mean = np.mean(Q[i, :])
        diagonal_row_means.append(row_mean)

    if diagonal_row_means:
        avg_row_mean = np.mean(diagonal_row_means)
        std_row_mean = np.std(diagonal_row_means)
        print(f"대각 요소를 포함하는 행들의 평균: {avg_row_mean:.6f}")
        print(f"표준편차: {std_row_mean:.6f}")
        print(f"최소값: {min(diagonal_row_means):.6f}")
        print(f"최대값: {max(diagonal_row_means):.6f}\n")

    # === 비대각 요소 검사 ===
    print("=" * 60)
    print("비대각 요소 검사 (Off-Diagonal Elements)")
    print("=" * 60)

    off_diagonal_row_means = []
    off_diagonal_col_means = []

    for i, j, weight in off_diagonal_elements:
        # i번째 행의 평균
        row_mean = np.mean(Q[i, :])
        off_diagonal_row_means.append(row_mean)

        # j번째 열의 평균
        col_mean = np.mean(Q[:, j])
        off_diagonal_col_means.append(col_mean)

    if off_diagonal_row_means:
        avg_row_mean = np.mean(off_diagonal_row_means)
        std_row_mean = np.std(off_diagonal_row_means)
        print(f"비대각 요소를 포함하는 행들의 평균: {avg_row_mean:.6f}")
        print(f"표준편차: {std_row_mean:.6f}")
        print(f"최소값: {min(off_diagonal_row_means):.6f}")
        print(f"최대값: {max(off_diagonal_row_means):.6f}\n")

    if off_diagonal_col_means:
        avg_col_mean = np.mean(off_diagonal_col_means)
        std_col_mean = np.std(off_diagonal_col_means)
        print(f"비대각 요소를 포함하는 열들의 평균: {avg_col_mean:.6f}")
        print(f"표준편차: {std_col_mean:.6f}")
        print(f"최소값: {min(off_diagonal_col_means):.6f}")
        print(f"최대값: {max(off_diagonal_col_means):.6f}\n")

    # === 전체 행렬 통계 ===
    print("=" * 60)
    print("전체 Q 행렬 통계")
    print("=" * 60)

    all_row_means = [np.mean(Q[i, :]) for i in range(n)]
    all_col_means = [np.mean(Q[:, j]) for j in range(n)]

    print(f"모든 행 평균의 평균: {np.mean(all_row_means):.6f}")
    print(f"모든 행 평균의 표준편차: {np.std(all_row_means):.6f}")
    print(f"모든 열 평균의 평균: {np.mean(all_col_means):.6f}")
    print(f"모든 열 평균의 표준편차: {np.std(all_col_means):.6f}\n")

    # === 0으로 수렴 판정 ===
    print("=" * 60)
    print("0으로의 수렴 판정")
    print("=" * 60)

    threshold = 1.0  # 평균이 이 값보다 작으면 0으로 수렴한다고 판정

    if abs(np.mean(all_row_means)) < threshold:
        print(f"✓ 모든 행의 평균이 0에 수렴합니다 (평균: {np.mean(all_row_means):.6f})")
    else:
        print(f"✗ 모든 행의 평균이 0에서 멀리 떨어져 있습니다 (평균: {np.mean(all_row_means):.6f})")

    if abs(np.mean(all_col_means)) < threshold:
        print(f"✓ 모든 열의 평균이 0에 수렴합니다 (평균: {np.mean(all_col_means):.6f})")
    else:
        print(f"✗ 모든 열의 평균이 0에서 멀리 떨어져 있습니다 (평균: {np.mean(all_col_means):.6f})")

    print()

if __name__ == "__main__":
    file_path = "qubo_results/qubo_00100_50.txt"
    check_zero_expectation(file_path)
