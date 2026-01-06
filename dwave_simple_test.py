

import neal
import dimod
import sys
import os
import glob

def load_qubo_from_file(filepath):
    """
    qubo_zero_expectation.py에서 생성한 엣지 리스트 파일(.txt)을 읽어서
    Q 딕셔너리와 target 문자열을 반환합니다.
    """
    Q = {}
    target = None
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith("# target"):
            parts = line.split(',')
            if len(parts) >= 2:
                target = parts[1].strip()
            continue
            
        if line.startswith("i,j,weight"):
            continue
            
        # Parse output like: 0,0,-12.77
        try:
            parts = line.split(',')
            if len(parts) == 3:
                i = int(parts[0])
                j = int(parts[1])
                weight = float(parts[2])
                Q[(i, j)] = weight
        except ValueError:
            pass
            
    return Q, target

# 1. 파일에서 Q 로드 또는 최신 파일 찾기
filename = None
if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    # 가장 최근 생성된 결과 파일 찾기
    list_of_files = glob.glob('qubo_results/*.txt')
    if list_of_files:
        filename = max(list_of_files, key=os.path.getctime)
        print(f"[알림] 별도의 파일 인자가 없어 최신 파일을 사용합니다: {filename}")

if filename and os.path.exists(filename):
    print(f"Loading QUBO from: {filename}")
    Q, target_solution = load_qubo_from_file(filename)
    if not Q:
        print("Error: Q matrix is empty.")
        sys.exit(1)
else:
    print("Error: 파일을 찾을 수 없습니다. python3 dwave_simple_test.py [파일경로] 형태로 실행해주세요.")
    # Fallback for demonstration if no file exists (legacy code removed for clarity)
    sys.exit(1)

# 2. 시뮬레이티드 어닐링 샘플러 초기화 (로컬 CPU 사용)
sampler = neal.SimulatedAnnealingSampler()

# 3. 문제 풀이 실행
sampleset = sampler.sample_qubo(Q, num_reads=100)

# 4. 결과 출력
print(f"{'상태 (Binary Sequence)':<45} | {'에너지 (Energy)':<15}")
print("-" * 65)

# N 계산
n_bits = max(max(k) for k in Q.keys()) + 1 if Q else 0

for sample, energy in sampleset.data(['sample', 'energy']):
    # Compact string representation for readability
    bit_string = "".join(str(sample.get(i, 0)) for i in range(n_bits))
    print(f"{bit_string:<45} | {energy:<15.4f}")

# 최적해 요약
best_sample = sampleset.first.sample
best_energy = sampleset.first.energy
best_bitstring = "".join(str(best_sample.get(i, 0)) for i in range(n_bits))

print("\n[최적해 요약]")
print(f"최저 에너지: {best_energy:.4f}")
print(f"최적 비트 조합: {best_bitstring}")

if target_solution:
    print(f"목표 비트 조합: {target_solution}")
    if best_bitstring == target_solution:
        print("✓ SUCCESS: 목표 해를 찾았습니다!")
    else:
        print("✗ FAILURE: 목표 해와 다릅니다.")