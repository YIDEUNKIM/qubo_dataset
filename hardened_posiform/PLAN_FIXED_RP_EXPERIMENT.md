# 구현 계획: R/P 고정 α 비교 실험

## 목적

기존 실험에서는 `create_qubo_hardened_posiform()`을 α마다 호출하여 R과 P를 매번 재생성한다. 동일 seed로 인해 결과적으로 같은 R, P가 생성되지만, 이는 seed 결정론에 의존하는 **암묵적** 보장이다.

본 실험에서는 R과 P를 **명시적으로 한 번만 생성**하고, α = {0, 0.001, 0.01, 0.1}만 바꿔서 결합한다. 이를 통해:

1. 동일한 R, P에서 α만의 순수한 영향을 관찰
2. 불필요한 재생성 제거 (특히 brute-force 풀이, posiform 생성)
3. 실험 의도가 코드에 명확히 드러남

## 현재 구조의 문제

```python
# 현재: α마다 전체 재생성
for alpha in [0, 0.001, 0.01, 0.1]:
    for run in range(num_instances):
        Q, info = create_qubo_hardened_posiform(
            n, ..., posiform_scale=alpha, seed=run * 53
        )
        # SA 실행
```

- `create_qubo_hardened_posiform` 내부에서 매번 R_i 생성 + brute-force 풀이 + posiform 생성
- seed가 같으므로 결과는 동일하지만, 코드만 보면 알 수 없음
- α 4개 × 인스턴스 10개 = 40회 호출 → 실제 필요한 건 10회

## 변경할 구조

```python
# 변경: R, P를 한 번만 생성 → α만 바꿔서 결합
for run in range(num_instances):
    # 1회만 생성
    R_list, target, P = generate_fixed_rp(n, coeff_type, seed=run * 53)

    for alpha in [0, 0.001, 0.01, 0.1]:
        # 결합만 수행
        Q = combine_rp(R_list, P, alpha)
        # SA 실행
```

## 구현 단계

### Step 1: R/P 분리 생성 함수

`qubo_posiform_hardened.py`의 `create_qubo_hardened_posiform()`을 분해하여, R과 P를 별도로 반환하는 함수를 새 실험 파일 내에 구현한다. 기존 생성기 코드는 수정하지 않는다.

```python
def generate_components(n, max_subgraph_size=15, coeff_type='lin2',
                        posiform_coeff_range=(1.0, 1.0), seed=None):
    """
    R_i 리스트와 P를 분리 생성.

    Returns:
        R_list: [R_0, R_1, ...] — 각 subgraph의 random QUBO dict
        target: str — 각 R_i 최적해를 concatenate한 비트스트링
        Q_posiform: dict — posiform QUBO (스케일링 전)
    """
```

내부 로직은 `create_qubo_hardened_posiform()`의 Step 1~5를 그대로 따르되, Step 6 (결합)은 수행하지 않는다.

### Step 2: 결합 함수

```python
def combine_rp(R_list, Q_posiform, alpha):
    """
    Q_final = Σ R_i + α × P

    Returns:
        Q_final: dict
    """
```

### Step 3: 실험 본체

```python
def run_fixed_rp_experiment():
    """
    실험 설정:
    - N = 500, max_subgraph_size = 15
    - coeff_type: lin2, lin20
    - alpha: [0, 0.001, 0.01, 0.1]
    - SA: sweep = [100, 500, 1000, 5000], num_reads = 50
    - 인스턴스: 10개

    각 인스턴스에서 R, P를 1회 생성 → 4개 α에 대해 결합 + SA 실행.
    """
```

### Step 4: 측정 항목

인스턴스별 + α별로 다음을 기록:

| 항목 | 설명 |
|---|---|
| target_rate | SA 해가 target과 정확히 일치하는 비율 (%) |
| hamming_med | SA 해와 target 사이 Hamming distance 중앙값 |
| hamming_avg | SA 해와 target 사이 Hamming distance 평균 |
| best_energy | SA가 찾은 최저 에너지 |
| target_energy | Q(target) |
| gs_rate | SA가 target 에너지에 도달한 비율 (%) |

### Step 5: 검증

R, P 고정이 올바른지 확인하는 sanity check:
- 각 인스턴스에서 α=0.1일 때 target이 ground state인지 (single-flip delta > 0)
- 동일 인스턴스의 α=0에서 Σ R_i(target) 에너지가 α에 관계없이 동일한지

### Step 6: 출력 형식

```
============================================================
Fixed R/P α 비교 실험
N=500, coeff_type=lin2, instances=10, reads=50
============================================================

[Instance 0] target=10110...
  α=0     | sw=1000: target=0.0%, hamming_med=28, gs=78.0%
  α=0.001 | sw=1000: target=8.0%, hamming_med=12, gs=8.0%
  α=0.01  | sw=1000: target=20.0%, hamming_med=5, gs=20.0%
  α=0.1   | sw=1000: target=100%, hamming_med=0, gs=100%

... (인스턴스별 반복)

[종합 - sweep=5000]
  α    | target_rate | hamming_med | hamming_avg | gs_rate
  0    |       0.0%  |          28 |        29.2 |  80.2%
  0.001|      13.4%  |           8 |         9.2 |  13.4%
  0.01 |      29.0%  |           4 |         5.1 |  29.0%
  0.1  |     100.0%  |           0 |         0.0 | 100.0%
```

## 파일

- 새 파일: `hardened_posiform/experiment_fixed_rp.py`
- 기존 파일 수정 없음

## 기존 실험과의 차이

| 항목 | 기존 실험 | 본 실험 |
|---|---|---|
| R, P 생성 | α마다 재생성 (seed로 동일성 보장) | 인스턴스당 1회 생성, 명시적 공유 |
| α 범위 | 실험마다 다름 | 0, 0.001, 0.01, 0.1 (고정) |
| 비교 방식 | α 간 비교 시 seed 동일성에 의존 | 동일 객체 참조로 보장 |
| 계산 효율 | α × instances회 생성 | instances회 생성 |

## 예상 결과

기존 실험과 수치적으로 동일한 결과가 나와야 한다 (같은 seed 사용 시). 차이가 있다면 기존 구현에 seed 결정론이 깨지는 버그가 있다는 의미.
