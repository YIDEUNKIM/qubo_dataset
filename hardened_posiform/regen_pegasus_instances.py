"""Regenerate pegasus16_lin2_{100,200}.pkl on the LATEST Pegasus chip
(Advantage_system6.4 as of 2026-05-20). One-shot script — runs the same pipeline
the notebook's instance-generation cell does, but with the solver name pinned
explicitly so we never accidentally pick up the older 4.1 chip.

Run from anywhere (uses absolute paths). Idempotent w.r.t. existing pkls: if a
target file already exists it is left alone.
"""

import os
import pickle
import random
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

import numpy as np
import networkx as nx
import dwave_networkx as dnx
from networkx.algorithms.community import kernighan_lin_bisection
from pysat.solvers import Minisat22
from dwave.system import DWaveSampler

TOKEN = os.environ.get("DWAVE_API_TOKEN")  # set via env var; never hardcode
SOLVER_NAME = "Advantage_system6.4"   # latest pegasus chip on this account
INST_TOPO = "pegasus"
INST_TOPO_SIZE = 16
INST_COEFF = "lin2"
INST_MAX_SUB = 16
INST_COUNTS = [100, 200]
INST_DIR = "/home/yideun/qubo_dataset/hardened_posiform/instances"
NUM_WORKERS = os.cpu_count()

COEFF_LIN2 = [-1, 1]
COEFF_LIN20 = [round(-1 + 0.1 * i, 1) for i in range(21)]


def kl_recursive_bisection(G, max_sub_size=15):
    nodes = set(G.nodes())
    if len(nodes) <= max_sub_size:
        return [nodes]
    sub = G.subgraph(nodes).copy()
    try:
        A, B = kernighan_lin_bisection(sub)
    except nx.NetworkXError:
        parts = []
        for comp in nx.connected_components(sub):
            parts.extend(kl_recursive_bisection(G.subgraph(comp).copy(), max_sub_size))
        return parts
    result = []
    for part in [A, B]:
        if len(part) <= max_sub_size:
            result.append(part)
        else:
            result.extend(kl_recursive_bisection(G.subgraph(part).copy(), max_sub_size))
    return result


def gen_hardware_random_qubo(G, variables, coeff_type="lin2"):
    coeffs = COEFF_LIN2 if coeff_type == "lin2" else COEFF_LIN20
    var_set = set(variables)
    Q = {}
    for v in variables:
        Q[(v, v)] = random.choice(coeffs)
    sub = G.subgraph(var_set)
    for i, j in sub.edges():
        lo, hi = min(i, j), max(i, j)
        Q[(lo, hi)] = random.choice(coeffs)
    return Q


def solve_subproblem_brute_force(Q, variables):
    n = len(variables)
    if n > 23:
        raise ValueError(f"Brute force 불가: subgraph size {n} > 23")
    var_list = sorted(variables)
    var_to_idx = {v: i for i, v in enumerate(var_list)}

    linear_map, quad = {}, []
    for (i, j), w in Q.items():
        if i == j:
            linear_map[var_to_idx[i]] = linear_map.get(var_to_idx[i], 0) + w
        else:
            quad.append((var_to_idx[i], var_to_idx[j], w))
    linear = [(idx, w) for idx, w in linear_map.items() if w != 0]
    quad = [(a, b, w) for a, b, w in quad if w != 0]

    best_energy = float("inf")
    best_bits = 0
    num_deg = 0
    all_gs_bits = []

    for bits in range(1 << n):
        energy = 0.0
        for idx, w in linear:
            if (bits >> (n - 1 - idx)) & 1:
                energy += w
        for a, b, w in quad:
            if ((bits >> (n - 1 - a)) & 1) and ((bits >> (n - 1 - b)) & 1):
                energy += w
        if energy < best_energy - 1e-12:
            best_energy = energy
            best_bits = bits
            num_deg = 1
            all_gs_bits = [bits]
        elif abs(energy - best_energy) < 1e-12:
            num_deg += 1
            all_gs_bits.append(bits)

    def bits_to_assignment(bits):
        return {v: (bits >> (n - 1 - i)) & 1 for i, v in enumerate(var_list)}

    return (
        bits_to_assignment(best_bits),
        best_energy,
        num_deg,
        [bits_to_assignment(b) for b in all_gs_bits],
    )


def posiform_planting_hardware(G, target_assignment, alpha):
    nodes = sorted(target_assignment.keys())
    n = len(nodes)
    node_to_idx = {v: i for i, v in enumerate(nodes)}
    edges = [(min(u, v), max(u, v)) for u, v in G.subgraph(nodes).edges()]
    if not edges:
        return {}

    Q_posiform = {}
    clauses_cnf = []
    all_tuples = [(0, 0), (0, 1), (1, 0), (1, 1)]
    max_clauses = 1000 * n

    def add_posiform_term(i, j):
        target_tuple = (target_assignment[i], target_assignment[j])
        wrong_tuples = [t for t in all_tuples if t != target_tuple]
        wi, wj = random.choice(wrong_tuples)
        lo, hi = min(i, j), max(i, j)
        if wi == 0 and wj == 0:
            Q_posiform[(i, i)] = Q_posiform.get((i, i), 0) - alpha
            Q_posiform[(j, j)] = Q_posiform.get((j, j), 0) - alpha
            Q_posiform[(lo, hi)] = Q_posiform.get((lo, hi), 0) + alpha
        elif wi == 0 and wj == 1:
            Q_posiform[(j, j)] = Q_posiform.get((j, j), 0) + alpha
            Q_posiform[(lo, hi)] = Q_posiform.get((lo, hi), 0) - alpha
        elif wi == 1 and wj == 0:
            Q_posiform[(i, i)] = Q_posiform.get((i, i), 0) + alpha
            Q_posiform[(lo, hi)] = Q_posiform.get((lo, hi), 0) - alpha
        else:
            Q_posiform[(lo, hi)] = Q_posiform.get((lo, hi), 0) + alpha

        lit_i = (node_to_idx[i] + 1) if wi == 0 else -(node_to_idx[i] + 1)
        lit_j = (node_to_idx[j] + 1) if wj == 0 else -(node_to_idx[j] + 1)
        clauses_cnf.append([lit_i, lit_j])

    def check_uniqueness():
        with Minisat22() as solver:
            for clause in clauses_cnf:
                solver.add_clause(clause)
            blocking = []
            for node in nodes:
                idx = node_to_idx[node] + 1
                blocking.append(-idx if target_assignment[node] == 1 else idx)
            solver.add_clause(blocking)
            return not solver.solve()

    check_interval = max(1, n // 4)
    for step in range(max_clauses):
        i, j = random.choice(edges)
        add_posiform_term(i, j)
        if (step + 1) % check_interval == 0:
            if check_uniqueness():
                Q_posiform = {k: v for k, v in Q_posiform.items() if abs(v) > 1e-15}
                return Q_posiform

    Q_posiform = {k: v for k, v in Q_posiform.items() if abs(v) > 1e-15}
    return Q_posiform


def _generate_single_instance(inst, G, max_sub, coeff):
    random.seed(inst * 53)
    parts = kl_recursive_bisection(G, max_sub)
    random_qubos = []
    target = {}
    total_deg = 1
    all_block_gs = []
    for part in parts:
        variables = sorted(part)
        R = gen_hardware_random_qubo(G, variables, coeff)
        assignment, energy, deg, gs_assignments = solve_subproblem_brute_force(R, variables)
        target.update(assignment)
        random_qubos.append(R)
        total_deg *= deg
        all_block_gs.append(gs_assignments)
    P = posiform_planting_hardware(G, target, 1.0)
    R_sum = {}
    for R in random_qubos:
        for k, v in R.items():
            R_sum[k] = R_sum.get(k, 0) + v
    t_energy_r = sum(
        w * target.get(i, 0) * (target.get(j, 0) if i != j else 1)
        if i != j else w * target.get(i, 0)
        for (i, j), w in R_sum.items()
    )
    t_energy_p = sum(
        w * target.get(i, 0) * (target.get(j, 0) if i != j else 1)
        if i != j else w * target.get(i, 0)
        for (i, j), w in P.items()
    )
    sorted_nodes = sorted(target.keys())
    target_str = "".join(str(target[nd]) for nd in sorted_nodes)
    return {
        "R_sum": R_sum,
        "P": P,
        "target": target,
        "target_str": target_str,
        "sorted_nodes": sorted_nodes,
        "n": len(target),
        "t_energy_r": t_energy_r,
        "t_energy_p": t_energy_p,
        "total_degeneracy": total_deg,
        "seed": inst * 53,
        "all_block_gs": all_block_gs,
        "partitions": [sorted(p) for p in parts],
    }


def _extract_inst_topo(inst):
    nodes, edges = set(), set()
    for (i, j) in list(inst["R_sum"].keys()) + list(inst["P"].keys()):
        nodes.add(i)
        if i != j:
            nodes.add(j)
            edges.add((min(i, j), max(i, j)))
    return nodes, edges


def main():
    print(f"Connecting to {SOLVER_NAME} ...")
    qpu = DWaveSampler(token=TOKEN, solver=SOLVER_NAME)
    G_qpu = qpu.to_networkx_graph()
    qpu_nodes = set(G_qpu.nodes())
    qpu_edges = set((min(u, v), max(u, v)) for u, v in G_qpu.edges())
    topo_type = qpu.properties.get("topology", {}).get("type", "?")
    if topo_type != INST_TOPO:
        raise RuntimeError(f"unexpected topology mismatch: {topo_type} != {INST_TOPO}")
    print(f"  solver={qpu.solver.name}, topology={topo_type}, "
          f"active=({len(qpu_nodes)}, {len(qpu_edges)})")

    extra_meta = {
        "qpu_solver": qpu.solver.name,
        "active_qubits": len(qpu_nodes),
        "active_couplers": len(qpu_edges),
        "qpu_compat_verified_at": "2026-05-20",
    }

    os.makedirs(INST_DIR, exist_ok=True)
    for count in INST_COUNTS:
        fname = f"instances_{INST_TOPO}{INST_TOPO_SIZE}_{INST_COEFF}_{count}.pkl"
        fpath = os.path.join(INST_DIR, fname)
        if os.path.exists(fpath):
            print(f"  [{fname}] 이미 존재 — skip")
            continue

        print(f"\n  [{fname}] generating {count} instances on {qpu.solver.name} ...")
        t0 = time.perf_counter()
        worker = partial(_generate_single_instance, G=G_qpu, max_sub=INST_MAX_SUB, coeff=INST_COEFF)
        all_instances = [None] * count
        done = 0
        with ProcessPoolExecutor(max_workers=NUM_WORKERS) as ex:
            futures = {ex.submit(worker, i): i for i in range(count)}
            for fut in as_completed(futures):
                idx = futures[fut]
                all_instances[idx] = fut.result()
                done += 1
                if done % 10 == 0 or done == count:
                    elapsed = time.perf_counter() - t0
                    print(f"    {done}/{count}  ({elapsed:.0f}s)")
        gen_t = time.perf_counter() - t0

        # strict QPU compat verification
        failures = []
        for idx, inst in enumerate(all_instances):
            nodes, edges = _extract_inst_topo(inst)
            if not nodes.issubset(qpu_nodes) or not edges.issubset(qpu_edges):
                failures.append((idx, len(nodes - qpu_nodes), len(edges - qpu_edges)))
        if failures:
            print(f"  ✗ verify failed: {len(failures)}/{count}")
            for ex in failures[:5]:
                print(f"    inst[{ex[0]}]: missing {ex[1]} nodes, {ex[2]} edges")
            raise RuntimeError("QPU compat verification failed — pkl not saved")
        print(f"  ✓ verify passed: {count}/{count}")

        save_obj = {
            "meta": {
                "topology": INST_TOPO,
                "topo_size": INST_TOPO_SIZE,
                "coeff": INST_COEFF,
                "max_sub_size": INST_MAX_SUB,
                "num_instances": count,
                "n": all_instances[0]["n"],
                "graph_nodes": G_qpu.number_of_nodes(),
                "graph_edges": G_qpu.number_of_edges(),
                "qpu_verified": True,
                **extra_meta,
            },
            "instances": all_instances,
        }
        with open(fpath, "wb") as f:
            pickle.dump(save_obj, f)
        sz = os.path.getsize(fpath) / 1024 / 1024
        print(f"  saved {fpath} ({sz:.1f} MB, gen {gen_t:.0f}s)")


if __name__ == "__main__":
    main()
