"""Offline 500-inst regeneration — extracts the active P16/Adv6.4 graph from the
existing 200-inst pkl (union of all instance nodes/edges = full active set, 5612n/40088e)
so no QPU call is needed. Same pipeline as regen_pegasus_instances.py.
"""
import os
import pickle
import random
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

import networkx as nx
from networkx.algorithms.community import kernighan_lin_bisection
from pysat.solvers import Minisat22

INST_DIR = "/home/yideun/qubo_dataset/hardened_posiform/instances"
SRC_PKL = os.path.join(INST_DIR, "instances_pegasus16_lin2_200.pkl")
OUT_PKL = os.path.join(INST_DIR, "instances_pegasus16_lin2_500.pkl")
INST_TOPO = "pegasus"
INST_TOPO_SIZE = 16
INST_COEFF = "lin2"
INST_MAX_SUB = 16
TARGET_COUNT = 500
NUM_WORKERS = os.cpu_count()

COEFF_LIN2 = [-1, 1]


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
    coeffs = COEFF_LIN2
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


def build_graph_from_pkl(pkl_path):
    """Reconstruct active QPU graph from union of all-instance R_sum/P edges."""
    with open(pkl_path, "rb") as f:
        d = pickle.load(f)
    G = nx.Graph()
    for inst in d["instances"]:
        for (i, j) in list(inst["R_sum"].keys()) + list(inst["P"].keys()):
            G.add_node(i)
            if i != j:
                G.add_node(j)
                G.add_edge(min(i, j), max(i, j))
    return G, d["meta"]


def main():
    if os.path.exists(OUT_PKL):
        sz = os.path.getsize(OUT_PKL)
        try:
            with open(OUT_PKL, "rb") as f:
                pickle.load(f)
            print(f"  [{OUT_PKL}] already exists and is valid ({sz/1024/1024:.1f} MB) — skipping")
            return
        except Exception as e:
            print(f"  [{OUT_PKL}] exists but corrupt ({sz/1024/1024:.1f} MB, {e}) — regenerating")

    print(f"Loading graph from {SRC_PKL} ...")
    G, src_meta = build_graph_from_pkl(SRC_PKL)
    print(f"  graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"  src qpu_solver: {src_meta.get('qpu_solver')}")

    t0 = time.perf_counter()
    worker = partial(_generate_single_instance, G=G, max_sub=INST_MAX_SUB, coeff=INST_COEFF)
    all_instances = [None] * TARGET_COUNT
    done = 0
    print(f"\nGenerating {TARGET_COUNT} instances on {NUM_WORKERS} workers ...")
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as ex:
        futures = {ex.submit(worker, i): i for i in range(TARGET_COUNT)}
        for fut in as_completed(futures):
            idx = futures[fut]
            all_instances[idx] = fut.result()
            done += 1
            if done % 25 == 0 or done == TARGET_COUNT:
                elapsed = time.perf_counter() - t0
                eta = elapsed / done * (TARGET_COUNT - done)
                print(f"    {done}/{TARGET_COUNT}  ({elapsed:.0f}s, ETA {eta:.0f}s)")
    gen_t = time.perf_counter() - t0

    qpu_nodes = set(G.nodes())
    qpu_edges = set((min(u, v), max(u, v)) for u, v in G.edges())
    failures = []
    for idx, inst in enumerate(all_instances):
        nodes = set()
        edges = set()
        for (i, j) in list(inst["R_sum"].keys()) + list(inst["P"].keys()):
            nodes.add(i)
            if i != j:
                nodes.add(j)
                edges.add((min(i, j), max(i, j)))
        if not nodes.issubset(qpu_nodes) or not edges.issubset(qpu_edges):
            failures.append((idx, len(nodes - qpu_nodes), len(edges - qpu_edges)))
    if failures:
        raise RuntimeError(f"QPU compat verification failed: {len(failures)}/{TARGET_COUNT}")
    print(f"  ✓ verify passed: {TARGET_COUNT}/{TARGET_COUNT}")

    save_obj = {
        "meta": {
            "topology": INST_TOPO,
            "topo_size": INST_TOPO_SIZE,
            "coeff": INST_COEFF,
            "max_sub_size": INST_MAX_SUB,
            "num_instances": TARGET_COUNT,
            "n": all_instances[0]["n"],
            "graph_nodes": G.number_of_nodes(),
            "graph_edges": G.number_of_edges(),
            "qpu_verified": True,
            "qpu_solver": src_meta.get("qpu_solver"),
            "active_qubits": G.number_of_nodes(),
            "active_couplers": G.number_of_edges(),
            "qpu_compat_verified_at": src_meta.get("qpu_compat_verified_at"),
            "regen_source": "offline from " + os.path.basename(SRC_PKL),
        },
        "instances": all_instances,
    }
    tmp_path = OUT_PKL + ".tmp"
    with open(tmp_path, "wb") as f:
        pickle.dump(save_obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp_path, OUT_PKL)
    sz = os.path.getsize(OUT_PKL) / 1024 / 1024
    print(f"\n  saved {OUT_PKL} ({sz:.1f} MB, gen {gen_t:.0f}s)")


if __name__ == "__main__":
    main()
