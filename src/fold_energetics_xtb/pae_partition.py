from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple


@dataclass(frozen=True)
class RigidDomain:
    domain_id: int
    residue_indices: List[int]
    mean_plddt: float
    intra_domain_mean_pae: float
    density_score: float


@dataclass(frozen=True)
class PAEPartitionResult:
    domains: List[RigidDomain]
    interface_cross_pae: Dict[Tuple[int, int], float]
    disordered_linker_residues: List[int]
    modularity_score: float


class PAERigidDomainPartitioner:
    """Partitions AlphaFold 3 Predicted Aligned Error (PAE) matrices into rigid structural domains."""

    def __init__(
        self,
        min_plddt: float = 70.0,
        max_pae_angstrom: float = 5.0,
        min_domain_size: int = 3,
        kernel_scale: float = 3.5,
    ) -> None:
        if min_plddt < 0.0 or min_plddt > 100.0:
            raise ValueError("min_plddt must be in [0, 100]")
        if max_pae_angstrom <= 0.0:
            raise ValueError("max_pae_angstrom must be positive")
        if min_domain_size < 1:
            raise ValueError("min_domain_size must be >= 1")
        self.min_plddt = min_plddt
        self.max_pae_angstrom = max_pae_angstrom
        self.min_domain_size = min_domain_size
        self.kernel_scale = max(kernel_scale, 1e-6)

    def partition(
        self,
        pae_matrix: Sequence[Sequence[float]],
        plddt_scores: Sequence[float],
    ) -> PAEPartitionResult:
        n = len(plddt_scores)
        if n == 0:
            raise ValueError("plddt_scores cannot be empty")
        if len(pae_matrix) != n or any(len(row) != n for row in pae_matrix):
            raise ValueError(f"pae_matrix must be square of shape ({n}, {n})")

        for row in pae_matrix:
            for val in row:
                if math.isnan(val) or val < 0.0:
                    raise ValueError("pae_matrix contains negative or NaN values")

        disordered = [i for i, score in enumerate(plddt_scores) if score < self.min_plddt]
        ordered_set = {i for i in range(n) if plddt_scores[i] >= self.min_plddt}

        if not ordered_set:
            return PAEPartitionResult(
                domains=[],
                interface_cross_pae={},
                disordered_linker_residues=disordered,
                modularity_score=0.0,
            )

        # Build symmetric similarity adjacency among ordered residues
        adj: Dict[int, Dict[int, float]] = {i: {} for i in ordered_set}
        total_edge_weight = 0.0
        for i in sorted(ordered_set):
            for j in sorted(ordered_set):
                if i >= j:
                    continue
                sym_pae = 0.5 * (pae_matrix[i][j] + pae_matrix[j][i])
                if sym_pae <= self.max_pae_angstrom:
                    w = math.exp(-(sym_pae * sym_pae) / (2.0 * self.kernel_scale * self.kernel_scale))
                    adj[i][j] = w
                    adj[j][i] = w
                    total_edge_weight += w

        # Connected components over rigid PAE edges
        visited: set[int] = set()
        raw_clusters: List[List[int]] = []
        for node in sorted(ordered_set):
            if node in visited:
                continue
            comp: List[int] = []
            queue = [node]
            visited.add(node)
            while queue:
                curr = queue.pop(0)
                comp.append(curr)
                for nxt in sorted(adj[curr].keys()):
                    if nxt not in visited:
                        visited.add(nxt)
                        queue.append(nxt)
            comp.sort()
            if len(comp) >= self.min_domain_size:
                raw_clusters.append(comp)
            else:
                disordered.extend(comp)

        disordered.sort()

        domains: List[RigidDomain] = []
        node_to_domain: Dict[int, int] = {}
        for dom_id, residues in enumerate(raw_clusters):
            for r in residues:
                node_to_domain[r] = dom_id
            mean_plddt = sum(plddt_scores[r] for r in residues) / len(residues)
            pair_paes: List[float] = []
            internal_edges = 0
            possible_edges = max(1, len(residues) * (len(residues) - 1) // 2)
            for idx_a, r_a in enumerate(residues):
                for r_b in residues[idx_a + 1 :]:
                    pair_paes.append(0.5 * (pae_matrix[r_a][r_b] + pae_matrix[r_b][r_a]))
                    if r_b in adj[r_a]:
                        internal_edges += 1
            intra_pae = sum(pair_paes) / len(pair_paes) if pair_paes else 0.0
            density = internal_edges / possible_edges
            domains.append(
                RigidDomain(
                    domain_id=dom_id,
                    residue_indices=residues,
                    mean_plddt=round(mean_plddt, 3),
                    intra_domain_mean_pae=round(intra_pae, 4),
                    density_score=round(density, 4),
                )
            )

        # Inter-domain cross-PAE summary
        interface_cross_pae: Dict[Tuple[int, int], float] = {}
        for idx_a in range(len(domains)):
            for idx_b in range(idx_a + 1, len(domains)):
                cross_vals = [
                    0.5 * (pae_matrix[r_a][r_b] + pae_matrix[r_b][r_a])
                    for r_a in domains[idx_a].residue_indices
                    for r_b in domains[idx_b].residue_indices
                ]
                interface_cross_pae[(idx_a, idx_b)] = round(sum(cross_vals) / len(cross_vals), 4)

        # Newman-Girvan modularity Q over weighted adjacency
        modularity = 0.0
        if total_edge_weight > 1e-9 and len(domains) > 0:
            two_m = 2.0 * total_edge_weight
            degrees = {u: sum(adj[u].values()) for u in ordered_set}
            for u in ordered_set:
                if u not in node_to_domain:
                    continue
                for v in ordered_set:
                    if v not in node_to_domain or node_to_domain[u] != node_to_domain[v]:
                        continue
                    a_uv = adj[u].get(v, 0.0)
                    modularity += a_uv - (degrees[u] * degrees[v]) / two_m
            modularity /= two_m

        return PAEPartitionResult(
            domains=domains,
            interface_cross_pae=interface_cross_pae,
            disordered_linker_residues=disordered,
            modularity_score=round(modularity, 5),
        )
