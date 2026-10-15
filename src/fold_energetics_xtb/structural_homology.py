from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

ALPHABET_3DI = "ACDEFGHIKLMNPQRSTVWY"


@dataclass(frozen=True)
class Residue3D:
    res_idx: int
    res_name: str
    ca_xyz: Tuple[float, float, float]
    cb_xyz: Tuple[float, float, float]


@dataclass(frozen=True)
class CatalyticTriadMatch:
    residues: Tuple[int, int, int]
    triad_type: str
    nucleophile_base_dist: float
    base_acid_dist: float
    catalytic_geometry_valid: bool


def _euclidean_dist(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


class Foldseek3DiHomologyScorer:
    """Encodes 3D pocket coordinates into Foldseek-like 3Di structural tokens and detects catalytic triads."""

    CANONICAL_TRIADS = {
        ("SER", "HIS", "ASP"): "Ser-His-Asp Serine Hydrolase Triad",
        ("CYS", "HIS", "ASN"): "Cys-His-Asn Cysteine Protease Triad",
        ("SER", "HIS", "GLU"): "Ser-His-Glu Esterase Triad",
    }

    def encode_3di_sequence(self, residues: Sequence[Residue3D]) -> str:
        if not residues:
            return ""
        tokens: List[str] = []
        n = len(residues)
        for i, res in enumerate(residues):
            ca_cb = _euclidean_dist(res.ca_xyz, res.cb_xyz)
            prev_dist = _euclidean_dist(residues[i - 1].ca_xyz, res.ca_xyz) if i > 0 else 3.8
            next_dist = _euclidean_dist(res.ca_xyz, residues[i + 1].ca_xyz) if i + 1 < n else 3.8
            # Find nearest non-adjacent tertiary contact
            min_tertiary = 20.0
            for j, other in enumerate(residues):
                if abs(i - j) >= 2:
                    d = _euclidean_dist(res.cb_xyz, other.cb_xyz)
                    if d < min_tertiary:
                        min_tertiary = d
            descriptor = (ca_cb * 1.7 + prev_dist * 0.9 + next_dist * 0.9 + min_tertiary * 1.3)
            bucket = int(round(descriptor * 2.5)) % len(ALPHABET_3DI)
            tokens.append(ALPHABET_3DI[bucket])
        return "".join(tokens)

    def align_3di(self, seq_a: str, seq_b: str, gap_penalty: float = -2.0) -> Tuple[float, float]:
        """Needleman-Wunsch global alignment returning (raw_score, normalized_identity)."""
        clean_a = seq_a.strip().upper()
        clean_b = seq_b.strip().upper()
        if not clean_a or not clean_b:
            return 0.0, 0.0
        if any(ch not in ALPHABET_3DI for ch in clean_a + clean_b):
            raise ValueError("Sequences must contain only valid 3Di tokens")

        m, n = len(clean_a), len(clean_b)
        dp = [[0.0] * (n + 1) for _ in range(m + 1)]
        for i in range(1, m + 1):
            dp[i][0] = i * gap_penalty
        for j in range(1, n + 1):
            dp[0][j] = j * gap_penalty

        for i in range(1, m + 1):
            ca = clean_a[i - 1]
            idx_a = ALPHABET_3DI.index(ca)
            for j in range(1, n + 1):
                cb = clean_b[j - 1]
                idx_b = ALPHABET_3DI.index(cb)
                if ca == cb:
                    sub = 3.5
                else:
                    diff = abs(idx_a - idx_b)
                    sub = 1.2 if diff <= 2 else -1.0
                dp[i][j] = max(
                    dp[i - 1][j - 1] + sub,
                    dp[i - 1][j] + gap_penalty,
                    dp[i][j - 1] + gap_penalty,
                )

        # Traceback optimal global alignment to compute true aligned sequence identity
        i, j = m, n
        matches = 0
        aligned_columns = 0
        while i > 0 or j > 0:
            if i > 0 and j > 0:
                ca = clean_a[i - 1]
                cb = clean_b[j - 1]
                if ca == cb:
                    sub = 3.5
                else:
                    diff = abs(ALPHABET_3DI.index(ca) - ALPHABET_3DI.index(cb))
                    sub = 1.2 if diff <= 2 else -1.0
                if abs(dp[i][j] - (dp[i - 1][j - 1] + sub)) < 1e-9:
                    if ca == cb:
                        matches += 1
                    aligned_columns += 1
                    i -= 1
                    j -= 1
                    continue
            if i > 0 and abs(dp[i][j] - (dp[i - 1][j] + gap_penalty)) < 1e-9:
                aligned_columns += 1
                i -= 1
            else:
                aligned_columns += 1
                j -= 1

        identity = matches / max(aligned_columns, 1)
        return round(dp[m][n], 4), round(identity, 4)

    def compute_distance_rmsd(self, pocket_a: Sequence[Residue3D], pocket_b: Sequence[Residue3D]) -> float:
        """Rotation/translation-invariant pairwise distance RMSD (dRMSD) between two pockets of equal length."""
        if len(pocket_a) != len(pocket_b):
            raise ValueError("Pockets must have equal residue count for dRMSD")
        n = len(pocket_a)
        if n < 2:
            return 0.0
        sq_err = 0.0
        pairs = 0
        for i in range(n):
            for j in range(i + 1, n):
                da = _euclidean_dist(pocket_a[i].ca_xyz, pocket_a[j].ca_xyz)
                db = _euclidean_dist(pocket_b[i].ca_xyz, pocket_b[j].ca_xyz)
                sq_err += (da - db) ** 2
                pairs += 1
        return round(math.sqrt(sq_err / max(pairs, 1)), 4)

    def detect_catalytic_triads(
        self,
        residues: Sequence[Residue3D],
        max_nu_his_dist: float = 4.8,
        max_his_acid_dist: float = 5.2,
    ) -> List[CatalyticTriadMatch]:
        matches: List[CatalyticTriadMatch] = []
        for nu in residues:
            nu_name = nu.res_name.upper()
            if nu_name not in ("SER", "CYS"):
                continue
            for his in residues:
                if his.res_name.upper() != "HIS" or his.res_idx == nu.res_idx:
                    continue
                d_nu_his = _euclidean_dist(nu.cb_xyz, his.cb_xyz)
                if not (2.0 <= d_nu_his <= max_nu_his_dist):
                    continue
                for acid in residues:
                    acid_name = acid.res_name.upper()
                    key = (nu_name, "HIS", acid_name)
                    if key not in self.CANONICAL_TRIADS or acid.res_idx in (nu.res_idx, his.res_idx):
                        continue
                    d_his_acid = _euclidean_dist(his.cb_xyz, acid.cb_xyz)
                    valid = 2.0 <= d_his_acid <= max_his_acid_dist
                    if valid:
                        matches.append(
                            CatalyticTriadMatch(
                                residues=(nu.res_idx, his.res_idx, acid.res_idx),
                                triad_type=self.CANONICAL_TRIADS[key],
                                nucleophile_base_dist=round(d_nu_his, 3),
                                base_acid_dist=round(d_his_acid, 3),
                                catalytic_geometry_valid=True,
                            )
                        )
        return matches
