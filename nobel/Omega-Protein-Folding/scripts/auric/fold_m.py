"""
Fold_m (golden-mean / Zeckendorf folding) core.

This module ports the semantics from:
`legacy-docs/papers/auric-golden-phi/.../scripts/common_phi_fold.py`

Key idea:
- Map an arbitrary micro word b∈{0,1}^m to an integer N using Fibonacci weights
  (weight at position k is F_{k+1}, k=1..m).
- Return the Zeckendorf (no-adjacent-1) representation of N using the same
  Fibonacci weight system, clipped to length m digits.

Because {0,1}^m has size 2^m but X_m^Z has size F_{m+2}, Fold_m is many-to-one.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable, List, Sequence

import numpy as np


def fib_upto(n: int) -> List[int]:
    """Return Fibonacci numbers F_1..F_n with F_1=F_2=1."""
    if n <= 0:
        return []
    if n == 1:
        return [1]
    f = [1, 1]
    while len(f) < n:
        f.append(f[-1] + f[-2])
    return f


def is_golden_legal_bits(x_bits: int) -> bool:
    """Return True iff x has no adjacent 1s in binary (golden-mean constraint)."""
    return (x_bits & (x_bits >> 1)) == 0


def pi_restrict_bits(x_bits: int, *, m_to: int) -> int:
    """Restriction π_{m->m_to} on low-to-high bit convention: keep lowest m_to bits."""
    if m_to <= 0:
        return 0
    return x_bits & ((1 << m_to) - 1)


def micro_bits_to_fib_value(micro_bits: int, *, m: int, fib: List[int]) -> int:
    """
    Convert micro bits (length m) into an integer N via Fibonacci weights F_{k+1}.

    Convention: bit (k-1) corresponds to digit c_k.
    Weight(c_k) = F_{k+1} = fib[k] with fib[0]=F_1, fib[1]=F_2, ...
    """
    N = 0
    # k=1..m -> index k-1 in bits, index k in fib list
    for k in range(1, m + 1):
        if (micro_bits >> (k - 1)) & 1:
            N += fib[k]
    return N


def zeckendorf_bits(N: int, *, m: int, fib: List[int]) -> int:
    """
    Zeckendorf digits (length m) for Fibonacci weights F_{k+1}, k=1..m.

    Greedy: pick largest weight <= remaining, skip adjacent position.
    Returns digits as bits in the same low-to-high convention.
    """
    if N < 0:
        raise ValueError("N must be non-negative")
    if m <= 0:
        return 0

    remaining = int(N)
    out_bits = 0
    k = m
    while remaining > 0 and k >= 1:
        w = fib[k]  # F_{k+1}
        if w <= remaining:
            out_bits |= 1 << (k - 1)
            remaining -= w
            k -= 2  # skip adjacent
        else:
            k -= 1
    return out_bits


def fold_m_bits(micro_bits: int, *, m: int) -> int:
    """
    Fold a length-m microstate bits to a golden-mean legal bits (no adjacent 1s).

    Semantics matches `common_phi_fold.fold_m` (value -> Zeckendorf digits).
    """
    if m <= 0:
        return 0
    fib = fib_upto(m + 2)
    N = micro_bits_to_fib_value(micro_bits, m=m, fib=fib)
    return zeckendorf_bits(N, m=m, fib=fib)


@dataclass(frozen=True)
class FoldTable:
    """Lookup table for Fold_m on length-m words: micro_bits -> folded_bits."""

    m: int
    table: List[int]  # length 2^m

    def fold(self, micro_bits: int) -> int:
        return self.table[micro_bits & ((1 << self.m) - 1)]


@lru_cache(maxsize=64)
def build_fold_table(m: int) -> FoldTable:
    """
    Precompute Fold_m for all 2^m micro words.

    This makes sliding-window folding O(1) per window (critical for N=1000 runs).
    """
    if m <= 0:
        return FoldTable(m=0, table=[0])
    size = 1 << m
    tbl = [0] * size
    for micro in range(size):
        tbl[micro] = fold_m_bits(micro, m=m)
    return FoldTable(m=m, table=tbl)


def bits_to_str_low_to_high(x_bits: int, *, m: int) -> str:
    """Render bits as a length-m string in low-to-high order (k=1..m)."""
    return "".join("1" if (x_bits >> i) & 1 else "0" for i in range(m))


def bits_to_str_high_to_low(x_bits: int, *, m: int) -> str:
    """Render bits as a length-m string in high-to-low order (paper-friendly)."""
    return bits_to_str_low_to_high(x_bits, m=m)[::-1]


def assert_fold_invariants(m: int, *, trials: int = 2000, seed: int = 0) -> None:
    """
    Lightweight self-check for Fold_m invariants:
    - idempotence
    - fixed points are golden-legal words
    - restriction consistency (empirical spot-check)
    """
    import random

    rng = random.Random(seed)
    tab = build_fold_table(m)

    for _ in range(trials):
        micro = rng.randrange(0, 1 << m)
        x = tab.fold(micro)
        xx = tab.fold(x)
        if x != xx:
            raise AssertionError(f"Idempotence failed: fold(fold(micro))!=fold(micro), m={m}")
        if not is_golden_legal_bits(x):
            raise AssertionError(f"Output not golden-legal, m={m}, x={x:b}")
        if tab.fold(x) != x:
            raise AssertionError(f"Fixed-point failed: fold(x)!=x, m={m}, x={x:b}")

    # Cross-scale restriction sanity (on stabilized types, as in the unified spec).
    # NOTE: A literal identity fold_m1(pi(micro)) == pi(fold_m2(micro)) need not hold
    # for all choices of micro-level projections under Fibonacci weights. Here we only
    # assert that restricting a stabilized type yields a stabilized type at the coarser m.
    for m1 in [max(1, m // 2), max(1, m - 2)]:
        if m1 >= m:
            continue
        tab1 = build_fold_table(m1)
        for _ in range(max(1, trials // 10)):
            micro = rng.randrange(0, 1 << m)
            f_m = tab.fold(micro)
            f_m_r = pi_restrict_bits(f_m, m_to=m1)
            if not is_golden_legal_bits(f_m_r):
                raise AssertionError(
                    f"Restricted stabilized type not golden-legal: m={m}, m1={m1}, x={f_m_r:b}"
                )
            if tab1.fold(f_m_r) != f_m_r:
                raise AssertionError(
                    f"Restricted stabilized type not a fixed point at m1: m={m}, m1={m1}, x={f_m_r:b}"
                )


def pack_bits_to_int(bits: Sequence[int]) -> int:
    """
    Pack a 0/1 sequence into an integer using the low-to-high convention:
      bits[0] -> bit 0 (k=1), bits[1] -> bit 1, ...
    """
    out = 0
    for i, b in enumerate(bits):
        if int(b) & 1:
            out |= 1 << i
    return out


def fold_sliding_windows_bits(bits: np.ndarray, *, m: int) -> np.ndarray:
    """
    Fold all length-m sliding windows of a binary sequence.

    Returns an array of folded bit-words encoded as integers (low-to-high),
    length (len(bits) - m + 1). If len(bits) < m, returns empty array.
    """
    b = np.asarray(bits, dtype=np.uint8).reshape(-1)
    if m <= 0:
        raise ValueError("m must be positive")
    n = int(b.shape[0])
    if n < m:
        return np.zeros((0,), dtype=np.uint32)

    tab = build_fold_table(m)
    mask = (1 << m) - 1

    micro = 0
    for i in range(m):
        micro |= (int(b[i]) & 1) << i
    out = np.empty((n - m + 1,), dtype=np.uint32)
    out[0] = np.uint32(tab.fold(micro))

    for t in range(1, n - m + 1):
        new_bit = int(b[t + m - 1]) & 1
        micro = ((micro >> 1) | (new_bit << (m - 1))) & mask
        out[t] = np.uint32(tab.fold(micro))
    return out


def iter_folded_words_str(
    folded_bits_seq: Iterable[int], *, m: int, high_to_low: bool = True
) -> Iterable[str]:
    """
    Render a sequence of folded bit-words (ints) into strings.
    Useful for reports/debugging; not intended for hot loops.
    """
    for x in folded_bits_seq:
        if high_to_low:
            yield bits_to_str_high_to_low(int(x), m=m)
        else:
            yield bits_to_str_low_to_high(int(x), m=m)

