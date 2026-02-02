from __future__ import annotations

import random

import numpy as np

from fold_m import assert_fold_invariants, fold_sliding_windows_bits  # noqa: E402
from metrics import block_entropy, entropy_rate_block_slope, runlength_1, type_entropy, type_support  # noqa: E402
from readout import rho_B_from_npath  # noqa: E402


def _rand_bits(n: int, rng: random.Random) -> np.ndarray:
    return np.asarray([rng.randrange(0, 2) for _ in range(n)], dtype=np.uint8)


def main() -> None:
    # Fold invariants
    for m in [4, 6, 8, 10]:
        assert_fold_invariants(m, trials=500, seed=0)

    rng = random.Random(0)

    # Sliding-window folding sanity
    bits = _rand_bits(200, rng)
    folded = fold_sliding_windows_bits(bits, m=6)
    assert folded.shape[0] == 200 - 6 + 1
    assert type_support(folded) >= 1
    assert np.isfinite(type_entropy(folded))

    # Metrics sanity
    rs = runlength_1(bits)
    assert rs.run1_max >= 0
    h4 = block_entropy(bits, 4)
    assert np.isfinite(h4)
    h_hat = entropy_rate_block_slope(bits, k_min=4, k_max=8)
    assert np.isfinite(h_hat)

    # rho_B determinism sanity
    n_path = np.zeros((20, 6), dtype=np.int32)
    for t in range(19):
        n_path[t + 1, 0] = n_path[t, 0] + 1  # cumulative +e0 steps
    b = rho_B_from_npath(n_path)
    assert b.shape[0] == 19 and int(np.sum(b)) == 19

    print("[auric.selfcheck] ok")


if __name__ == "__main__":
    main()

