# rho(k) for the greedy duel sequence Du(k,0)

k players, vectors c of length J=8. The lexicographically smallest c picks next
(ties to lowest index). Pick e adds (-1)^j C(e,j) to c[j].

## Files
- `duel.c`: exact 256-bit engine plus periodicity analyser (one k per process).
- `run.py`: runs k=3..50 on all cores, writes `results.csv` and `chart.png` after each k.
- `reference.py`: literal pure-Python implementation (any J), used for validation.
- `validate.py`: k=2 vs Thue-Morse, C vs Python pick-for-pick (k=2..12), known (mu, lambda).

## Method
- c[0] counts picks, so play goes in rounds of k. Within a round only the picker's
  vector changes and it leaves the candidate set, so each round's order is the sort
  of c[1..7] at the round start. No priority queue needed.
- Periodicity: KMP prefix function over the reversed sequence of rounds finds the
  longest periodic tail, then the exact pick-level lambda (smallest divisor of the
  round period that works) and mu (backward scan).
- Accept only if the tail is at least 5 lambda long (first cycle plus 4 repeats) and
  the same (mu, lambda) is found at two checkpoints where the second is at least
  double the first.
- Any full tie between two vectors (after round 0, where all start at zero) marks the
  k as `full tie ... (J too small)` and is never reported as settled.

## Run
    gcc -O3 -march=native -o duel duel.c
    python3 validate.py
    python3 run.py --cap 2e8                 # all k
    python3 run.py --only 42,50 --cap 1e9    # recompute some k, merge into results.csv
