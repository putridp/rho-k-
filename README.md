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
- A checkpoint proposes a candidate when the tail is at least 5 lambda long. The
  candidate is then verified pick by pick against its first cycle. It is accepted
  only when the tail after mu reaches max(10 (mu + lambda), 1000 lambda) picks
  (`repeats_seen` in results.csv). A mismatch prints BREAK and analysis resumes at
  the next checkpoint. Short-tail candidates do break: k=16 (486295, 22528) breaks
  at pick 2188192 after 75 repeats and k=24 (607570, 25344) at pick 9091779 after
  334, both confirmed against the pure-Python reference.
- Verification reads only the first cycle, so it runs past the stored prefix
  (`--store`, about 1.1 bytes per pick per running k) up to `--cap` with no extra memory.
- Any full tie between two vectors (after round 0, where all start at zero) marks the
  k as `full tie ... (J too small)` and is never reported as settled.

## Drift certificate (`certify.py`)
Proves a loop permanent or predicts its exact break pick, from one loop of picks and
the exact state at its start (dumped by `duel` with STATE_AT).
- For the loop read from a round start x, F_p[m] = sum of C(y, m) over p's picks. By
  Vandermonde, after r loops c_p[j] gains (-1)^j sum_i A_i(r) F_p[j-i], A_0(r) = r.
  So for a pair (p, q), layers below the first unbalanced layer m (F_p[m] != F_q[m])
  never change, and layer m drifts by exactly (-1)^m (F_p[m] - F_q[m]) per loop.
- Each round's order is the sort at its start, so checking every adjacent pair of
  every round covers every pick. A pair decided below m is permanent; decided at or
  above m it is permanent if the drift favours the winner, else it flips once the
  gap is used up. All pairs permanent proves the loop never breaks.
- Checked against every false loop seen so far (k=15, 16, 24, 48: breaks after 13 to
  618 repeats): the predicted break pick equals the engine's BREAK pick every time.

    python3 certify.py 15 30591 8100     # false loop: breaks after about 25 repeats, pick 239670
    python3 certify.py --settled         # every settled k -> certificates.csv

## Run
    gcc -O3 -march=native -o duel duel.c
    python3 validate.py
    python3 run.py --cap 2e10 --store 2e9          # all k (k=42 needs about 1.5e10 picks)
    python3 run.py --only 42,50 --cap 2e10         # recompute some k, merge into results.csv
