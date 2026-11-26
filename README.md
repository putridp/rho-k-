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
- Settle rule: a checkpoint proposes a candidate when the tail is at least 5 lambda
  long, and the candidate is settled only when `certify.py` proves it permanent.
  Repeat counts are no proof: k=39's (1002830, 216216) held 19668 repeats before
  breaking. `repeats_seen` (engine repeats at detection) is a diagnostic only.
- Jump-ahead: when the certificate predicts a break, it also gives the exact state at
  the start of the break round (closed form below). The engine resumes there with
  the picks just before it as a prefix (so a new loop's mu can start before the jump
  point), instead of simulating the false loop. `false_loops` in results.csv lists
  every candidate that broke, with its exact break pick.
- Width: c[7] grows like e^8 / (8! k); the engine guards it and run.py switches to the
  320-bit build (`duel320`, `-DL=5`) for runs past about 2e10 picks.
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
- Checked against every false loop seen so far: the predicted break pick equals the
  engine's BREAK pick every time (k=15, 16, 24, 39, 48), and the jumped-to state equals
  direct simulation (validate.py).
- Full run: all k=3..50 settled by certificate. 84 false loops were found and jumped
  over on the way, 33 of them for k=48, which settles at (8945831100, 4478976) after
  candidates that broke after 10 to 618 repeats. k=39 settles at (4253556902, 216216).

    python3 certify.py 15 30591 8100     # false loop: breaks after about 25 repeats, pick 239670
    python3 certify.py --settled         # every settled k -> certificates.csv

## Run
    gcc -O3 -march=native -o duel duel.c && gcc -O3 -march=native -DL=5 -o duel320 duel.c
    python3 validate.py
    python3 run.py                         # all k, certificate settle rule (about 15 min on 4 cores)
    python3 run.py --only 39,48            # recompute some k, merge into results.csv
