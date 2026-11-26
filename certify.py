"""Drift certificate for a found loop (mu, lambda) of Du(k,0).

Usage:
    python3 certify.py K MU LAMBDA          # one loop
    python3 certify.py --settled            # every settled k in results.csv -> certificates.csv

Output: "permanent" (a proof that the loop repeats forever) or "breaks after about N
repeats" with the exact predicted break pick.

Method (all exact integers):
- Loop length Lam = lcm(lambda, k) picks, so each player picks Lam/k times per loop.
  Work from round starts x in [T, T + Lam), T = first round start >= mu.
- For the loop read from round start x, F_p[m] = sum of C(y, m) over p's picks y in
  [x, x + Lam). By Vandermonde, after r more loops
      c_p[j](x + r Lam) = c_p[j](x) + (-1)^j sum_i A_i(r) F_p[j - i],
  A_i(r) = sum_{r' < r} C(r' Lam, i), A_0(r) = r.
  Layer m "balances" for players p, q when F_p[m] = F_q[m]. If layers < m balance,
  layers < m of c_p - c_q never change, and layer m changes by exactly
  d = (-1)^m (F_p[m] - F_q[m]) per loop (the per-repeat drift).
- A round's order is the sort of the players at its start and lex order is total, so
  the loop holds at repeat r iff every adjacent pair (p, q) of every round still has
  p < q. For each pair: j* = layer that decides it now, m = first unbalanced layer.
    j* < m                 -> permanent (decided in a frozen layer)
    j* > m (layer m tied)  -> layer m decides from the next loop: permanent if d < 0,
                              else breaks at repeat 1
    j* = m                 -> permanent if d < 0, else flips once the gap
                              c_p[m] - c_q[m] (< 0) is eaten by d: repeat ceil(gap / d)
  "permanent" for every pair proves by induction that the loop never breaks.
- The earliest flip gives the break repeat R; the exact break pick comes from sorting
  that round at repeat R with the closed form above.
- Diagnostic: the round-start phase where the most leading layers balance for all
  players at once, and the drift spread of the first unbalanced layer there.
"""
import csv
import math
import os
import subprocess
import sys
import tempfile
import time
from collections import Counter
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "duel")
J = 8
SIGN = [(-1) ** j for j in range(J)]


def dump_loop(k, T, Lam):
    """Run the C engine to pick T + Lam; return (state c at T, picks [T, T + Lam))."""
    cap = T + Lam
    with tempfile.TemporaryDirectory() as td:
        subprocess.run([BIN, str(k), str(cap), str(cap), str(Lam), str(cap), str(T)],
                       capture_output=True, text=True, check=True, cwd=td)
        state = []
        for line in open(os.path.join(td, "state_dump.txt")):
            row = []
            for h in line.split():
                v = int(h, 16)
                row.append(v - (1 << 256) if v >> 255 else v)
            state.append(row)
        picks = open(os.path.join(td, "picks_dump.bin"), "rb").read()
    assert len(state) == k and len(picks) == Lam, (len(state), len(picks))
    return state, picks


def binoms(y):
    return [math.comb(y, m) for m in range(J)]


def step(B):
    """C(y, .) -> C(y + 1, .)"""
    for m in range(J - 1, 0, -1):
        B[m] += B[m - 1]


def A_coeffs(i, Lam):
    """a_t with C(r' Lam, i) = sum_t a_t C(r', t); then A_i(r) = sum_t a_t C(r, t + 1)."""
    vals = [math.comb(r * Lam, i) for r in range(i + 1)]
    a = []
    while vals:
        a.append(vals[0])
        vals = [vals[n + 1] - vals[n] for n in range(len(vals) - 1)]
    return a


def state_at_repeat(c, F, r, Acoef):
    A = [sum(a * math.comb(r, t + 1) for t, a in enumerate(Acoef[i])) for i in range(J)]
    return [[c[p][j] + SIGN[j] * sum(A[i] * F[p][j - i] for i in range(j + 1)) for j in range(J)]
            for p in range(len(c))]


def less(vp, p, vq, q):
    return (vp, p) < (vq, q)


def check_pair(cp, cq, Fp, Fq, p, q, Acoef):
    """Return (break_repeat or None, j_star, m, d)."""
    dc = [cp[j] - cq[j] for j in range(J)]
    js = next((j for j in range(J) if dc[j]), None)
    if not ((js is not None and dc[js] < 0) or (js is None and p < q)):
        raise RuntimeError(f"state does not match loop order for pair {p},{q}")
    dF = [Fp[m] - Fq[m] for m in range(J)]
    m = next((i for i in range(J) if dF[i]), None)
    if m is None:
        return None, js, None, 0
    d = SIGN[m] * dF[m]
    if js is not None and js < m:
        return None, js, m, d
    if js is None or js > m:
        return (None if d < 0 else 1), js, m, d
    if d < 0:
        return None, js, m, d
    gap = -dc[m]  # > 0
    r0, rem = divmod(gap, d)
    if rem:
        return r0 + 1, js, m, d
    # exact tie in layer m at repeat r0: higher layers decide there
    v = state_at_repeat([cp, cq], [Fp, Fq], r0, Acoef)
    return (r0 + 1 if less(v[0], p, v[1], q) else r0), js, m, d


def certify(k, mu, lam, verbose=True):
    t0 = time.time()
    Lam = lam * k // math.gcd(lam, k)
    T = -(-mu // k) * k
    c, w = dump_loop(k, T, Lam)
    Acoef = [A_coeffs(i, Lam) for i in range(J)]

    # F for the loop read from T
    F = [[0] * J for _ in range(k)]
    B = binoms(T)
    for i in range(Lam):
        Fp = F[w[i]]
        for m in range(J):
            Fp[m] += B[m]
        step(B)

    By, ByL = binoms(T), binoms(T + Lam)
    best_R, best = None, None
    jstar_hist, m_hist, case = Counter(), Counter(), Counter()
    best_phase = (-1, None, None)  # (balanced layers, phase, drift spread of next layer)
    for t in range(Lam // k):
        x = T + t * k
        rnd = w[t * k:(t + 1) * k]
        # global balance diagnostic at this phase
        bal = 0
        for m in range(J):
            if len({F[p][m] for p in range(k)}) == 1:
                bal += 1
            else:
                break
        if bal > best_phase[0]:
            spread = None
            if bal < J:
                vals = [SIGN[bal] * F[p][bal] for p in range(k)]
                spread = max(vals) - min(vals)
            best_phase = (bal, x, spread)
        # adjacent pairs of this round
        for a in range(k - 1):
            p, q = rnd[a], rnd[a + 1]
            r, js, m, d = check_pair(c[p], c[q], F[p], F[q], p, q, Acoef)
            jstar_hist[js] += 1
            m_hist[m] += 1
            if r is None:
                case["permanent"] += 1
            else:
                case["flips"] += 1
                if best_R is None or r < best_R:
                    best_R = r
                    best = (x, [row[:] for row in c], [row[:] for row in F], rnd)
                    binding = f"players {p}<{q} at pick {x + a}: decided in layer {js}, gap {c[p][js] - c[q][js] if js is not None else 0}, layer {m} drifts {d:+d} per loop"
        # advance through the round
        for y in range(k):
            p = rnd[y]
            cp, Fp = c[p], F[p]
            for j in range(J):
                cp[j] += SIGN[j] * By[j]
                Fp[j] += ByL[j] - By[j]
            step(By)
            step(ByL)

    res = {"k": k, "mu": mu, "lambda": lam, "loop_picks": Lam, "pairs": sum(case.values()),
           "deciding_layers": dict(sorted(jstar_hist.items(), key=lambda kv: (kv[0] is None, kv[0]))),
           "first_unbalanced_layers": dict(sorted(m_hist.items(), key=lambda kv: (kv[0] is None, kv[0]))),
           "best_phase": best_phase[1], "balanced_layers": best_phase[0],
           "drift_spread": best_phase[2], "seconds": round(time.time() - t0, 1)}
    if best_R is None:
        res.update(verdict="permanent", break_pick="", repeats="")
    else:
        x, cx, Fx, rnd = best
        v = state_at_repeat(cx, Fx, best_R, Acoef)
        order = sorted(range(k), key=lambda p: (v[p], p))
        pos = next(a for a in range(k) if order[a] != rnd[a])
        bp = x + best_R * Lam + pos
        res.update(verdict=f"breaks after about {(bp - mu) // lam} repeats", break_pick=bp,
                   repeats=(bp - mu) // lam, binding=binding)
    if verbose:
        report(res)
    return res


def report(r):
    print(f"k={r['k']} loop (mu {r['mu']}, lambda {r['lambda']}): {r['verdict'].upper()}"
          + (f" (predicted break at pick {r['break_pick']})" if r["break_pick"] != "" else ""))
    if r.get("binding"):
        print(f"  first flip: {r['binding']}")
    print(f"  phase with most balanced layers: pick {r['best_phase']}, layers 0..{r['balanced_layers'] - 1} "
          f"balance for all players, layer {r['balanced_layers']} drift spread {r['drift_spread']}")
    print(f"  {r['pairs']} adjacent-pair decisions; deciding layer histogram {r['deciding_layers']}; "
          f"first unbalanced layer per pair {r['first_unbalanced_layers']}  [{r['seconds']}s]", flush=True)


def _job(args):
    return certify(*args, verbose=False)


def main():
    if sys.argv[1:] == ["--settled"]:
        rows = [r for r in csv.DictReader(open(os.path.join(HERE, "results.csv"))) if r["status"] == "settled"]
        jobs = sorted(((int(r["k"]), int(r["mu"]), int(r["lambda"])) for r in rows),
                      key=lambda j: -j[2] * j[0])  # longest first
        out = []
        with Pool(os.cpu_count()) as pool:
            for res in pool.imap_unordered(_job, jobs):
                report(res)
                out.append(res)
        fields = ["k", "mu", "lambda", "verdict", "break_pick", "repeats", "loop_picks", "pairs",
                  "best_phase", "balanced_layers", "drift_spread", "deciding_layers", "first_unbalanced_layers",
                  "binding"]
        with open(os.path.join(HERE, "certificates.csv"), "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            wr.writeheader()
            for res in sorted(out, key=lambda r: r["k"]):
                wr.writerow(res)
    else:
        k, mu, lam = (int(float(a)) for a in sys.argv[1:4])
        certify(k, mu, lam)


if __name__ == "__main__":
    main()
