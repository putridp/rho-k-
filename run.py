"""Compute rho(k) for Du(k,0), k = KMIN..KMAX, one k per process, smallest k first.

Usage: python3 run.py [--cap 1e12] [--store 2e9] [--kmin 3] [--kmax 50] [--workers N]

Settle rule: a candidate loop (mu, lambda) is settled only when certify.py proves it
permanent. Per k: the engine runs until its first candidate, certify.py proves it or
predicts the exact break, and on a break the engine jumps straight to the break round
using the exact state from the closed form, then looks for the next candidate.
--cap is the absolute pick limit per k, --store the picks kept per engine run (about
1.1 bytes each). 320-bit engine (duel320) is used automatically past about 2e10 picks.
repeats_seen (engine repeats seen at detection) is a diagnostic only.
Writes results.csv and chart.png after every finished k.
"""
import argparse
import csv
import math
import os
import subprocess
import tempfile
from fractions import Fraction
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import certify

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "duel")
FIELDS = ["k", "mu", "lambda", "lead_rounds", "cycle_rounds", "picks_checked", "repeats_seen",
          "status", "ties", "first_tie", "false_loops"]


def fmt_cap(n):
    e = int(math.floor(math.log10(n)))
    m = n / 10 ** e
    return f"{m:g}e{e}"


def find_candidate(k, cap, store, resume):
    """Run the engine from pick 0 or a jump state until its first candidate loop."""
    with tempfile.TemporaryDirectory() as td:
        start = resume[0] - resume[1] if resume else 0
        end = min(cap, start + store)
        binary, extra = certify.resume_args(k, resume, end, td)
        out = subprocess.run([binary, str(k), str(int(end)), "1024", "0", str(int(store)), "-1", "0"] + extra,
                             capture_output=True, text=True, check=True, cwd=td).stdout
    return dict(kv.split("=") for kv in out.split("RESULT", 1)[1].split())


def solve(args):
    k, cap, store = args
    resume, history = None, []
    row = {"k": k, "mu": "", "lambda": "", "lead_rounds": "", "cycle_rounds": "", "repeats_seen": ""}
    while True:
        r = find_candidate(k, cap, store, resume)
        row.update(picks_checked=int(r["picks"]), ties=r["ties"],
                   first_tie=r["first_tie"] if r["first_tie"] != "-1" else "")
        if r["status"] != "candidate":
            last = f" (last loop {history[-1]})" if history else ""
            row["status"] = {"unsettled": f"no loop found by pick {r['picks']}{last}",
                             "overflow": f"stopped at {r['picks']} (overflow guard){last}",
                             "tie_J_too_small": f"full tie at pick {r['first_tie']} (J too small)"}[r["status"]]
            break
        mu, lam = int(r["mu"]), int(r["lambda"])
        cert = certify.certify(k, mu, lam, verbose=False, resume=resume)
        print(f"k={k} candidate ({mu}, {lam}) seen {r['repeats']} repeats: {cert['verdict']}"
              + (f", break at pick {cert['break_pick']}" if cert["break_pick"] != "" else ""), flush=True)
        if cert["verdict"] == "permanent":
            cyc = Fraction(lam, k)
            row.update(mu=mu, **{"lambda": lam}, lead_rounds=math.ceil(mu / k), repeats_seen=int(r["repeats"]),
                       cycle_rounds=str(cyc) if cyc.denominator > 1 else cyc.numerator, status="settled")
            break
        history.append(f"({mu}, {lam}) broke at pick {cert['break_pick']} after {cert['repeats']} repeats")
        row["picks_checked"] = cert["break_pick"]
        if cert["break_pick"] >= cap:
            row["status"] = f"not settled by {fmt_cap(cap)} (last loop {history[-1]})"
            break
        resume = cert["jump"]
    row["false_loops"] = "; ".join(history)
    return row


def write_csv(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for k in sorted(rows):
            w.writerow(rows[k])


def draw(rows, path, kmin, kmax):
    is_good = lambda r: r["status"] == "settled" or r["status"].startswith("certified permanent")
    ok = [rows[k] for k in sorted(rows) if is_good(rows[k])]
    bad = [rows[k] for k in sorted(rows) if not is_good(rows[k])]
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=150)
    ink, muted, grid = "#1f1f1f", "#6b6b6b", "#e6e6e6"
    blue, orange, gray = "#257de0", "#f47230", "#8a8a8a"
    if ok:
        # NaN at unsettled k so the lines break there instead of bridging the gap
        ks = sorted(rows)
        good = lambda k: is_good(rows[k])
        lead = [max(int(rows[k]["lead_rounds"]), 0.5) if good(k) else float("nan") for k in ks]
        cyc = [float(Fraction(str(rows[k]["cycle_rounds"]))) if good(k) else float("nan") for k in ks]
        ax.plot(ks, lead, "-o", color=blue, lw=2, ms=6, mec="white", mew=1.5, label="lead-in (rounds)")
        ax.plot(ks, cyc, "-s", color=orange, lw=2, ms=6, mec="white", mew=1.5, label="cycle (rounds)")
    ax.set_yscale("log")
    lo, hi = ax.get_ylim() if ok else (1, 1e6)
    top = hi * 3
    if bad:
        ax.scatter([r["k"] for r in bad], [top] * len(bad), marker="^", s=70, color=gray,
                   edgecolor="white", linewidth=1.5, zorder=3, label="not settled / flagged")
        for r in bad:
            ax.annotate(str(r["k"]), (r["k"], top), textcoords="offset points", xytext=(0, 8),
                        ha="center", fontsize=8, color=muted)
    ax.set_ylim(min(lo, 0.5), top * 4)
    ax.set_xlim(kmin - 1, kmax + 1)
    ax.set_xlabel("k (players)", color=ink)
    ax.set_ylabel("rounds (log scale)", color=ink)
    ax.set_title("ρ(k) for the greedy duel Du(k,0): lead-in and cycle length", color=ink, loc="left")
    ax.grid(True, which="major", color=grid, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(grid)
    ax.tick_params(colors=muted)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=float, default=1e12)
    ap.add_argument("--store", type=float, default=2e9)
    ap.add_argument("--kmin", type=int, default=3)
    ap.add_argument("--kmax", type=int, default=50)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--out", default=HERE)
    ap.add_argument("--chart-only", action="store_true", help="redraw chart.png from results.csv")
    ap.add_argument("--only", default="", help="comma list of k to recompute, merged into existing results.csv")
    a = ap.parse_args()
    for b, flags in ((BIN, []), (BIN + "320", ["-DL=5"])):
        if not os.path.exists(b):
            subprocess.run(["gcc", "-O3", "-march=native", *flags, "-o", b, os.path.join(HERE, "duel.c")], check=True)
    csv_path, png_path = os.path.join(a.out, "results.csv"), os.path.join(a.out, "chart.png")
    rows = {}
    if a.only and os.path.exists(csv_path):  # keep earlier rows, recompute only the listed k
        for r in csv.DictReader(open(csv_path)):
            for f in ("k", "mu", "lambda", "lead_rounds", "picks_checked", "repeats_seen"):
                r[f] = int(r[f]) if r.get(f) else ""
            rows[r["k"]] = r
    if a.chart_only:
        rows = {}
        for r in csv.DictReader(open(csv_path)):
            r["k"] = int(r["k"])
            rows[r["k"]] = r
        draw(rows, png_path, a.kmin, a.kmax)
        return
    ks = [int(x) for x in a.only.split(",")] if a.only else range(a.kmin, a.kmax + 1)
    jobs = [(k, a.cap, a.store) for k in sorted(ks)]
    with Pool(a.workers) as pool:
        for row in pool.imap_unordered(solve, jobs):
            rows[row["k"]] = row
            write_csv(rows, csv_path)
            draw(rows, png_path, a.kmin, a.kmax)
            print(f"k={row['k']:>2} mu={row['mu']} lambda={row['lambda']} "
                  f"rho=({row['lead_rounds']}, {row['cycle_rounds']}) {row['status']}"
                  + (f"; false loops: {row['false_loops']}" if row["false_loops"] else ""), flush=True)


if __name__ == "__main__":
    main()
