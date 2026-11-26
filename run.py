"""Compute rho(k) for Du(k,0), k = KMIN..KMAX, one k per process, smallest k first.

Usage: python3 run.py [--cap 2e8] [--kmin 3] [--kmax 50] [--workers N]
Writes results.csv and chart.png after every finished k.
"""
import argparse
import csv
import math
import os
import subprocess
from fractions import Fraction
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "duel")
FIELDS = ["k", "mu", "lambda", "lead_rounds", "cycle_rounds", "picks_checked", "status",
          "ties", "first_tie"]


def fmt_cap(n):
    e = int(math.floor(math.log10(n)))
    m = n / 10 ** e
    return f"{m:g}e{e}"


def solve(args):
    k, cap = args
    out = subprocess.run([BIN, str(k), str(int(cap))], capture_output=True, text=True, check=True).stdout
    r = dict(kv.split("=") for kv in out.split("RESULT", 1)[1].split())
    mu, lam, picks = int(r["mu"]), int(r["lambda"]), int(r["picks"])
    row = {"k": k, "picks_checked": picks, "ties": r["ties"],
           "first_tie": r["first_tie"] if r["first_tie"] != "-1" else ""}
    if r["status"] == "settled":
        cyc = Fraction(lam, k)
        row.update(mu=mu, **{"lambda": lam}, lead_rounds=math.ceil(mu / k),
                   cycle_rounds=str(cyc) if cyc.denominator > 1 else cyc.numerator, status="settled")
    else:
        cand = f" (last candidate mu {r['cand_mu']} lambda {r['cand_lambda']})" if r["cand_mu"] != "-1" else ""
        status = {"unsettled": f"not settled by {fmt_cap(cap)}{cand}",
                  "overflow": f"stopped at {picks} (256-bit guard)",
                  "tie_J_too_small": f"full tie at pick {r['first_tie']} (J too small)"}[r["status"]]
        row.update(mu="", **{"lambda": ""}, lead_rounds="", cycle_rounds="", status=status)
    return row


def write_csv(rows, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for k in sorted(rows):
            w.writerow(rows[k])


def draw(rows, path, kmin, kmax):
    ok = [rows[k] for k in sorted(rows) if rows[k]["status"] == "settled"]
    bad = [rows[k] for k in sorted(rows) if rows[k]["status"] != "settled"]
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=150)
    ink, muted, grid = "#1f1f1f", "#6b6b6b", "#e6e6e6"
    blue, orange, gray = "#257de0", "#f47230", "#8a8a8a"
    if ok:
        # NaN at unsettled k so the lines break there instead of bridging the gap
        ks = sorted(rows)
        good = lambda k: rows[k]["status"] == "settled"
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
    ap.add_argument("--cap", type=float, default=2e8)
    ap.add_argument("--kmin", type=int, default=3)
    ap.add_argument("--kmax", type=int, default=50)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--out", default=HERE)
    ap.add_argument("--chart-only", action="store_true", help="redraw chart.png from results.csv")
    ap.add_argument("--only", default="", help="comma list of k to recompute, merged into existing results.csv")
    a = ap.parse_args()
    if not os.path.exists(BIN):
        subprocess.run(["gcc", "-O3", "-march=native", "-o", BIN, os.path.join(HERE, "duel.c")], check=True)
    csv_path, png_path = os.path.join(a.out, "results.csv"), os.path.join(a.out, "chart.png")
    rows = {}
    if a.only and os.path.exists(csv_path):  # keep earlier rows, recompute only the listed k
        for r in csv.DictReader(open(csv_path)):
            for f in ("k", "mu", "lambda", "lead_rounds", "picks_checked"):
                r[f] = int(r[f]) if r[f] else ""
            rows[r["k"]] = r
    if a.chart_only:
        rows = {}
        for r in csv.DictReader(open(csv_path)):
            r["k"] = int(r["k"])
            rows[r["k"]] = r
        draw(rows, png_path, a.kmin, a.kmax)
        return
    ks = [int(x) for x in a.only.split(",")] if a.only else range(a.kmin, a.kmax + 1)
    jobs = [(k, a.cap) for k in sorted(ks)]
    with Pool(a.workers) as pool:
        for row in pool.imap_unordered(solve, jobs):
            rows[row["k"]] = row
            write_csv(rows, csv_path)
            draw(rows, png_path, a.kmin, a.kmax)
            print(f"k={row['k']:>2} mu={row['mu']} lambda={row['lambda']} "
                  f"rho=({row['lead_rounds']}, {row['cycle_rounds']}) {row['status']}", flush=True)


if __name__ == "__main__":
    main()
