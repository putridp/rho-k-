"""Validation: C engine vs literal Python reference, k=2 Thue-Morse, known (mu, lambda)."""
import os
import subprocess
import sys
import tempfile

from reference import du_sequence, naive_mu_lambda, thue_morse

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "duel")
# k=24 was (607570, 25344) under the old 2-checkpoint rule; that cycle breaks at pick
# 9091779 after 334 repeats (confirmed with reference.py). Value below holds 1000 repeats.
KNOWN = {3: (47, 54), 4: (46, 192), 9: (14468, 4860), 10: (5374, 1600), 24: (50652467, 532224)}


def run_c(k, cap, start=1024, dump=0, cwd=None, want_err=False):
    p = subprocess.run([BIN, str(k), str(cap), str(start), str(dump)],
                       capture_output=True, text=True, check=True, cwd=cwd)
    line = [l for l in p.stdout.splitlines() if l.startswith("RESULT")][0]
    r = dict(kv.split("=") for kv in line.split()[1:])
    return (r, p.stderr) if want_err else r


def main():
    ok = True
    # 1. k=2 against Thue-Morse (python reference, several J)
    for J in (8, 10, 12):
        n = 2 ** J + 8
        picks, ties = du_sequence(2, n, J)
        tm = thue_morse(n)
        first_dev = next((i for i in range(n) if picks[i] != tm[i]), None)
        print(f"k=2 J={J}: matches Thue-Morse for first {first_dev if first_dev is not None else n} picks, "
              f"first full tie at pick {ties[0] if ties else None}")
        ok &= first_dev == 2 ** J and ties and ties[0] == 2 ** J

    # 2. C engine vs literal reference, picks identical
    with tempfile.TemporaryDirectory() as td:
        for k in range(2, 13):
            n = 6000 if k < 8 else 3000
            ref, _ = du_sequence(k, n)
            run_c(k, n, start=10 ** 9, dump=n, cwd=td)
            got = list(open(os.path.join(td, "picks_dump.bin"), "rb").read())
            same = got == ref[: len(got)]
            print(f"C vs reference k={k}: {len(got)} picks identical={same}")
            ok &= same

    # 3. known (mu, lambda); also naive cross-check of the analyser for small k
    for k, (mu, lam) in KNOWN.items():
        r = run_c(k, 1e9)
        got = (int(r["mu"]), int(r["lambda"]))
        print(f"k={k}: got {got} expected {(mu, lam)} status={r['status']} picks={r['picks']} "
              f"repeats={r['repeats']} ties={r['ties']}")
        ok &= got == (mu, lam) and r["status"] == "settled" and int(r["repeats"]) >= 1000

    # 3b. k=16 has a false cycle (486295, 22528) that breaks at pick 2188192 after 75
    #     repeats. The engine must catch the break and not report it.
    r, err = run_c(16, 2e8, want_err=True)
    got = (int(r["mu"]), int(r["lambda"]))
    caught = "BREAK k=16 mu=486295 lambda=22528 at pick 2188192 after 75 repeats" in err
    print(f"k=16: got {got} status={r['status']} repeats={r['repeats']} false cycle caught={caught}")
    ok &= got != (486295, 22528) and r["status"] == "settled" and caught
    for k in (3, 4):
        picks, _ = du_sequence(k, 3000)
        print(f"naive analyser k={k}: {naive_mu_lambda(picks)}")

    # 4. k=2 in C: J=8 cannot hold Thue-Morse past 2^8 picks, so the engine must
    #    flag a full tie at pick 256 and refuse to call the result settled.
    r = run_c(2, 1e6)
    print(f"k=2 C engine 1e6 picks: status={r['status']} first_tie={r['first_tie']}")
    ok &= r["status"] == "tie_J_too_small" and r["first_tie"] == "256"

    print("VALIDATION", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
