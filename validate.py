"""Validation: C engine vs literal Python reference, k=2 Thue-Morse, known (mu, lambda)."""
import os
import subprocess
import sys
import tempfile

from reference import du_sequence, naive_mu_lambda, thue_morse

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "duel")
KNOWN = {3: (47, 54), 4: (46, 192), 9: (14468, 4860), 10: (5374, 1600), 24: (607570, 25344)}


def run_c(k, cap, start=1024, dump=0, cwd=None):
    out = subprocess.run([BIN, str(k), str(cap), str(start), str(dump)],
                         capture_output=True, text=True, check=True, cwd=cwd).stdout
    line = [l for l in out.splitlines() if l.startswith("RESULT")][0]
    return dict(kv.split("=") for kv in line.split()[1:])


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
        r = run_c(k, 2e8)
        got = (int(r["mu"]), int(r["lambda"]))
        print(f"k={k}: got {got} expected {(mu, lam)} status={r['status']} picks={r['picks']} ties={r['ties']}")
        ok &= got == (mu, lam)
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
