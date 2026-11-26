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

    # 5. drift certificate: false k=15 loop must predict the engine's break pick exactly,
    #    the k=3 and real k=15 loops must be certified permanent.
    import certify
    r = certify.certify(15, 30591, 8100, verbose=False)
    print(f"certify k=15 (30591, 8100): {r['verdict']}, break pick {r['break_pick']} (engine: 239670)")
    ok &= r["break_pick"] == 239670
    for k, mu, lam in ((3, 47, 54), (15, 330741, 48600)):
        r = certify.certify(k, mu, lam, verbose=False)
        print(f"certify k={k} ({mu}, {lam}): {r['verdict']}")
        ok &= r["verdict"] == "permanent"

    # 6. k=39's (1002830, 216216) held 1000 repeats but is false: the certificate must
    #    predict the engine's break at pick 4253557401 (after 19668 repeats).
    r = certify.certify(39, 1002830, 216216, verbose=False)
    print(f"certify k=39 (1002830, 216216): {r['verdict']}, break pick {r['break_pick']} (engine: 4253557401)")
    ok &= r["break_pick"] == 4253557401

    # 7. jump-ahead: the closed-form state at the k=16 break round must equal direct
    #    simulation, and resuming the engine there (256- and 320-bit) must reproduce
    #    the directly simulated picks.
    r = certify.certify(16, 486295, 22528, verbose=False)
    X, P, v, prefix = r["jump"]
    n = 10 ** 6
    with tempfile.TemporaryDirectory() as td:
        subprocess.run([BIN, "16", str(X + n), "1e12", str(n), str(X + n), str(X)], check=True,
                       capture_output=True, cwd=td)
        ref_state = [[certify.from_hex(h) for h in l.split()] for l in open(os.path.join(td, "state_dump.txt"))]
        ref = open(os.path.join(td, "picks_dump.bin"), "rb").read()
    ok &= v == ref_state
    print(f"jump-ahead k=16 to pick {X}: closed-form state == simulated: {v == ref_state}")
    for binary, L in ((BIN, 4), (BIN + "320", 5)):
        with tempfile.TemporaryDirectory() as td:
            certify.write_resume(os.path.join(td, "resume"), X, P, v, prefix, L)
            subprocess.run([binary, "16", str(X + n), "1e12", str(n), str(X + n), str(X), "1000",
                            os.path.join(td, "resume")], check=True, capture_output=True, cwd=td)
            got = open(os.path.join(td, "picks_dump.bin"), "rb").read()
        print(f"  resumed {os.path.basename(binary)}: {len(got)} picks identical to direct run: {got == ref}")
        ok &= got == ref

    # 8. full driver on k=39: false loop certified as breaking, jump, real loop certified.
    import run
    row = run.solve((39, 1e12, 2e9))
    print(f"run.solve k=39: ({row['mu']}, {row['lambda']}) {row['status']}; false loops: {row['false_loops']}")
    ok &= (row["mu"], row["lambda"], row["status"]) == (4253556902, 216216, "settled")

    print("VALIDATION", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
