"""Pure-Python reference for the greedy duel sequence Du(k,0).

Literal implementation of the definition (lex-min scan every pick, exact
integers, any J). Slow, used only to validate the C engine and the k=2 case.
"""
from math import comb


def du_sequence(k, n_picks, J=8):
    """Return (picks, ties). ties = list of pick indices (after round 0) where
    the chosen player's vector was fully equal to another candidate's."""
    c = [[0] * J for _ in range(k)]
    picks, ties = [], []
    for e in range(n_picks):
        best = 0
        for i in range(1, k):
            if c[i] < c[best]:
                best = i
        if e >= k and any(c[i] == c[best] for i in range(k) if i != best):
            ties.append(e)
        picks.append(best)
        for j in range(J):
            c[best][j] += (-1) ** j * comb(e, j)
    return picks, ties


def thue_morse(n):
    return [bin(i).count("1") & 1 for i in range(n)]


def naive_mu_lambda(s, min_cycles=5):
    """Smallest lambda (then smallest mu) with s[i]==s[i+lambda] for i>=mu
    and len(s)-mu >= min_cycles*lambda. Brute force, small inputs only."""
    n = len(s)
    for lam in range(1, n // min_cycles + 1):
        i = n - lam - 1
        while i >= 0 and s[i] == s[i + lam]:
            i -= 1
        mu = i + 1
        if n - mu >= min_cycles * lam:
            return mu, lam
    return None
