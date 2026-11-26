/*
 * Greedy duel sequence Du(k,0), J = 8, exact signed integers of 64*L bits
 * (L = 4 by default; build with -DL=5 for 320 bits, see run.py's duel320).
 *
 * Usage: duel K CAP_PICKS [START_ROUNDS] [DUMP_PICKS] [STORE_PICKS] [STATE_AT]
 *             [MIN_REPEATS] [RESUME]
 *
 * CAP_PICKS is the absolute pick to stop at. STORE_PICKS is how many picks are kept
 * for analysis (about 1.1 bytes each), counted from the storage start.
 *
 * MIN_REPEATS (default 1000) is the repeat floor of the built-in settle rule; a large
 * value keeps verifying a candidate up to CAP_PICKS. 0 stops at the first candidate
 * (status=candidate): run.py then settles it with certify.py.
 *
 * RESUME: jump-ahead start. RESUME.txt holds "X P", then K lines of c[0..7] and one
 * line of C(X, 0..7), all hex two's complement with 16*L digits; RESUME.bin holds the
 * P picks before X. Play starts at pick X; storage (and mu) starts at X - P.
 *
 * STATE_AT (a multiple of K): write every player's c[0..7] at that pick to
 * state_dump.txt (hex, two's complement), and make picks_dump.bin hold picks
 * [STATE_AT, STATE_AT + DUMP_PICKS) instead of [0, DUMP_PICKS). Used by certify.py.
 *
 * Every player has c[0] = number of picks taken, so the lex-min rule makes the
 * game proceed in rounds of k picks (each player once). Within a round only
 * the picker's vector changes, and it leaves the candidate set (c[0] + 1), so
 * the order of a round is exactly the sort of c[1..7] at the round start
 * (ties to lowest index). We sort once per round with insertion sort seeded by
 * the previous round's order.
 *
 * Periodicity is analysed at doubling checkpoints (in rounds). A candidate
 * (mu, lambda) is then verified pick by pick against its first cycle; it is
 * settled only when the tail after mu reaches max(10 (mu + lambda), 1000 lambda)
 * picks. A mismatch prints BREAK and analysis resumes at the next checkpoint.
 * Verification needs no extra memory, so it continues past STORE_PICKS (the
 * stored prefix, about 1.1 bytes per pick) up to CAP_PICKS.
 * Output is a single line "RESULT k=.. mu=.. lambda=.. picks=.. status=.. ...
 * repeats=.. break=.." plus CHECK / VERIFY / BREAK lines on stderr.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define J 8
#ifndef L
#define L 4 /* 4 x 64 = 256 bits */
#endif
#define MIN_TAIL_CYCLES 5 /* candidate detection: tail >= first cycle + 4 repeats */
#define MIN_REPEATS 1000 /* settled: tail >= 1000 lambda ... */
#define MIN_TAIL_FACTOR 10 /* ... and tail >= 10 (mu + lambda) */

typedef struct { uint64_t w[L]; } i256;

static inline void add256(i256 *a, const i256 *b) {
    unsigned __int128 carry = 0;
    for (int i = 0; i < L; i++) {
        carry += (unsigned __int128)a->w[i] + b->w[i];
        a->w[i] = (uint64_t)carry;
        carry >>= 64;
    }
}
static inline void sub256(i256 *a, const i256 *b) {
    uint64_t borrow = 0;
    for (int i = 0; i < L; i++) {
        uint64_t x = a->w[i], y = b->w[i];
        uint64_t d = x - y - borrow;
        borrow = (x < y) || (x - y < borrow);
        a->w[i] = d;
    }
}
static inline int cmp256(const i256 *a, const i256 *b) {
    int64_t ha = (int64_t)a->w[L - 1], hb = (int64_t)b->w[L - 1];
    if (ha != hb) return ha < hb ? -1 : 1;
    for (int i = L - 2; i >= 0; i--)
        if (a->w[i] != b->w[i]) return a->w[i] < b->w[i] ? -1 : 1;
    return 0;
}
/* overflow guard: top 16 bits must be pure sign extension */
static inline int near_overflow(const i256 *a) {
    uint64_t top = a->w[L - 1] >> 48;
    return !(top == 0 || top == 0xFFFF);
}
/* accumulator guard: |a| < 2^253 (top 3 bits equal). c[J-1] of a player grows like
 * e^8 / (8! k), far faster than C(e, J-1), so it is the one that can wrap. */
static inline int acc_near_overflow(const i256 *a) {
    uint64_t top = a->w[L - 1] >> 61;
    return !(top == 0 || top == 7);
}
/* bit length of |a| */
static int bits256(const i256 *a) {
    i256 x = *a;
    if ((int64_t)x.w[L - 1] < 0) { i256 z; memset(&z, 0, sizeof z); sub256(&z, &x); x = z; }
    for (int i = L - 1; i >= 0; i--)
        if (x.w[i]) return 64 * i + 64 - __builtin_clzll(x.w[i]);
    return 0;
}

static int K;
static i256 (*C)[J];

/* lex compare c[1..J-1] (c[0] equal for all players at round start) */
static inline int cmp_player(int a, int b) {
    for (int j = 1; j < J; j++) {
        int r = cmp256(&C[a][j], &C[b][j]);
        if (r) return r;
    }
    return 0;
}

/* ---------- periodicity analysis ---------- */
static int32_t *pi_buf;

/* returns 1 and fills mu, lam if a periodic tail with >= MIN_TAIL_CYCLES is found */
static int analyse(const uint8_t *s, int64_t M, int64_t *mu_out, int64_t *lam_out) {
    int64_t N = M * K;
    /* prefix function over reversed rounds; round i of reversed = round M-1-i */
#define RND(i) (s + (M - 1 - (i)) * (int64_t)K)
    pi_buf[0] = 0;
    for (int64_t i = 1; i < M; i++) {
        int32_t q = pi_buf[i - 1];
        while (q > 0 && memcmp(RND(i), RND(q), K) != 0) q = pi_buf[q - 1];
        if (memcmp(RND(i), RND(q), K) == 0) q++;
        pi_buf[i] = q;
    }
    /* per(Lr) = Lr - pi[Lr-1] is non-decreasing in Lr */
    int64_t Lstar = 0;
    for (int64_t Lr = M; Lr >= 1; Lr--) {
        int64_t per = Lr - pi_buf[Lr - 1];
        if (MIN_TAIL_CYCLES * per <= Lr) { Lstar = Lr; break; }
    }
    if (!Lstar) return 0;
    int64_t Lam = Lstar - pi_buf[Lstar - 1];
    int64_t Lmax = Lstar;
    while (Lmax < M && (Lmax + 1) - pi_buf[Lmax] == Lam) Lmax++;
#undef RND
    /* pick level: smallest divisor d of K*Lam that is a period on the tail */
    int64_t P = (int64_t)K * Lam, start = N - Lmax * K, lam = P;
    for (int64_t d = 1; d < P; d++) {
        if (P % d) continue;
        int ok = 1;
        for (int64_t i = start; i + d < N; i++)
            if (s[i] != s[i + d]) { ok = 0; break; }
        if (ok) { lam = d; break; }
    }
    int64_t i = N - lam - 1;
    while (i >= 0 && s[i] == s[i + lam]) i--;
    int64_t mu = i + 1;
    if (N - mu < MIN_TAIL_CYCLES * lam) return 0;
    *mu_out = mu; *lam_out = lam;
    return 1;
}

static void print_hex(FILE *f, const i256 *a) {
    for (int i = L - 1; i >= 0; i--) fprintf(f, "%016llx", (unsigned long long)a->w[i]);
}
static int read_hex(FILE *f, i256 *a) {
    char buf[16 * L + 1];
    if (fscanf(f, "%s", buf) != 1 || strlen(buf) != 16 * L) return 0;
    for (int i = 0; i < L; i++) {
        char part[17];
        memcpy(part, buf + 16 * (L - 1 - i), 16); part[16] = 0;
        a->w[i] = strtoull(part, NULL, 16);
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: duel K CAP_PICKS [START_ROUNDS] [DUMP_PICKS] [STORE_PICKS] [STATE_AT] [MIN_REPEATS] [RESUME]\n"); return 2; }
    K = atoi(argv[1]);
    int64_t cap = (int64_t)atof(argv[2]);
    int64_t start_rounds = argc > 3 ? (int64_t)atof(argv[3]) : 1024;
    int64_t dump = argc > 4 ? (int64_t)atof(argv[4]) : 0;
    int64_t store = argc > 5 ? (int64_t)atof(argv[5]) : cap;
    int64_t state_at = argc > 6 ? (int64_t)atof(argv[6]) : -1;
    int64_t min_repeats = argc > 7 ? (int64_t)atof(argv[7]) : MIN_REPEATS;
    const char *resume = argc > 8 ? argv[8] : NULL;
    if (K < 2 || K > 255) { fprintf(stderr, "K out of range\n"); return 2; }

    C = calloc(K, sizeof(*C));
    int *order = malloc(sizeof(int) * K);
    i256 B[J];
    memset(B, 0, sizeof B);
    B[0].w[0] = 1; /* C(0,0) */
    int64_t X = 0, P = 0;
    FILE *rf = NULL;
    if (resume) {
        char path[4096];
        snprintf(path, sizeof path, "%s.txt", resume);
        rf = fopen(path, "r");
        long long x_, p_;
        if (!rf || fscanf(rf, "%lld %lld", &x_, &p_) != 2) { fprintf(stderr, "bad resume file\n"); return 2; }
        X = x_; P = p_;
        if (X % K || P % K || P > X) { fprintf(stderr, "resume X, P must be multiples of K, P <= X\n"); return 2; }
        for (int i = 0; i < K; i++)
            for (int j = 0; j < J; j++)
                if (!read_hex(rf, &C[i][j])) { fprintf(stderr, "bad resume state\n"); return 2; }
        for (int j = 0; j < J; j++)
            if (!read_hex(rf, &B[j])) { fprintf(stderr, "bad resume binomials\n"); return 2; }
        fclose(rf);
    }
    const int64_t sb = X - P; /* absolute pick of storage index 0 */
    if (store <= 0 || store > cap - sb) store = cap - sb;
    if (store < P) store = P;
    int64_t store_rounds = store / K;
    int64_t N_store = store_rounds * K;
    uint8_t *s = malloc(N_store > 0 ? N_store : 1);
    pi_buf = malloc(sizeof(int32_t) * (store_rounds + 1));
    if (!s || !pi_buf || !C || !order) { fprintf(stderr, "out of memory\n"); return 3; }
    if (P) {
        char path[4096];
        snprintf(path, sizeof path, "%s.bin", resume);
        FILE *pf = fopen(path, "rb");
        if (!pf || (int64_t)fread(s, 1, P, pf) != P) { fprintf(stderr, "bad resume prefix\n"); return 2; }
        fclose(pf);
    }
    for (int i = 0; i < K; i++) order[i] = i;

    int64_t ties = 0, first_tie = -1, e = X;
    int64_t next_check = start_rounds;
    while (next_check <= P / K) next_check *= 2;
    if (next_check > store_rounds) next_check = store_rounds;
    int settled = 0, overflow = 0;
    /* candidate under verification: every pick e must equal s[v_mu - sb + (e - v_mu) % v_lam] */
    int verifying = 0;
    int64_t v_mu = -1, v_lam = -1, v_target = 0, v_idx = 0, v_end = 0, brk = -1;

    while (e + K <= cap) {
        if (e == state_at) {
            FILE *f = fopen("state_dump.txt", "w");
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < J; j++) { if (j) fputc(' ', f); print_hex(f, &C[i][j]); }
                fprintf(f, "\n");
            }
            fclose(f);
        }
        /* insertion sort of order[] by c[1..], ties to lowest index */
        for (int a = 1; a < K; a++) {
            int x = order[a], b = a - 1;
            while (b >= 0) {
                int cr = cmp_player(order[b], x);
                if (cr > 0 || (cr == 0 && order[b] > x)) { order[b + 1] = order[b]; b--; }
                else break;
            }
            order[b + 1] = x;
        }
        if (e > 0)
            for (int a = 1; a < K; a++)
                if (cmp_player(order[a - 1], order[a]) == 0) {
                    ties++;
                    if (first_tie < 0) first_tie = e + a - 1;
                }
        for (int pos = 0; pos < K; pos++, e++) {
            int p = order[pos];
            if (e - sb < N_store) s[e - sb] = (uint8_t)p;
            if (verifying) {
                if ((uint8_t)p != s[v_mu - sb + v_idx]) {
                    verifying = 0; brk = e; v_end = e;
                    fprintf(stderr, "BREAK k=%d mu=%lld lambda=%lld at pick %lld after %lld repeats\n", K,
                            (long long)v_mu, (long long)v_lam, (long long)e, (long long)((e - v_mu) / v_lam));
                    fflush(stderr);
                } else if (++v_idx == v_lam) v_idx = 0;
            }
            for (int j = 0; j < J; j++) {
                if (j & 1) sub256(&C[p][j], &B[j]); else add256(&C[p][j], &B[j]);
            }
            if (acc_near_overflow(&C[p][J - 1])) overflow = 1;
            for (int j = J - 1; j >= 1; j--) add256(&B[j], &B[j - 1]);
        }
        if (near_overflow(&B[J - 1])) overflow = 1;
        if (overflow) { settled = 0; break; }
        if (verifying) v_end = e;
        if (verifying && e >= v_target) { settled = 1; break; }

        int64_t M = (e - sb) / K; /* stored rounds so far */
        if ((M == next_check || M == store_rounds) && M <= store_rounds) {
            if (verifying) {
                fprintf(stderr, "VERIFY k=%d picks=%lld mu=%lld lambda=%lld repeats=%lld target=%lld\n", K,
                        (long long)e, (long long)v_mu, (long long)v_lam,
                        (long long)((e - v_mu) / v_lam), (long long)v_target);
            } else {
                int64_t mu, lam;
                int found = analyse(s, M, &mu, &lam);
                if (found) mu += sb;
                fprintf(stderr, "CHECK k=%d picks=%lld found=%d mu=%lld lambda=%lld\n", K,
                        (long long)e, found, found ? (long long)mu : -1LL, found ? (long long)lam : -1LL);
                if (found) {
                    int64_t a = MIN_TAIL_FACTOR * (mu + lam), b = min_repeats * lam;
                    v_mu = mu; v_lam = lam; v_target = min_repeats ? mu + (a > b ? a : b) : e;
                    v_idx = (e - mu) % lam; v_end = e; verifying = 1; brk = -1; /* break= refers to this candidate */
                    if (e >= v_target) { settled = 1; break; }
                }
            }
            fflush(stderr);
            next_check = 2 * M < store_rounds ? 2 * M : store_rounds;
        }
        /* past storage, only a live candidate can make progress */
        if (M >= store_rounds && !verifying) break;
    }
    if (dump) {
        FILE *f = fopen("picks_dump.bin", "wb");
        int64_t from = state_at > 0 ? state_at : sb;
        int64_t to = from + dump;
        if (to > e) to = e;
        if (to > sb + N_store) to = sb + N_store;
        if (from >= sb && to > from) fwrite(s + (from - sb), 1, to - from, f);
        fclose(f);
    }
    int64_t reps = v_lam > 0 ? (v_end - v_mu) / v_lam : 0;
    int maxbits = 0;
    for (int i = 0; i < K; i++)
        for (int j = 0; j < J; j++) { int b = bits256(&C[i][j]); if (b > maxbits) maxbits = b; }
    const char *status = ties ? "tie_J_too_small"
                       : settled ? (min_repeats ? "settled" : "candidate")
                       : overflow ? "overflow" : "unsettled";
    printf("RESULT k=%d mu=%lld lambda=%lld picks=%lld status=%s ties=%lld first_tie=%lld cand_mu=%lld cand_lambda=%lld repeats=%lld break=%lld bits=%d width=%d\n",
           K, settled ? (long long)v_mu : -1LL, settled ? (long long)v_lam : -1LL, (long long)e, status,
           (long long)ties, (long long)first_tie, (long long)v_mu, (long long)v_lam,
           (long long)reps, (long long)brk, maxbits, 64 * L);
    return 0;
}
