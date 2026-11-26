/*
 * Greedy duel sequence Du(k,0), J = 8, exact 256-bit signed integers.
 *
 * Usage: duel K CAP_PICKS [START_ROUNDS] [DUMP_PICKS]
 *
 * Every player has c[0] = number of picks taken, so the lex-min rule makes the
 * game proceed in rounds of k picks (each player once). Within a round only
 * the picker's vector changes, and it leaves the candidate set (c[0] + 1), so
 * the order of a round is exactly the sort of c[1..7] at the round start
 * (ties to lowest index). We sort once per round with insertion sort seeded by
 * the previous round's order.
 *
 * Periodicity is analysed at doubling checkpoints (in rounds). Output is a
 * single line "RESULT k=.. mu=.. lambda=.. picks=.. status=.. ties=.. first_tie=.."
 * plus "CHECK" lines per checkpoint on stderr.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define J 8
#define L 4 /* 4 x 64 = 256 bits */
#define MIN_TAIL_CYCLES 5 /* tail >= first cycle + 4 repeats */

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

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: duel K CAP_PICKS [START_ROUNDS] [DUMP_PICKS]\n"); return 2; }
    K = atoi(argv[1]);
    int64_t cap = (int64_t)atof(argv[2]);
    int64_t start_rounds = argc > 3 ? (int64_t)atof(argv[3]) : 1024;
    int64_t dump = argc > 4 ? (int64_t)atof(argv[4]) : 0;
    if (K < 2 || K > 255) { fprintf(stderr, "K out of range\n"); return 2; }

    int64_t cap_rounds = cap / K;
    int64_t N_cap = cap_rounds * K;
    uint8_t *s = malloc(N_cap);
    pi_buf = malloc(sizeof(int32_t) * (cap_rounds + 1));
    C = calloc(K, sizeof(*C));
    int *order = malloc(sizeof(int) * K);
    if (!s || !pi_buf || !C || !order) { fprintf(stderr, "out of memory\n"); return 3; }
    for (int i = 0; i < K; i++) order[i] = i;

    i256 B[J];
    memset(B, 0, sizeof B);
    B[0].w[0] = 1; /* C(0,0) */

    int64_t ties = 0, first_tie = -1, e = 0;
    int64_t next_check = start_rounds < cap_rounds ? start_rounds : cap_rounds;
    int have_prev = 0; int64_t prev_mu = 0, prev_lam = 0, prev_M = 0;
    int settled = 0; int64_t res_mu = -1, res_lam = -1, last_mu = -1, last_lam = -1;
    int overflow = 0;

    for (int64_t r = 0; r < cap_rounds; r++) {
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
        if (r > 0)
            for (int a = 1; a < K; a++)
                if (cmp_player(order[a - 1], order[a]) == 0) {
                    ties++;
                    if (first_tie < 0) first_tie = e + a - 1;
                }
        for (int pos = 0; pos < K; pos++, e++) {
            int p = order[pos];
            s[e] = (uint8_t)p;
            for (int j = 0; j < J; j++) {
                if (j & 1) sub256(&C[p][j], &B[j]); else add256(&C[p][j], &B[j]);
            }
            for (int j = J - 1; j >= 1; j--) add256(&B[j], &B[j - 1]);
        }
        if (near_overflow(&B[J - 1])) { overflow = 1; cap_rounds = r + 1; }

        if (r + 1 == next_check || r + 1 == cap_rounds) {
            int64_t M = r + 1, mu, lam;
            int found = analyse(s, M, &mu, &lam);
            fprintf(stderr, "CHECK k=%d picks=%lld found=%d mu=%lld lambda=%lld\n", K,
                    (long long)(M * K), found, found ? (long long)mu : -1LL, found ? (long long)lam : -1LL);
            fflush(stderr);
            if (found) { last_mu = mu; last_lam = lam; } else { last_mu = last_lam = -1; }
            /* confirmation must come from a true doubling, not a short final step at the cap */
            if (found && have_prev && mu == prev_mu && lam == prev_lam && M >= 2 * prev_M) {
                settled = 1; res_mu = mu; res_lam = lam; break;
            }
            have_prev = found; prev_mu = mu; prev_lam = lam; prev_M = M;
            next_check = 2 * M < cap_rounds ? 2 * M : cap_rounds;
        }
    }
    if (dump) {
        FILE *f = fopen("picks_dump.bin", "wb");
        fwrite(s, 1, dump < e ? dump : e, f);
        fclose(f);
    }
    printf("RESULT k=%d mu=%lld lambda=%lld picks=%lld status=%s ties=%lld first_tie=%lld cand_mu=%lld cand_lambda=%lld\n",
           K, (long long)res_mu, (long long)res_lam, (long long)e,
           ties ? "tie_J_too_small" : (settled ? "settled" : (overflow ? "overflow" : "unsettled")),
           (long long)ties, (long long)first_tie, (long long)last_mu, (long long)last_lam);
    return 0;
}
