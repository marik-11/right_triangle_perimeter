#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef long long ll;
typedef unsigned long long u64;
typedef __int128 i128;
#ifndef DIRECT_THR
#define DIRECT_THR 600
#endif
#ifndef D1_MUL
#define D1_MUL 2.5
#endif

static u64 isqrt64(u64 n) {
    u64 x = (u64)sqrt((double)n);
    while (x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

static ll inv_mod(ll b, ll a) {
    ll r0 = a, r1 = b % a, s0 = 0, s1 = 1;
    while (r1) {
        ll q = r0 / r1, t = r0 - q * r1;
        r0 = r1; r1 = t;
        t = s0 - q * s1;
        s0 = s1; s1 = t;
    }
    s0 %= a;
    if (s0 < 0) s0 += a;
    return s0;
}

static inline int feasible(ll m, ll x, ll y, ll q, ll p) {
    ll yy = y - p;
    if (yy <= 0) return 0;
    return (u64)(x + q) * (u64)yy > (u64)m;
}

static ll hull_sum(ll m, ll a, ll b) {
    if (a > b) return 0;
    ll x = a, y = m / x + 1;
    i128 total = 0;
    ll q1 = 1, p1 = 0, q2 = 0, p2 = 1;
    ll sq_[256], sp_[256], sk_[256];
    int sd = 0;
    while (x < b) {
        ll room = b - x;
        for (;;) {
            if (q1 + q2 > room) break;
            if (feasible(m, x, y, q1 + q2, p1 + p2)) {
                ll klo = 1, khi = 2;
                for (;;) {
                    ll qq = q1 + khi * q2;
                    if (qq > room || !feasible(m, x, y, qq, p1 + khi * p2)) break;
                    klo = khi; khi <<= 1;
                }
                while (khi - klo > 1) {
                    ll mid = (klo + khi) >> 1;
                    ll qq = q1 + mid * q2;
                    if (qq <= room && feasible(m, x, y, qq, p1 + mid * p2)) klo = mid; else khi = mid;
                }
                sq_[sd] = q2; sp_[sd] = p2; sk_[sd] = klo; sd++;
                q1 += klo * q2; p1 += klo * p2;
            } else {
                ll jmax = q1 ? (room - q2) / q1 : 0;
                if (jmax < 1) break;
                ll jlo = 0, jhi = 1;
                for (;;) {
                    if (jhi > jmax) break;
                    ll xj = x + q2 + jhi * q1;
                    int climb = p1 == 0 || (i128)m * q1 > (i128)p1 * xj * (xj + q1);
                    if (!climb) break;
                    jlo = jhi; jhi <<= 1;
                }
                if (jhi > jmax + 1) jhi = jmax + 1;
                while (jhi - jlo > 1) {
                    ll mid = (jlo + jhi) >> 1;
                    ll xj = x + q2 + mid * q1;
                    int climb = p1 == 0 || (i128)m * q1 > (i128)p1 * xj * (xj + q1);
                    if (climb) jlo = mid; else jhi = mid;
                }
                ll jc = jhi < jmax ? jhi : jmax;
                if (jc < 1 || !feasible(m, x, y, q2 + jc * q1, p2 + jc * p1)) break;
                ll slo = 0, shi = jc;
                while (shi - slo > 1) {
                    ll mid = (slo + shi) >> 1;
                    if (feasible(m, x, y, q2 + mid * q1, p2 + mid * p1)) shi = mid; else slo = mid;
                }
                q2 += slo * q1; p2 += slo * p1;
            }
        }
        ll da = q1, db = p1;
        ll tlo = 0, thi = 1;
        for (;;) {
            ll nx = x + thi * da;
            if (nx > b || !feasible(m, x, y, thi * da, thi * db)) break;
            tlo = thi; thi <<= 1;
        }
        while (thi - tlo > 1) {
            ll mid = (tlo + thi) >> 1;
            ll nx = x + mid * da;
            if (nx <= b && feasible(m, x, y, mid * da, mid * db)) tlo = mid; else thi = mid;
        }
        ll t = tlo;
        total += (i128)t * da * (y - 1) - (i128)da * db * t * (t - 1) / 2 - (i128)t * (da - 1) * (db - 1) / 2;
        x += t * da; y -= t * db;
        q2 = da; p2 = db;
        room = b - x;
        if (sd == 0) { q1 = 1; p1 = 0; q2 = 0; p2 = 1; continue; }
        for (;;) {
            if (sd == 0) { q1 = 1; p1 = 0; break; }
            q1 -= sq_[sd - 1]; p1 -= sp_[sd - 1];
            if (--sk_[sd - 1] == 0) sd--;
            if (q1 <= room && feasible(m, x, y, q1, p1)) break;
            q2 = q1; p2 = p1;
        }
    }
    return (ll)(total + m / b);
}

static ll region(ll m) {
    if (m < 6) return 0;
    ll a1 = (1 + (ll)isqrt64(1 + 8 * (u64)m)) / 4;
    ll a2 = ((ll)isqrt64(1 + 4 * (u64)m) - 1) / 2;
    i128 total = (i128)a1 * (a1 - 1) / 2;
    if (a2 > a1) {
        ll s = 0;
        if (a2 - a1 < DIRECT_THR) { for (ll a = a1 + 1; a <= a2; a++) s += m / a; }
        else s = hull_sum(m, a1 + 1, a2);
        total += (i128)s - ((i128)a2 * (a2 + 1) - (i128)a1 * (a1 + 1)) / 2;
    }
    return (ll)total;
}

static int *pm;
static signed char *mu;
static int *dl;
static ll dc;

static ll oddsqfree(ll z) {
    ll d_max = (ll)isqrt64((u64)z);
    ll d1 = (ll)(D1_MUL * cbrt((double)z)) + 1;
    if (d1 > d_max) d1 = d_max;
    ll q = 0;
    for (ll i = 0; i < dc; i++) {
        ll d = dl[i];
        if (d > d1) break;
        ll odd = (z / (d * d) + 1) / 2;
        q += mu[d] > 0 ? odd : -odd;
    }
    ll ds = d1 + 1;
    while (ds <= d_max) {
        ll v = z / (ds * ds);
        ll de = (ll)isqrt64((u64)(z / v));
        if (de > d_max) de = d_max;
        q += ((v + 1) / 2) * (ll)(pm[de] - pm[ds - 1]);
        ds = de + 1;
    }
    return q;
}

int main(int argc, char **argv) {
    u64 n = strtoull(argv[1], NULL, 10);
    double ratio = argc > 2 ? atof(argv[2]) : 1.0;
    int part = argc > 3 ? atoi(argv[3]) : 0;
    u64 lo_arg = argc > 4 ? strtoull(argv[4], NULL, 10) : 1;
    u64 hi_arg = argc > 5 ? strtoull(argv[5], NULL, 10) : ~0ULL;
    u64 x = n / 2;
    if (x > 3000000000000000000ULL) { fprintf(stderr, "N too large\n"); return 1; }
    if (x < 6) { printf("0\n"); return 0; }
    u64 sq = isqrt64(x);
    u64 y = (u64)(ratio * (double)sq);
    if (y < 1) y = 1;
    u64 rmax = x / (y + 1);
    u64 dmax = sq + 2;
    mu = malloc(dmax + 1);
    for (u64 i = 0; i <= dmax; i++) mu[i] = 1;
    char *comp = calloc(dmax + 1, 1);
    for (u64 p = 2; p <= dmax; p++) {
        if (comp[p]) continue;
        for (u64 k = p; k <= dmax; k += p) { if (k > p) comp[k] = 1; mu[k] = -mu[k]; }
        if (p * p <= dmax) for (u64 k = p * p; k <= dmax; k += p * p) mu[k] = 0;
    }
    free(comp);
    pm = malloc(sizeof(int) * (dmax + 1));
    pm[0] = 0;
    dl = malloc(sizeof(int) * (dmax / 2 + 2) * 1 + 8 * (dmax / 8 + 1));
    dc = 0;
    for (u64 d = 1; d <= dmax; d++) {
        pm[d] = pm[d - 1] + ((d & 1) ? mu[d] : 0);
        if ((d & 1) && mu[d]) dl[dc++] = (int)d;
    }
    if (part == 3) { printf("%llu %llu\n", y, rmax); return 0; }
    unsigned char *s = NULL;
    unsigned short *g = NULL;
    if (part != 1) {
        s = malloc(rmax + 1);
        for (u64 i = 0; i <= rmax; i++) s[i] = (i & 1);
        s[0] = 0;
        for (u64 p = 3; p * p <= rmax; p += 2)
            for (u64 k = p * p; k <= rmax; k += p * p) s[k] = 0;
    }
    if (part != 2) {
        g = calloc(y + 1, sizeof(unsigned short));
        for (u64 a = 1; a * (a + 1) <= y; a++)
            for (u64 b = a + 1; b < 2 * a && a * b <= y; b++) g[a * b]++;
    }
    i128 total = 0;
    u64 k_lo = part == 2 ? 1 : lo_arg, k_hi = part == 2 ? 0 : (hi_arg < y ? hi_arg : y);
    if (part == 0) { k_lo = 1; k_hi = y; }
    for (u64 k = k_lo; k <= k_hi; k++) {
        if (!g[k]) continue;
        total += (i128)oddsqfree((ll)(x / k)) * g[k];
    }
    ll ry = region((ll)y);
    u64 r_lo = part == 2 ? lo_arg : 1, r_hi = part == 1 ? 0 : (part == 2 ? (hi_arg < rmax ? hi_arg : rmax) : rmax);
    if (part == 0) { r_lo = 1; r_hi = rmax; }
    for (u64 r = r_lo; r <= r_hi; r++) {
        if (!s[r]) continue;
        total += (i128)region((ll)(x / r)) - (i128)ry;
    }
    u64 hi_part = (u64)(total / 10000000000000000000ULL);
    u64 lo_part = (u64)(total % 10000000000000000000ULL);
    if (hi_part) printf("%llu%019llu\n", hi_part, lo_part); else printf("%llu\n", lo_part);
    return 0;
}
