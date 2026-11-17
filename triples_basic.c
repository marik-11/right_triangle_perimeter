#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef unsigned long long u64;
typedef __int128 i128;

static u64 isqrt64(u64 n) {
    u64 x = (u64)sqrtl((long double)n);
    while (x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

static u64 region(u64 m) {
    if (m < 6) return 0;
    u64 a1 = (1 + isqrt64(1 + 8 * m)) / 4;
    while (a1 * (2 * a1 - 1) > m) a1--;
    while ((a1 + 1) * (2 * a1 + 1) <= m) a1++;
    u64 a2 = (isqrt64(1 + 4 * m) - 1) / 2;
    while (a2 * (a2 + 1) > m) a2--;
    while ((a2 + 1) * (a2 + 2) <= m) a2++;
    u64 total = a1 * (a1 - 1) / 2;
    if (a2 > a1) {
        u64 s = 0;
        for (u64 a = a1 + 1; a <= a2; a++) s += m / a;
        total += s - (a2 * (a2 + 1) - a1 * (a1 + 1)) / 2;
    }
    return total;
}

int main(int argc, char **argv) {
    u64 n = strtoull(argv[1], NULL, 10);
    double ratio = argc > 2 ? atof(argv[2]) : 0.3;
    int part = argc > 3 ? atoi(argv[3]) : 0;
    u64 lo_arg = argc > 4 ? strtoull(argv[4], NULL, 10) : 1;
    u64 hi_arg = argc > 5 ? strtoull(argv[5], NULL, 10) : ~0ULL;
    u64 x = n / 2;
    if (x < 6) { printf("0\n"); return 0; }
    u64 sq = isqrt64(x);
    u64 y = (u64)(ratio * (double)sq);
    if (y < 1) y = 1;
    u64 rmax = x / (y + 1);
    u64 dmax = sq + 2;
    signed char *mu = malloc(dmax + 1);
    for (u64 i = 0; i <= dmax; i++) mu[i] = 1;
    char *comp = calloc(dmax + 1, 1);
    for (u64 p = 2; p <= dmax; p++) {
        if (comp[p]) continue;
        for (u64 k = p; k <= dmax; k += p) { if (k > p) comp[k] = 1; mu[k] = -mu[k]; }
        if (p * p <= dmax) for (u64 k = p * p; k <= dmax; k += p * p) mu[k] = 0;
    }
    u64 *dl = malloc(sizeof(u64) * (dmax + 1));
    u64 dc = 0;
    for (u64 d = 1; d <= dmax; d += 2) if (mu[d]) dl[dc++] = d;
    unsigned char *s = malloc(rmax + 1);
    for (u64 i = 0; i <= rmax; i++) s[i] = (i & 1);
    s[0] = 0;
    for (u64 p = 3; p * p <= rmax; p += 2)
        for (u64 k = p * p; k <= rmax; k += p * p) s[k] = 0;
    unsigned int *g = calloc(y + 1, sizeof(unsigned int));
    for (u64 a = 1; a * (a + 1) <= y; a++)
        for (u64 b = a + 1; b < 2 * a && a * b <= y; b++) g[a * b]++;
    i128 total = 0;
    u64 k_lo = part == 2 ? 1 : lo_arg, k_hi = part == 2 ? 0 : (hi_arg < y ? hi_arg : y);
    if (part == 0) { k_lo = 1; k_hi = y; }
    for (u64 k = k_lo; k <= k_hi; k++) {
        if (!g[k]) continue;
        u64 z = x / k;
        i128 q = 0;
        for (u64 i = 0; i < dc; i++) {
            u64 d = dl[i];
            u64 dd = d * d;
            if (dd > z) break;
            u64 odd = (z / dd + 1) / 2;
            q += mu[d] > 0 ? (i128)odd : -(i128)odd;
        }
        total += q * g[k];
    }
    u64 ry = region(y);
    if (part == 3) { printf("%llu %llu\n", y, rmax); return 0; }
    u64 r_lo = part == 1 ? 1 : (part == 2 ? lo_arg : 1), r_hi = part == 1 ? 0 : (part == 2 ? (hi_arg < rmax ? hi_arg : rmax) : rmax);
    if (part == 0) { r_lo = 1; r_hi = rmax; }
    for (u64 r = r_lo; r <= r_hi; r++) {
        if (!s[r]) continue;
        total += (i128)region(x / r) - (i128)ry;
    }
    char buf[64];
    int len = 0;
    u64 hi_part = (u64)(total / 10000000000000000000ULL);
    u64 lo_part = (u64)(total % 10000000000000000000ULL);
    if (hi_part) printf("%llu%019llu\n", hi_part, lo_part); else printf("%llu\n", lo_part);
    (void)buf; (void)len;
    return 0;
}
