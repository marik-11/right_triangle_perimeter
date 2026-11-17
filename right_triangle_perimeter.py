import sys
from math import isqrt
import numpy as np

CHUNK = 1 << 23
DIRECT = 600


def sum_floor(m, lo, hi):
    total = 0
    for start in range(lo, hi + 1, CHUNK):
        stop = min(hi, start + CHUNK - 1)
        arr = np.arange(start, stop + 1, dtype=np.int64)
        total += int((np.int64(m) // arr).sum())
    return total


def feasible(m, x, y, q, p):
    return y - p > 0 and (x + q) * (y - p) > m


def hull_sum_root(m, a, b):
    if a > b:
        return 0
    x = a
    y = m // x + 1
    total = 0
    while x < b:
        room = b - x
        q1, p1, q2, p2 = 1, 0, 0, 1
        while q1 + q2 <= room:
            if feasible(m, x, y, q1 + q2, p1 + p2):
                k_lo, k_hi = 1, 2
                while q1 + k_hi * q2 <= room and feasible(m, x, y, q1 + k_hi * q2, p1 + k_hi * p2):
                    k_lo, k_hi = k_hi, k_hi << 1
                while k_hi - k_lo > 1:
                    mid = (k_lo + k_hi) >> 1
                    if q1 + mid * q2 <= room and feasible(m, x, y, q1 + mid * q2, p1 + mid * p2):
                        k_lo = mid
                    else:
                        k_hi = mid
                q1, p1 = q1 + k_lo * q2, p1 + k_lo * p2
            else:
                jmax = (room - q2) // q1
                if jmax < 1:
                    break

                def climbing(j):
                    xj = x + q2 + j * q1
                    return p1 == 0 or m * q1 > p1 * xj * (xj + q1)

                j_lo, j_hi = 0, 1
                while j_hi <= jmax and climbing(j_hi):
                    j_lo, j_hi = j_hi, j_hi << 1
                j_hi = min(j_hi, jmax + 1)
                while j_hi - j_lo > 1:
                    mid = (j_lo + j_hi) >> 1
                    if climbing(mid):
                        j_lo = mid
                    else:
                        j_hi = mid
                jc = min(j_hi, jmax)
                if jc < 1 or not feasible(m, x, y, q2 + jc * q1, p2 + jc * p1):
                    break
                s_lo, s_hi = 0, jc
                while s_hi - s_lo > 1:
                    mid = (s_lo + s_hi) >> 1
                    if feasible(m, x, y, q2 + mid * q1, p2 + mid * p1):
                        s_hi = mid
                    else:
                        s_lo = mid
                q2, p2 = q2 + s_lo * q1, p2 + s_lo * p1
        t_lo, t_hi = 0, 1
        while x + t_hi * q1 <= b and feasible(m, x, y, t_hi * q1, t_hi * p1):
            t_lo, t_hi = t_hi, t_hi << 1
        while t_hi - t_lo > 1:
            mid = (t_lo + t_hi) >> 1
            if x + mid * q1 <= b and feasible(m, x, y, mid * q1, mid * p1):
                t_lo = mid
            else:
                t_hi = mid
        t = t_lo
        total += t * q1 * (y - 1) - q1 * p1 * t * (t - 1) // 2 - t * (q1 - 1) * (p1 - 1) // 2
        x += t * q1
        y -= t * p1
    return total + m // b


def hull_sum(m, a, b):
    if a > b:
        return 0
    x = a
    y = m // x + 1
    total = 0
    runs = []
    q1, p1, q2, p2 = 1, 0, 0, 1
    while x < b:
        room = b - x
        while q1 + q2 <= room:
            if feasible(m, x, y, q1 + q2, p1 + p2):
                k_lo, k_hi = 1, 2
                while q1 + k_hi * q2 <= room and feasible(m, x, y, q1 + k_hi * q2, p1 + k_hi * p2):
                    k_lo, k_hi = k_hi, k_hi << 1
                while k_hi - k_lo > 1:
                    mid = (k_lo + k_hi) >> 1
                    if q1 + mid * q2 <= room and feasible(m, x, y, q1 + mid * q2, p1 + mid * p2):
                        k_lo = mid
                    else:
                        k_hi = mid
                runs.append([q2, p2, k_lo])
                q1, p1 = q1 + k_lo * q2, p1 + k_lo * p2
            else:
                jmax = (room - q2) // q1
                if jmax < 1:
                    break

                def climbing(j):
                    xj = x + q2 + j * q1
                    return p1 == 0 or m * q1 > p1 * xj * (xj + q1)

                j_lo, j_hi = 0, 1
                while j_hi <= jmax and climbing(j_hi):
                    j_lo, j_hi = j_hi, j_hi << 1
                j_hi = min(j_hi, jmax + 1)
                while j_hi - j_lo > 1:
                    mid = (j_lo + j_hi) >> 1
                    if climbing(mid):
                        j_lo = mid
                    else:
                        j_hi = mid
                jc = min(j_hi, jmax)
                if jc < 1 or not feasible(m, x, y, q2 + jc * q1, p2 + jc * p1):
                    break
                s_lo, s_hi = 0, jc
                while s_hi - s_lo > 1:
                    mid = (s_lo + s_hi) >> 1
                    if feasible(m, x, y, q2 + mid * q1, p2 + mid * p1):
                        s_hi = mid
                    else:
                        s_lo = mid
                q2, p2 = q2 + s_lo * q1, p2 + s_lo * p1
        t_lo, t_hi = 0, 1
        while x + t_hi * q1 <= b and feasible(m, x, y, t_hi * q1, t_hi * p1):
            t_lo, t_hi = t_hi, t_hi << 1
        while t_hi - t_lo > 1:
            mid = (t_lo + t_hi) >> 1
            if x + mid * q1 <= b and feasible(m, x, y, mid * q1, mid * p1):
                t_lo = mid
            else:
                t_hi = mid
        t = t_lo
        total += t * q1 * (y - 1) - q1 * p1 * t * (t - 1) // 2 - t * (q1 - 1) * (p1 - 1) // 2
        x += t * q1
        y -= t * p1
        q2, p2 = q1, p1
        room = b - x
        if not runs:
            q1, p1, q2, p2 = 1, 0, 0, 1
            continue
        while True:
            if not runs:
                q1, p1 = 1, 0
                break
            hq, hp, k = runs[-1]
            q1 -= hq
            p1 -= hp
            if k == 1:
                runs.pop()
            else:
                runs[-1][2] = k - 1
            if q1 <= room and feasible(m, x, y, q1, p1):
                break
            q2, p2 = q1, p1
    return total + m // b


def region(m, hull=True):
    if m < 6:
        return 0
    a1 = (1 + isqrt(1 + 8 * m)) // 4
    a2 = (isqrt(1 + 4 * m) - 1) // 2
    total = a1 * (a1 - 1) // 2
    if a2 > a1:
        if hull and a2 - a1 >= DIRECT:
            s = hull_sum(m, a1 + 1, a2)
        else:
            s = sum_floor(m, a1 + 1, a2)
        total += s - (a2 * (a2 + 1) - a1 * (a1 + 1)) // 2
    return total


def mobius(limit):
    mu = np.ones(limit + 1, dtype=np.int8)
    comp = np.zeros(limit + 1, dtype=bool)
    for p in range(2, limit + 1):
        if comp[p]:
            continue
        comp[2 * p::p] = True
        mu[p::p] *= -1
        if p * p <= limit:
            mu[p * p::p * p] = 0
    return mu


class Squarefree:
    def __init__(self, dmax):
        mu = mobius(dmax)
        odd = np.arange(dmax + 1, dtype=np.int64) & 1
        self.prefix = np.cumsum(mu.astype(np.int64) * odd)
        d = np.arange(1, dmax + 1, 2, dtype=np.int64)
        d = d[mu[d] != 0]
        self.d = d
        self.sq = d * d
        self.sign = mu[d].astype(np.int64)

    def count_basic(self, z):
        cnt = int(np.searchsorted(self.sq, z, side="right"))
        return int(np.dot(self.sign[:cnt], (z // self.sq[:cnt] + 1) // 2))

    def count(self, z):
        dm = isqrt(z)
        d1 = min(dm, int(2.5 * z ** (1.0 / 3.0)) + 1)
        cnt = int(np.searchsorted(self.d, d1, side="right"))
        q = int(np.dot(self.sign[:cnt], (z // self.sq[:cnt] + 1) // 2))
        if dm > d1:
            top = z // ((d1 + 1) ** 2)
            v = np.arange(1, top + 2, dtype=np.int64)
            qq = z // v
            e = np.sqrt(qq.astype(np.float64)).astype(np.int64)
            e = np.where(e * e > qq, e - 1, e)
            e = np.where((e + 1) * (e + 1) <= qq, e + 1, e)
            e = np.maximum(np.minimum(e, dm), d1)
            q += int(np.dot((v[:-1] + 1) // 2, self.prefix[e[:-1]] - self.prefix[e[1:]]))
        return q


def triples(n, ratio=1.5, method="hull"):
    x = n // 2
    if x < 6:
        return 0
    sq = isqrt(x)
    y = max(1, int(ratio * sq))
    rmax = x // (y + 1)
    sf = Squarefree(sq + 2)
    g = np.zeros(y + 1, dtype=np.int64)
    a = 1
    while a * (a + 1) <= y:
        hi = min(2 * a - 1, y // a)
        if hi > a:
            g[a * np.arange(a + 1, hi + 1, dtype=np.int64)] += 1
        a += 1
    count = sf.count if method == "hull" else sf.count_basic
    total = 0
    for k in np.nonzero(g)[0].tolist():
        total += count(x // k) * int(g[k])
    s = (np.arange(rmax + 1) & 1).astype(bool)
    s[0] = False
    p = 3
    while p * p <= rmax:
        s[p * p::p * p] = False
        p += 2
    ry = region(y, False)
    for r in np.nonzero(s)[0].tolist():
        total += region(x // r, method == "hull") - ry
    return total


if __name__ == "__main__":
    print(triples(int(sys.argv[1])))
