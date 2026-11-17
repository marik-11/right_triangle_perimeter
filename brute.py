from math import isqrt


def triples_upto(n):
    total = 0
    for a in range(1, n // 3 + 1):
        for b in range(a + 1, (n - a) // 2 + 1):
            c2 = a * a + b * b
            c = isqrt(c2)
            if c * c == c2 and a + b + c <= n:
                total += 1
    return total


def mobius_sqfree_odd(limit):
    s = [True] * (limit + 1)
    s[0] = False
    for p in range(2, isqrt(limit) + 1):
        for k in range(p * p, limit + 1, p * p):
            s[k] = False
    return [s[i] and i % 2 == 1 for i in range(limit + 1)]


def region(m):
    total = 0
    a = 1
    while a * (a + 1) <= m:
        total += min(2 * a - 1, m // a) - a
        a += 1
    return total


def formula(n):
    x = n // 2
    odd_sqfree = mobius_sqfree_odd(max(x, 1))
    return sum(region(x // r) for r in range(1, x + 1) if odd_sqfree[r])


if __name__ == "__main__":
    for n in range(1, 300):
        assert triples_upto(n) == formula(n), n
    print([triples_upto(10 ** k) for k in range(1, 4)], [formula(10 ** k) for k in range(1, 5)])
