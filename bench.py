import sys
import time
import subprocess
from math import gcd
from right_triangle_perimeter import triples


def enumeration(limit):
    total = 0
    m = 2
    while 2 * m * (m + 1) <= limit:
        for n in range(1 + m % 2, m, 2):
            p = 2 * m * (m + n)
            if p > limit:
                break
            if gcd(m, n) == 1:
                total += limit // p
        m += 1
    return total


def timed(f, *args):
    t = time.perf_counter()
    r = f(*args)
    return time.perf_counter() - t, r


def c_run(binary, n):
    out = subprocess.run(["./" + binary, str(n), "1.5"], capture_output=True, text=True).stdout.strip()
    return int(out)


subprocess.run(["gcc", "-O3", "-march=native", "-o", "triples", "triples.c", "-lm"], check=True)
subprocess.run(["gcc", "-O3", "-march=native", "-o", "triples_basic", "triples_basic.c", "-lm"], check=True)
print("n enumeration python_basic python_hull c_basic c_hull")
for e in map(int, sys.argv[1:]):
    n = 10 ** e
    row = [str(e)]
    row.append("%.3f" % timed(enumeration, n)[0] if e <= 8 else "-")
    if e <= 12:
        row.append("%.3f" % timed(triples, n, 1.5, "basic")[0])
        row.append("%.3f" % timed(triples, n, 1.5, "hull")[0])
    else:
        row += ["-", "-"]
    tb, vb = timed(c_run, "triples_basic", n)
    th, vh = timed(c_run, "triples", n)
    assert vb == vh
    row += ["%.3f" % tb, "%.3f" % th]
    print(" ".join(row), flush=True)
