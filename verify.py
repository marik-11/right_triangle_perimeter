import random
import subprocess
from math import isqrt
from brute import triples_upto
from right_triangle_perimeter import triples, hull_sum, hull_sum_root, sum_floor

REFERENCE = [0, 17, 325, 4858, 64741, 808950, 9706567, 113236940, 1294080089, 14557915466]

subprocess.run(["gcc", "-O3", "-march=native", "-o", "triples", "triples.c", "-lm"], check=True)
subprocess.run(["gcc", "-O3", "-march=native", "-o", "triples_basic", "triples_basic.c", "-lm"], check=True)


def c_value(binary, n, ratio=1.5):
    return int(subprocess.run(["./" + binary, str(n), str(ratio)], capture_output=True, text=True, check=True).stdout)


bad = 0
for m in range(1, 400):
    r = isqrt(m)
    for a in range(1, r + 1):
        for b in range(a, r + 1):
            d = sum_floor(m, a, b)
            if hull_sum(m, a, b) != d or hull_sum_root(m, a, b) != d:
                bad += 1
for n in range(1, 400):
    v = triples_upto(n)
    if v != triples(n) or v != triples(n, 1.5, "basic") or v != c_value("triples", n) or v != c_value("triples_basic", n):
        bad += 1
for k, v in enumerate(REFERENCE, 1):
    if c_value("triples", 10 ** k) != v or c_value("triples_basic", 10 ** k) != v:
        bad += 1
    if k <= 9 and triples(10 ** k) != v:
        bad += 1
random.seed(1)
for _ in range(40):
    n = random.randrange(10 ** 5, 10 ** 10)
    if not triples(n) == triples(n, 1.5, "basic") == c_value("triples", n) == c_value("triples_basic", n):
        bad += 1
for _ in range(60):
    n = random.randrange(10 ** 10, 10 ** 13)
    if c_value("triples", n) != c_value("triples_basic", n):
        bad += 1
for ratio in (0.5, 1.5, 3.0):
    n = random.randrange(10 ** 11, 10 ** 12)
    if c_value("triples", n, ratio) != c_value("triples", n, 1.0):
        bad += 1
print("ok" if bad == 0 else "FAIL %d" % bad)
