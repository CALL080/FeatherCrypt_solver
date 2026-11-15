# Experiments for the paper: full 16-round FeatherCrypt (r = 17), 64-bit data path.
# The round keys are 17 independent uniformly random 64-bit words.  The attacks do not use the key
# schedule, so no particular key schedule (64-bit or 80-bit key, round constants) has to be fixed.
import math
import sys
import time
import numpy as np

print(f"Python {sys.version.split()[0]}, NumPy {np.__version__}")

SBOX = [0xC, 5, 6, 0xB, 9, 0, 0xA, 0xD, 3, 0xE, 0xF, 8, 4, 7, 1, 2]
R = 17
rng = np.random.default_rng(20261115)
U = np.uint64


def perm_bits(x, perm):
    y = 0
    for i in range(64):
        if (x >> i) & 1:
            y |= 1 << perm[i]
    return y


def make_tables(perm):
    t = np.zeros((8, 256), dtype=np.uint64)
    for k in range(8):
        for b in range(256):
            s = SBOX[b & 15] | (SBOX[b >> 4] << 4)
            t[k, b] = perm_bits(s << (8 * k), perm)
    return t


PERM = [(R * i) % 64 for i in range(64)]
SPT = make_tables(PERM)


def encrypt(p, rk, spt=SPT, rounds=16):
    x = p ^ rk[0]
    for r in range(1, rounds + 1):
        y = np.zeros_like(x)
        for k in range(8):
            y ^= spt[k][((x >> U(8 * k)) & U(255)).astype(np.intp)]
        x = y ^ rk[r]
    return x


def extract(x, c):  # 16-bit component V_c of a 64-bit word
    v = np.zeros(x.shape, dtype=np.uint64)
    for m in range(4):
        v |= ((x >> U(4 * c + 16 * m)) & U(15)) << U(4 * m)
    return v


def embed(v, c):
    x = np.zeros(v.shape, dtype=np.uint64)
    for m in range(4):
        x |= ((v >> U(4 * m)) & U(15)) << U(4 * c + 16 * m)
    return x


MASK = [int(embed(np.array([0xFFFF], dtype=np.uint64), c)[0]) for c in range(4)]
ALL = (1 << 64) - 1
rand64 = lambda n: rng.integers(0, 1 << 64, size=n, dtype=np.uint64)
rk = rand64(17)

print("== A. two-query distinguisher (full 16 rounds)")
N = 100000
p = rand64(N)
ok = 0
for c in range(4):
    d = embed(rng.integers(1, 1 << 16, size=N, dtype=np.uint64), c)
    dc = encrypt(p, rk) ^ encrypt(p ^ d, rk)
    ok += int(np.count_nonzero((dc & U(ALL ^ MASK[c])) == 0))
print(f"pairs whose ciphertext difference stays inside the active component: {ok}/{4 * N}")

print("== B. decryption of any ciphertext with 2^16 chosen plaintexts")
t0 = time.perf_counter()
x = np.arange(1 << 16, dtype=np.uint64)
queries = sum(embed(x, c) for c in range(4))  # P(x) = (x, x, x, x)
ct = encrypt(queries, rk)
T = []
for c in range(4):
    inv = np.zeros(1 << 16, dtype=np.uint16)  # 2^16 entries of 16 bits
    fc = extract(ct, c)
    assert len(np.unique(fc)) == 1 << 16  # each F_c is a permutation of 16-bit values
    inv[fc.astype(np.intp)] = x.astype(np.uint16)
    T.append(inv)
t1 = time.perf_counter()
M = 1000000
secret_p = rand64(M)
target_c = encrypt(secret_p, rk)
t2 = time.perf_counter()
rec = np.zeros(M, dtype=np.uint64)
for c in range(4):
    rec |= embed(T[c][extract(target_c, c).astype(np.intp)].astype(np.uint64), c)
t3 = time.perf_counter()
print(f"table construction: {t1 - t0:.3f} s for 2^16 queries; "
      f"the four tables occupy {sum(t.nbytes for t in T)} bytes (4 x 2^16 entries of 16 bits = 512 KiB)")
print(f"decrypted {int(np.count_nonzero(rec == secret_p))}/{M} random ciphertexts correctly in {t3 - t2:.3f} s")

print("== C. known-plaintext recovery: number of uniform random plaintexts until all four tables are complete")
TR = 200
need = []
for _ in range(TR):
    pl = rand64(1 << 21)
    worst = 0
    for c in range(4):
        vals, first = np.unique(extract(pl, c), return_index=True)
        assert len(vals) == 1 << 16
        worst = max(worst, int(first.max()) + 1)
    need.append(worst)
need = np.array(need)
N99 = math.ceil(65536 * math.log(4 * 65536 / 0.01))  # bound of the paper for failure probability 0.01
print(f"trials={TR}: mean 2^{np.log2(need.mean()):.2f}, median 2^{np.log2(np.median(need)):.2f}, "
      f"max 2^{np.log2(need.max()):.2f}; fraction within the bound N = {N99} (2^{math.log2(N99):.1f}): "
      f"{np.mean(need <= N99):.3f}")

print("== D. avalanche of the r = 17 instance: one plaintext bit flipped, 16 rounds")
N = 200000
p = rand64(N)
bit = rng.integers(0, 64, size=N).astype(np.uint64)
dc = encrypt(p, rk) ^ encrypt(p ^ (U(1) << bit), rk)
hw = np.array([bin(int(v)).count("1") for v in dc])
print(f"changed ciphertext bits: mean {hw.mean():.2f} of 64 ({100 * hw.mean() / 64:.1f}%), max {hw.max()} (bound 16 = 25%)")

print("== E. bit numbering: bit 0 = MSB instead of LSB")
perm_msb = [63 - ((R * (63 - q)) % 64) for q in range(64)]
cls = {}
for q in range(64):
    cls.setdefault((q // 4) % 4, set()).add((perm_msb[q] // 4) % 4)
print("S-box index mod 4 is preserved under MSB-first numbering:", all(v == {k} for k, v in cls.items()))
spt_m = make_tables(perm_msb)
p = rand64(N)
d = embed(rng.integers(1, 1 << 16, size=N, dtype=np.uint64), 2)
dc = encrypt(p, rk, spt_m) ^ encrypt(p ^ d, rk, spt_m)
print("two-query distinguisher under MSB-first numbering:", int(np.count_nonzero((dc & U(ALL ^ MASK[2])) == 0)), "/", N)

print("== F. decryption from known plaintexts only (no chosen query)")
kp = rand64(N99)
kc = encrypt(kp, rk)
TK, filled = [], 0
for c in range(4):
    inv = np.zeros(1 << 16, dtype=np.uint16)
    seen = np.zeros(1 << 16, dtype=bool)
    idx = extract(kc, c).astype(np.intp)
    inv[idx] = extract(kp, c).astype(np.uint16)
    seen[idx] = True
    filled += int(seen.sum())
    TK.append(inv)
M2 = 100000
fresh_p = rand64(M2)
fresh_c = encrypt(fresh_p, rk)
rec = np.zeros(M2, dtype=np.uint64)
for c in range(4):
    rec |= embed(TK[c][extract(fresh_c, c).astype(np.intp)].astype(np.uint64), c)
print(f"{N99} known pairs filled {filled} of {4 * 65536} table entries; "
      f"decrypted {int(np.count_nonzero(rec == fresh_p))}/{M2} fresh ciphertexts correctly")
