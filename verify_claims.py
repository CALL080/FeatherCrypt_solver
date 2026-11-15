# Checks for the paper: FeatherCrypt with r = 17.
# Round keys are drawn at random: the key schedule is underspecified in the
# original paper and none of the checks below depend on it.
import random, math, sys
from itertools import product

SBOX = [0xC, 5, 6, 0xB, 9, 0, 0xA, 0xD, 3, 0xE, 0xF, 8, 4, 7, 1, 2]
INV = [SBOX.index(i) for i in range(16)]
M64 = (1 << 64) - 1
R = 17


def perm_bits(x, r=R):
    y = 0
    for i in range(64):
        if (x >> i) & 1:
            y |= 1 << ((r * i) % 64)
    return y


def slayer(x):
    y = 0
    for j in range(16):
        y |= SBOX[(x >> (4 * j)) & 15] << (4 * j)
    return y


# byte-indexed SP tables for speed
SPT = [[perm_bits(slayer(b << (8 * k)) & (0xFF << (8 * k))) for b in range(256)] for k in range(8)]
# slayer(b<<8k) puts S(0)=C in the other nibbles; mask keeps only this byte's two nibbles.


def enc64(p, rks):
    x = p ^ rks[0]
    for r in range(1, 17):
        y = 0
        for k in range(8):
            y ^= SPT[k][(x >> (8 * k)) & 0xFF]
        x = y ^ rks[r]
    return x


def enc64_slow(p, rks):
    x = p ^ rks[0]
    for r in range(1, 17):
        x = perm_bits(slayer(x)) ^ rks[r]
    return x


def extract(x, c):  # component V_c -> 16-bit local value (nibble m = S_{c+4m})
    v = 0
    for m in range(4):
        v |= ((x >> (4 * (c + 4 * m))) & 15) << (4 * m)
    return v


def embed(v, c):
    x = 0
    for m in range(4):
        x |= ((v >> (4 * m)) & 15) << (4 * (c + 4 * m))
    return x


def perm16(v):  # local bit (m,b) -> ((m+b)%4, b)
    y = 0
    for m in range(4):
        for b in range(4):
            if (v >> (4 * m + b)) & 1:
                y |= 1 << (4 * ((m + b) % 4) + b)
    return y


def s16(v):
    return sum(SBOX[(v >> (4 * m)) & 15] << (4 * m) for m in range(4))


SP16 = [perm16(s16(v)) for v in range(1 << 16)]
P16 = [perm16(v) for v in range(1 << 16)]

rng = random.Random(20261115)

print("== 0. sanity: fast vs slow 64-bit implementation")
for _ in range(200):
    rks = [rng.getrandbits(64) for _ in range(17)]
    p = rng.getrandbits(64)
    assert enc64(p, rks) == enc64_slow(p, rks)
print("ok")

print("== 1. Theorem 1 (four independent 16-bit components), 2000 random (P,K)")
for _ in range(2000):
    rks = [rng.getrandbits(64) for _ in range(17)]
    p = rng.getrandbits(64)
    c64 = enc64(p, rks)
    for c in range(4):
        x = extract(p, c) ^ extract(rks[0], c)
        for r in range(1, 17):
            x = SP16[x] ^ extract(rks[r], c)
        assert x == extract(c64, c), "decomposition fails"
print("ok: E_K = F0 x F1 x F2 x F3 holds in all trials")

print("== 2. P commutes with every nibble rotation (RotL by 4k bits)")
rot = lambda x, s: ((x << s) | (x >> (64 - s))) & M64 if s else x
assert all(perm_bits(rot(1 << i, 4 * k)) == rot(perm_bits(1 << i), 4 * k) for i in range(64) for k in range(16))
print("ok for all 16 rotations -> every trail has 16 rotated copies")

print("== 3. pre-filter: fraction of random pairs with dP=DC1 passing the '56-bit' zero filter")
DP1 = 0x000B000700000000
keep = sum(1 << b for b in (0, 2, 17, 19, 32, 34, 49, 51))
zero56 = M64 ^ keep
rks = [rng.getrandbits(64) for _ in range(17)]
N = 40000
hit = 0
for _ in range(N):
    p = rng.getrandbits(64)
    d = enc64(p, rks) ^ enc64(p ^ DP1, rks)
    assert d & ~embed(0xFFFF, 0) & M64 == 0  # 48 bits are zero with probability 1
    if d & zero56 == 0:
        hit += 1
print(f"pass rate = {hit}/{N} = 2^{math.log2(hit / N):.2f}   (a filter on 56 bits would give 2^-56; only 8 bits are effective -> 2^-8)")

print("== 4. right pairs of DC1 inside V0, exhaustive over all 2^15 pairs, per key")
DC1 = [0x000B000700000000, 0x0000000900000000, 0x0000000000000004, 0x0000000400000001,
       0x00000000000A0005, 0x0000000000000009, 0x0000000400000000, 0x0000000100000004,
       0x000A000500000000, 0x0000000900000000, 0x0000000000000004, 0x0000000400000001,
       0x00000000000A0005, 0x0000000000000009, 0x0000000400000000, 0x0000000100000004]
assert all(d == embed(extract(d, 0), 0) for d in DC1)
loc = [extract(d, 0) for d in DC1]
z16 = extract(zero56, 0)
KEYS = 60
tot_char = tot_diff = tot_filt = 0
for _ in range(KEYS):
    k = [rng.getrandbits(16) for _ in range(17)]
    states = []
    x = [v ^ k[0] for v in range(1 << 16)]
    follows = [True] * (1 << 16)
    for r in range(1, 17):
        x = [SP16[v] ^ k[r] for v in x]
        if r <= 15:
            d = loc[r]
            a = loc[0]
            for v in range(1 << 16):
                if follows[v] and (x[v] ^ x[v ^ a]) != d:
                    follows[v] = False
            if r == 15:
                tot_diff += sum(1 for v in range(1 << 16) if (x[v] ^ x[v ^ loc[0]]) == loc[15]) // 2
    tot_char += sum(follows) // 2
    tot_filt += sum(1 for v in range(1 << 16) if (x[v] ^ x[v ^ loc[0]]) & z16 == 0) // 2
print(f"keys={KEYS}: pairs following the whole 15-round characteristic: total {tot_char}")
print(f"          pairs with dX15 = DC1 output (differential): mean {tot_diff / KEYS:.3f} per key (16-bit random perm: 0.5)")
print(f"          pairs passing the ciphertext filter: mean {tot_filt / KEYS:.1f} of 32768 per key (2^15 * 2^-8 = 128)")

print("== 5. exact minimum weights by DP inside one component (independent of STP)")
INF = 10 ** 9
DDT = [[0] * 16 for _ in range(16)]
for a in range(16):
    for x in range(16):
        DDT[a][SBOX[x] ^ SBOX[x ^ a]] += 1
W = [[(4 - int(math.log2(DDT[a][b]))) if DDT[a][b] else INF for b in range(16)] for a in range(16)]
OUTS = [[(b, W[a][b]) for b in range(16) if DDT[a][b]] for a in range(16)]


def round_dp(f, cnt):
    # S-layer, one nibble at a time (min-plus, with path counts), then local permutation
    for m in range(4):
        g = [INF] * 65536
        gc = [0] * 65536
        sh = 4 * m
        for v in range(65536):
            fv = f[v]
            if fv >= INF:
                continue
            a = (v >> sh) & 15
            base = v & ~(15 << sh)
            cv = cnt[v]
            for b, w in OUTS[a]:
                u = base | (b << sh)
                t = fv + w
                if t < g[u]:
                    g[u] = t
                    gc[u] = cv
                elif t == g[u]:
                    gc[u] += cv
        f, cnt = g, gc
    g = [INF] * 65536
    gc = [0] * 65536
    for v in range(65536):
        g[P16[v]] = f[v]
        gc[P16[v]] = cnt[v]
    return g, gc


f = [0] * 65536
cnt = [1] * 65536
f[0] = INF
cnt[0] = 0
paper = [2, 4, 8, 12, 16, 19, 23, 27, 31, 34, 38, 42, 46, 49, 53, 57]
for r in range(1, 17):
    f, cnt = round_dp(f, cnt)
    w = min(f)
    n = sum(c for fv, c in zip(f, cnt) if fv == w)
    flag = "OK" if w == paper[r - 1] else "MISMATCH"
    print(f"R={r:2d}  W*={w:2d}  paper={paper[r - 1]:2d}  {flag}   optimal trails in one component: {n}  (x4 in the cipher: {4 * n})")

print("== 6. all 32 admissible multipliers: independent components of the S-box graph")
for r in range(1, 64, 2):
    parent = list(range(16))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i in range(64):
        a, b = find(i // 4), find(((r * i) % 64) // 4)
        parent[a] = b
    comps = {}
    for j in range(16):
        comps.setdefault(find(j), []).append(j)
    sizes = sorted(len(v) for v in comps.values())
    order = next(k for k in range(1, 65) if pow(r, k, 64) == 1)
    if len(comps) > 1:
        print(f"r={r:2d} order={order:2d}: {len(comps)} components of {sizes[0] * 4} bits")
print("(multipliers not listed: one 64-bit component)")

