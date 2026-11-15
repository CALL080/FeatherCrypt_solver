# Which findings survive when the multiplier r is changed?  (r=17 is only "e.g." in the original paper)
# E1 diffusion (bit dependency), E2 avalanche vs. the designers' claims, E3 differential upper bounds, E4 one-round linear claim.
import math, random, sys
SBOX = [0xC, 5, 6, 0xB, 9, 0, 0xA, 0xD, 3, 0xE, 0xF, 8, 4, 7, 1, 2]
M64 = (1 << 64) - 1
INF = 10 ** 9
rng = random.Random(1)
DDT = [[0] * 16 for _ in range(16)]
for a in range(16):
    for x in range(16):
        DDT[a][SBOX[x] ^ SBOX[x ^ a]] += 1
OUTS = [[(b, 4 - int(math.log2(DDT[a][b]))) for b in range(16) if DDT[a][b]] for a in range(16)]
par = lambda v: bin(v).count("1") & 1
LAT = [[sum(1 if par(a & x) == par(b & SBOX[x]) else -1 for x in range(16)) for b in range(16)] for a in range(16)]
mx = max(abs(LAT[a][b]) for a in range(1, 16) for b in range(1, 16))
print("E4 S-box: max |correlation| =", mx, "/16 -> one-round linear trail with squared correlation 2^%d" % (2 * int(math.log2(mx / 16))))
print("   S-box: every output bit depends on every input bit:",
      all(any((SBOX[x] ^ SBOX[x ^ (1 << i)]) >> o & 1 for x in range(16)) for i in range(4) for o in range(4)))
NIB = 0x1111111111111111
def nact(x):
    return bin((x | x >> 1 | x >> 2 | x >> 3) & NIB).count("1")
def study(r):
    perm = [(r * i) % 64 for i in range(64)]
    def P(x):
        y = 0
        i = 0
        while x:
            if x & 1:
                y |= 1 << perm[i]
            x >>= 1
            i += 1
        return y
    PT = [[P(b << (4 * j)) for b in range(16)] for j in range(16)]
    # components
    par_ = list(range(16))
    def find(a):
        while par_[a] != a:
            a = par_[a]
        return a
    for i in range(64):
        par_[find(i // 4)] = find(perm[i] // 4)
    comps = len({find(j) for j in range(16)})
    # E1 dependency: dep[i] = mask of state bits influenced by plaintext bit i
    dep = [1 << i for i in range(64)]
    full = None
    dep3 = None
    for R in range(1, 33):
        nd = []
        for m in dep:
            s = 0
            for j in range(16):
                if (m >> (4 * j)) & 15:
                    s |= PT[j][15]
            nd.append(s)
        dep = nd
        if R == 3:
            dep3 = sum(bin(m).count("1") for m in dep) / 64 / 64
        if full is None and all(m == M64 for m in dep):
            full = R
    # E2 avalanche, 16 rounds, random round keys, one plaintext bit flipped
    SPT = [[PT[2 * k][SBOX[b & 15]] ^ PT[2 * k + 1][SBOX[b >> 4]] for b in range(256)] for k in range(8)]
    def enc(p, rk):
        x = p ^ rk[0]
        for q in range(1, 17):
            y = 0
            for k in range(8):
                y ^= SPT[k][(x >> (8 * k)) & 255]
            x = y ^ rk[q]
        return x
    hd = 0
    T = 4000
    for _ in range(T):
        rk = [rng.getrandbits(64) for _ in range(17)]
        p = rng.getrandbits(64)
        hd += bin(enc(p, rk) ^ enc(p ^ (1 << rng.randrange(64)), rk)).count("1")
    aval = hd / T / 64
    # E3 upper bound on the minimum differential weight: trails with <=2 active S-boxes in every round
    edges = {}
    for j1 in range(16):
        for a1 in range(1, 16):
            s1 = a1 << (4 * j1)
            e = [(PT[j1][b], w) for b, w in OUTS[a1] if nact(PT[j1][b]) <= 2]
            if e:
                edges[s1] = e
            for j2 in range(j1 + 1, 16):
                for a2 in range(1, 16):
                    e = []
                    for b1, w1 in OUTS[a1]:
                        y1 = PT[j1][b1]
                        for b2, w2 in OUTS[a2]:
                            y = y1 ^ PT[j2][b2]
                            if nact(y) <= 2:
                                e.append((y, w1 + w2))
                    if e:
                        edges[s1 | a2 << (4 * j2)] = e
    f = {s: 0 for s in edges}
    ub = []
    for R in range(1, 17):
        g = {}
        for s, fs in f.items():
            for n, w in edges.get(s, ()):
                t = fs + w
                if t < g.get(n, INF):
                    g[n] = t
        f = g
        ub.append(min(g.values()) if g else None)
    order = next(k for k in range(1, 65) if pow(r, k, 64) == 1)
    return order, comps, dep3, full, aval, ub
print("E1-E3   r ord comps dep3%% fullR aval16%%  UB on W*(R) for R=4,8,12,15,16")
for r in range(1, 64, 2):
    o, c, d3, fl, av, ub = study(r)
    print(f"r={r:2d} {o:3d} {c:3d}  {100 * d3:5.1f} {str(fl):>5s}  {100 * av:5.1f}   {ub[3]} {ub[7]} {ub[11]} {ub[14]} {ub[15]}")
    sys.stdout.flush()
