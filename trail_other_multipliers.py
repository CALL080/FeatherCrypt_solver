# 16-round differential characteristics for multipliers that do NOT split the cipher.
# Search: dynamic program over characteristics with at most two active S-boxes in every round
# (an upper bound on the minimum weight).  The trail is printed and every transition is re-checked
# with the DDT by a separate routine.  Usage: python trail_other_multipliers.py [r ...]
import math, sys

SBOX = [0xC, 5, 6, 0xB, 9, 0, 0xA, 0xD, 3, 0xE, 0xF, 8, 4, 7, 1, 2]
DDT = [[0] * 16 for _ in range(16)]
for a in range(16):
    for x in range(16):
        DDT[a][SBOX[x] ^ SBOX[x ^ a]] += 1
OUTS = [[(b, 4 - int(math.log2(DDT[a][b]))) for b in range(16) if DDT[a][b]] for a in range(16)]
NIB = 0x1111111111111111
INF = 10 ** 9
ROUNDS = 16


def nact(x):
    return bin((x | x >> 1 | x >> 2 | x >> 3) & NIB).count("1")


def search(r):
    perm = [(r * i) % 64 for i in range(64)]

    def P(x):
        return sum(1 << perm[i] for i in range(64) if (x >> i) & 1)

    PT = [[P(b << (4 * j)) for b in range(16)] for j in range(16)]
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
    parents = []
    for _ in range(ROUNDS):
        g, par = {}, {}
        for s, fs in f.items():
            for n, w in edges.get(s, ()):
                if fs + w < g.get(n, INF):
                    g[n] = fs + w
                    par[n] = s
        f = g
        parents.append(par)
    end = min(f, key=f.get)
    trail = [end]
    for par in reversed(parents):
        trail.append(par[trail[-1]])
    trail.reverse()
    return trail, f[end], perm


def check(trail, perm):
    """Independent re-check: undo P, then look up every S-box transition in the DDT."""
    inv = [0] * 64
    for i in range(64):
        inv[perm[i]] = i
    ws = []
    for a, b in zip(trail, trail[1:]):
        y = sum(((b >> o) & 1) << inv[o] for o in range(64))
        w = 0
        for j in range(16):
            ai, bi = (a >> 4 * j) & 15, (y >> 4 * j) & 15
            assert DDT[ai][bi] > 0, "impossible transition"
            if ai:
                w += 4 - int(math.log2(DDT[ai][bi]))
        ws.append(w)
    return ws


for r in [int(v) for v in sys.argv[1:]] or [43, 3, 31, 63]:
    trail, w, perm = search(r)
    ws = check(trail, perm)
    assert sum(ws) == w
    print(f"r={r}: 16-round characteristic of weight {w}; round weights {ws}")
    print("   input  0x%016X\n   output 0x%016X" % (trail[0], trail[-1]))
    print("   states " + " ".join("%016X" % s for s in trail))
