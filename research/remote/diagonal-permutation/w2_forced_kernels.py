"""Equivariant Hilbert functions of the Jacobian ring R = S/J(F) for F smooth even sextic with symmetry group
G = S3 x <tau> (simultaneous permutations of the pairs (X_i,Y_i); tau: Y -> -Y) and the forced kernels of multiplication by a
Y-odd, (01)-anti-invariant target A (such as A = X0^3 Y0^3 - X1^3 Y1^3) per (Y-parity, (01)-eigenvalue) sector.
Koszul series: sum_d tr(g | R_d) t^d = det(1 - t^5 g) / det(1 - t g)  (regular sequence of six quintics).
Reproduces: k_9 = 21 (sigma-invariant, both Y-parities), k_10^+ = 170 + 124 = 294, k_10^- = 152 + 109 (+1 from D8) = 262."""
import numpy as np
N = 30
def series(num, den):
    s = np.zeros(N + 1); s[0] = 1
    mul = lambda a, b: np.convolve(a, b)[:N + 1]
    for f in num: s = mul(s, np.array(f + [0] * (N + 1 - len(f))))
    for f in den:
        inv = np.zeros(N + 1); inv[0] = 1; f = np.array(f + [0] * (N + 1 - len(f)))
        for k in range(1, N + 1): inv[k] = -sum(f[j] * inv[k - j] for j in range(1, k + 1))
        s = mul(s, inv)
    return s
m5 = [1, 0, 0, 0, 0, -1]; p5 = [1, 0, 0, 0, 0, 1]; m1 = [1, -1]; p1 = [1, 1]
m3 = [1, 0, 0, -1]; p3 = [1, 0, 0, 1]; m15 = [1] + [0] * 14 + [-1]; p15 = [1] + [0] * 14 + [1]
H = {'e': series([m5] * 6, [m1] * 6), 't': series([m5] * 3 + [p5] * 3, [m1] * 3 + [p1] * 3),
     's': series([m5] * 4 + [p5] * 2, [m1] * 4 + [p1] * 2), 'st': series([m5] * 3 + [p5] * 3, [m1] * 3 + [p1] * 3),
     'c': series([m15] * 2, [m3] * 2), 'ct': series([m15, p15], [m3, p3])}
def sector(d, ypar, sig):
    """dim of the (Y-parity ypar (0 even / 1 odd), sigma-eigenvalue sig (+1/-1)) part of R_d."""
    yt = 1 if ypar == 0 else -1
    return int(round((H['e'][d] + yt * H['t'][d] + sig * (H['s'][d] + yt * H['st'][d])) / 4))
if __name__ == '__main__':
    print('d: R_d^{+,s+} R_d^{+,s-} R_d^{-,s+} R_d^{-,s-}')
    for d in range(0, 19):
        print(d, sector(d, 0, 1), sector(d, 0, -1), sector(d, 1, 1), sector(d, 1, -1))
    print('forced kernel of [A] (A Y-odd, sigma-anti): R_d^{e,s} -> R_{d+6}^{1-e,-s}')
    for d in range(4, 13):
        f = {(e, s): max(0, sector(d, e, s) - sector(d + 6, 1 - e, -s)) for e in (0, 1) for s in (1, -1)}
        print(f'd={d}: even D: sigma+ {f[(0,1)]}, sigma- {f[(0,-1)]} (total {f[(0,1)]+f[(0,-1)]}) | odd D: sigma+ {f[(1,1)]}, sigma- {f[(1,-1)]} (total {f[(1,1)]+f[(1,-1)]})')
