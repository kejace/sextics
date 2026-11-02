"""Branch-and-kill for a homogeneous quadric ideal I (JSON from d7mixed.py or family exports):
 * c_a^2 in I_2  or  c_a^3 in I_3  =>  c_a = 0 on V(I);
 * a monomial c_a c_b (in I_2) or c_a c_b c_c (in I_3) with distinct variables => V subset union of the coordinate hyperplanes
   => branch.  Each child substitutes one variable = 0 (quadrics lose those terms).  Memoized on the killed set.
Monomials of I_k are read off the rref (singleton rows).  Everything is exact linear algebra mod p (python-flint).
Usage: branchkill.py quads.json out_prefix [max_nodes] [use_I3=1]"""
import sys, json, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import nmod_mat

def load(fn):
    J = json.load(open(fn))
    quads = [{(k, l): v for k, l, v in cf} for cf in J['quads']]
    chi = [tuple(c) for c in J['chi']]
    return quads, chi, J.get('names')

def chi_add(a, b):
    return ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2, (a[2] + b[2]) % 6, (a[3] + b[3]) % 2)

class BK:
    def __init__(self, quads, chi, p, out_prefix, max_nodes=2000, use_I3=True, log=print):
        self.Q0, self.chi, self.p, self.out, self.max_nodes, self.use_I3, self.log = quads, chi, p, out_prefix, max_nodes, use_I3, log
        self.n = len(chi)
        self.seen = {}
        self.nodes = 0
        self.residuals = []
        self.t0 = time.time()
    def restrict(self, killed):
        Q = []
        for cf in self.Q0:
            cf2 = {(k, l): v for (k, l), v in cf.items() if k not in killed and l not in killed}
            if cf2: Q.append(cf2)
        return Q
    def I2_monomials(self, Q, alive):
        p = self.p
        pairs = [(a, b) for i, a in enumerate(alive) for b in alive[i:]]
        pidx = {pr: i for i, pr in enumerate(pairs)}
        if not Q: return [], 0, pairs
        flat = [0] * (len(Q) * len(pairs))
        for r, cf in enumerate(Q):
            for pr, v in cf.items(): flat[r * len(pairs) + pidx[pr]] = v
        M = nmod_mat(len(Q), len(pairs), flat, p)
        R, rk = M.rref()
        monos = []
        npairs = len(pairs)
        for i in range(rk):
            j = 0
            while int(R[i, j]) == 0: j += 1
            # singleton row <=> no other nonzero after pivot
            single = True
            for jj in range(j + 1, npairs):
                if int(R[i, jj]): single = False; break
            if single: monos.append(pairs[j])
        return monos, rk, pairs
    def I3_monomials(self, Q, alive):
        p = self.p
        chi_q = [chi_add(self.chi[k], self.chi[l]) for cf in Q for (k, l) in [next(iter(cf))]]
        blocks = {}
        for r, cq in enumerate(chi_q):
            for v in alive:
                blocks.setdefault(chi_add(cq, self.chi[v]), []).append((v, r))
        monos = []
        quot = 0
        for cc, lst in blocks.items():
            cols = {}
            rows = []
            for (v, r) in lst:
                row = {}
                for (k, l), val in Q[r].items():
                    key = tuple(sorted((v, k, l)))
                    j = cols.setdefault(key, len(cols))
                    row[j] = (row.get(j, 0) + val) % p
                rows.append(row)
            inv = {i: k for k, i in cols.items()}
            nc = len(cols)
            flat = [0] * (len(rows) * nc)
            for r, row in enumerate(rows):
                for j, val in row.items(): flat[r * nc + j] = val
            M = nmod_mat(len(rows), nc, flat, p)
            R, rk = M.rref()
            quot += nc - rk
            for i in range(rk):
                j = 0
                while int(R[i, j]) == 0: j += 1
                single = True
                for jj in range(j + 1, nc):
                    if int(R[i, jj]): single = False; break
                if single: monos.append(inv[j])
        return monos, quot
    def solve(self, killed=frozenset(), depth=0):
        """returns True if V(I restricted) = {0} is certified on this branch, False if stuck (residual exported)."""
        if killed in self.seen: return self.seen[killed]
        self.nodes += 1
        if self.nodes > self.max_nodes:
            self.residuals.append(sorted(killed)); return False
        killed = set(killed)
        while True:
            alive = [a for a in range(self.n) if a not in killed]
            if not alive:
                self.seen[frozenset(killed)] = True; return True
            Q = self.restrict(killed)
            monos2, rk2, pairs = self.I2_monomials(Q, alive)
            sq = {a for (a, b) in monos2 if a == b}
            if sq:
                killed |= sq; continue
            monos3, quot3 = ([], None)
            if self.use_I3:
                monos3, quot3 = self.I3_monomials(Q, alive)
                cubes = {m[0] for m in monos3 if m[0] == m[1] == m[2]}
                if cubes:
                    killed |= cubes; continue
            break
        alive = [a for a in range(self.n) if a not in killed]
        # branching candidates: monomials with 2 distinct vars first
        cands = [set(m) for m in monos2 if m[0] != m[1]] + [set(m) for m in monos3 if len(set(m)) >= 2]
        cands = [c for c in cands if all(a in alive for a in c)]
        self.log(f'node {self.nodes} depth {depth} killed {len(killed)} alive {len(alive)} quads {len(Q)} rank2 {rk2} I2monos {len(monos2)} I3monos {len(monos3)} quot3 {quot3} cands {len(cands)} ({time.time()-self.t0:.0f}s)')
        if not cands:
            self.residuals.append(sorted(killed))
            # export residual system for an external Groebner run
            tag = f'{self.out}.residual{len(self.residuals)}'
            kidx = {k: i for i, k in enumerate(alive)}
            eqs = []
            for cf in Q:
                terms = [(f'{v}*c{k}^2' if k == l else f'{v}*c{k}*c{l}') for (k, l), v in cf.items()]
                eqs.append('+'.join(terms))
            with open(tag + '.ms', 'w') as f:
                f.write(','.join(f'c{k}' for k in alive) + '\n' + str(self.p) + '\n' + ',\n'.join(eqs) + '\n')
            json.dump({'killed': sorted(killed), 'alive': alive}, open(tag + '.json', 'w'))
            self.seen[frozenset(killed)] = False
            return False
        best = min(cands, key=len)
        ok = True
        for a in sorted(best):
            r = self.solve(frozenset(killed | {a}), depth + 1)
            ok = ok and r
        self.seen[frozenset(killed)] = ok
        return ok

if __name__ == '__main__':
    fn, out = sys.argv[1], sys.argv[2]
    max_nodes = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
    use_I3 = (sys.argv[4] != '0') if len(sys.argv) > 4 else True
    J = json.load(open(fn))
    p = J.get('p') or int(os.environ.get('HODGE_P', '2147483629'))
    quads, chi, names = load(fn)
    bk = BK(quads, chi, p, out, max_nodes=max_nodes, use_I3=use_I3, log=lambda s: print(s, file=sys.stderr, flush=True))
    ok = bk.solve()
    res = {'certified_only_origin': ok, 'nodes': bk.nodes, 'n_residuals': len(bk.residuals), 'residuals': bk.residuals[:20], 'seconds': round(time.time() - bk.t0, 1)}
    print(json.dumps(res))
    json.dump(res, open(out + '.result.json', 'w'))
