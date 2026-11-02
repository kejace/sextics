"""Degree-6 relaxation restricted to unions of character sectors (lens excl, lane cycle-q2).
For a set T of sectors, V_T = relaxation variety of denominators supported on T.  Reports, per T, the affine dimension
of V_T and (when T contains the trivial sector) whether V_T meets {c_a != 0} for a nontrivial variable (saturation).
usage: python3 wave2-excl/subsets6.py --size 2|3 [--batch 400]
"""
import json, subprocess, re, collections, itertools, time, argparse, os, sys
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--size', type=int, default=2); ap.add_argument('--batch', type=int, default=400); ap.add_argument('--prime', type=int, default=1000000007)
    a = ap.parse_args(); t0 = time.time()
    Q = json.load(open(f'{w2lib.PACKET}/hess6_p{a.prime}.quadrics.json'))
    keys = [tuple(k) for k in Q['basis_keys']]; p = Q['prime']
    bysec = collections.defaultdict(list)
    for i, k in enumerate(keys): bysec[k].append(i)
    secs = sorted(bysec); triv = (0, 0, 0)
    quads = [({tuple(map(int, kk.split(','))): v for kk, v in q['terms'].items()}) for q in Q['quadrics']]
    vars_ = ','.join(f'c{i+1}' for i in range(62))
    subsets = list(itertools.combinations(secs, a.size))
    results = []
    for start in range(0, len(subsets), a.batch):
        chunk = subsets[start:start + a.batch]
        lines = [f'LIB "elim.lib";\nring r = {p}, ({vars_}), dp;']
        for T in chunk:
            sidx = set().union(*[set(bysec[k]) for k in T])
            polys = []
            for terms in quads:
                tt = {ab: v for ab, v in terms.items() if ab[0] in sidx and ab[1] in sidx}
                if tt: polys.append('+'.join(f'{v}*c{x+1}*c{y+1}' if x != y else f'{v}*c{x+1}^2' for (x, y), v in sorted(tt.items())))
            others = [f'c{j+1}' for j in range(62) if j not in sidx]
            lines.append('ideal I = ' + (',\n'.join(polys) if polys else '0') + (',' + ','.join(others) if others else '') + ';')
            tag = '|'.join(''.join(map(str, k)) for k in T)
            lines.append(f'ideal G = std(I); "RES {tag} DIM", dim(G);')
            if triv in T:
                nontriv = [f'c{i+1}' for k in T if k != triv for i in bysec[k]]
                lines.append('poly f = ' + '*'.join(nontriv) + ';')
                lines.append('def S = sat(I, f); ideal Sat; if (typeof(S) == "list") { Sat = S[1]; } else { Sat = S; }')
                lines.append(f'"SAT {tag} DIM", dim(std(Sat));')
        lines.append('quit;')
        path = f'/tmp/subsets6_{a.size}_{start}.sing'
        open(path, 'w').write('\n'.join(lines) + '\n')
        r = subprocess.run(['Singular', '-q', path], capture_output=True, text=True, timeout=3000)
        res = re.findall(r'RES (\S+) DIM (-?\d+)', r.stdout); sat = dict(re.findall(r'SAT (\S+) DIM (-?\d+)', r.stdout))
        for tag, dm in res:
            results.append({'T': tag, 'dim': int(dm), 'sat_dim': int(sat[tag]) if tag in sat else None})
        os.remove(path)
        print(f'{start + len(chunk)}/{len(subsets)} done [{time.time()-t0:.0f}s]', flush=True)
    dims = collections.Counter(r['dim'] for r in results)
    nontriv_bad = [r for r in results if r['sat_dim'] is None and r['dim'] > 0]
    triv_bad = [r for r in results if r['sat_dim'] is not None and r['sat_dim'] >= 0]
    out = {'prime': p, 'size': a.size, 'n': len(subsets), 'dims': dict(dims), 'nontrivial_sets_with_nonempty_variety': nontriv_bad[:50],
           'trivial_containing_sets_with_points_off_V0': triv_bad[:50], 'n_nontriv_bad': len(nontriv_bad), 'n_triv_bad': len(triv_bad), 'seconds': round(time.time() - t0, 1)}
    json.dump({'summary': out, 'results': results}, open(f'{w2lib.PACKET}/sectors/subsets_d6_size{a.size}_p{p}.json', 'w'))
    verdict = (f'all {len(subsets)} unions of {a.size} sectors: nontrivial unions have empty relaxation variety' if not nontriv_bad else f'{len(nontriv_bad)} nontrivial unions have NONEMPTY relaxation variety') + \
              ('; unions containing the trivial sector have no points outside V0 (saturation empty)' if not triv_bad else f'; {len(triv_bad)} trivial-containing unions have points with a nontrivial coordinate nonzero')
    w2lib.log(exp=f'subsets6-size{a.size}', d=6, prime=p, result=verdict, rank=dict(dims), seconds=out['seconds'], job=os.environ.get('W2JOB', 'subsets6'))
    print(json.dumps(out, indent=1))

if __name__ == '__main__':
    main()
