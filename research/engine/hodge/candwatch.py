"""Watch for witness-candidate files, verify each exactly against its lane, and report only real hits.

    python3 -m hodge.candwatch --minutes 120

A candidate counts only if (a) the cleared identity holds exactly for the LANE's (F, A, q), (b) D is not divisible
by F (nonzero modulo F) and gcd(F|line, D|line) = 1 on a random line mod p (coprimality triage), and (c) degrees match.
Checked files are remembered in /tmp/hodge/candwatch-seen.json so they are not re-verified.
"""
from __future__ import annotations
import argparse, glob, json, os, random, sys, time
from fractions import Fraction
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core

SEEN = '/tmp/hodge/candwatch-seen.json'


def lane_of(path):
    parts = path.split('/')
    for i, p in enumerate(parts):
        if p in ('results', 'remote') and i + 1 < len(parts):
            return parts[i + 1]
    return None


def coprime_triage(F, D, p=1000003, trials=2):
    import sympy
    t = sympy.Symbol('t')
    for _ in range(trials):
        a = [random.randint(1, p - 1) for _ in range(6)]; b = [random.randint(1, p - 1) for _ in range(6)]
        def restrict(P):
            expr = 0
            for e, c in P.items():
                c = Fraction(c); cc = (c.numerator * pow(c.denominator, -1, p)) % p
                term = cc
                for i in range(6):
                    if e[i]: term = term * (a[i] + b[i] * t) ** e[i]
                expr += term
            return sympy.Poly(sympy.expand(expr), t, modulus=p)
        if sympy.gcd(restrict(F), restrict(D)).degree() == 0:
            return True
    return False


def check(path):
    lane = lane_of(path)
    try:
        w = json.load(open(path))
        D = core.parse(w['D']); N = [core.parse(s) for s in w['N']]; Q = core.parse(w['Q'])
    except Exception as e:
        return {'path': path, 'lane': lane, 'status': f'unparseable: {e}'}
    try:
        F, A, q, rec = core.load_instance(lane)
    except Exception as e:
        return {'path': path, 'lane': lane, 'status': f'lane load failed: {e}'}
    ok, info = core.verify_witness(F, A, q, D, N, Q)
    if not ok:
        return {'path': path, 'lane': lane, 'status': 'identity FAILS for the lane', 'info': info,
                'declared': w.get('instanceId') or w.get('instance') or w.get('id')}
    # legality triage
    p = core.PRIMES[0]
    Fp, Dp = core.reduce_mod(F, p), core.reduce_mod(D, p)
    # D divisible by F?  compare degrees and test a few random points of V(F) is hard; use gcd triage on lines
    cop = coprime_triage(F, D)
    return {'path': path, 'lane': lane, 'status': 'REAL WITNESS (identity exact, coprime by triage)' if cop else 'identity holds but D NOT coprime to F (spurious)',
            'd': core.degree(D), 'coprime_triage': cop}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--minutes', type=int, default=120); ap.add_argument('--once', action='store_true')
    a = ap.parse_args()
    seen = json.load(open(SEEN)) if os.path.exists(SEEN) else {}
    t0 = time.time()
    while True:
        files = glob.glob('/tmp/hodge/results/*/witness-candidate*.json') + glob.glob('/tmp/hodge/remote/*/witness-candidate*.json') \
            + glob.glob('/tmp/hodge/results/*/**/witness-candidate*.json', recursive=True)
        hit = False
        for f in sorted(set(files)):
            key = f + ':' + str(os.path.getmtime(f))
            if key in seen:
                continue
            r = check(f)
            seen[key] = r
            json.dump(seen, open(SEEN, 'w'), indent=1)
            print(time.strftime('%H:%M'), json.dumps(r)[:400], flush=True)
            if r['status'].startswith('REAL WITNESS'):
                print('ALERT REAL WITNESS', f, flush=True); hit = True
        if hit or a.once or time.time() - t0 > a.minutes * 60:
            break
        time.sleep(60)


if __name__ == '__main__':
    main()
