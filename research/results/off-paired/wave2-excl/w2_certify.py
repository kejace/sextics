"""Driver: build sector systems (mod p), reconstruct rational data, and produce EXACT emptiness certificates for every chart.
usage: python3 w2_certify.py --d 8 --chars all|00,20 --job NAME [--prime P] [--tsub all]
Logs one record per chart: result 'PROVEN-empty(exact certificate)' or 'no-certificate(modp-only)'."""
import sys, os, json, time, argparse
sys.path.insert(0, '/tmp/hodge/results/off-paired/wave2-excl')
import w2lib, w2_bilinear as WB, w2_cert as WC, w2_fast
from w2lib import char_of, cadd, monos_char
from hodge.core import monomials

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True); ap.add_argument('--chars', default='all'); ap.add_argument('--job', required=True)
    ap.add_argument('--prime', type=int, default=WB.PM); ap.add_argument('--tsub', default='all'); ap.add_argument('--charts', default='all')
    a = ap.parse_args()
    F, A, q = w2lib.load(); cA = char_of(A)
    chars = [(x, y) for x in range(6) for y in range(6)] if a.chars == 'all' else [(int(s[0]), int(s[1])) for s in a.chars.split(',')]
    recfile = f'{w2lib.PK}/receipts_{a.job}.jsonl'
    for c in chars:
        t0 = time.time()
        sysd, info = WB.build_fast(F, A, q, a.d, c, a.prime, 'sector', a.tsub, addJ=False, verbose=False)
        base = {'exp': 'bilinear-cert', 'job': a.job, 'd': a.d, 'char': c, 'tsub': a.tsub, 'prime': a.prime, 'family': f'all D in sector ({a.tsub})',
                'nt': info.get('nt'), 'ns': info.get('ns'), 'neq': info.get('neq')}
        if sysd is None:
            rec = dict(base, result='empty-t-space(PROVEN: modular colon dim 0 bounds the rational one)', seconds=round(time.time() - t0, 1))
            w2lib.log(rec); open(recfile, 'a').write(json.dumps(rec, default=str) + '\n'); print(json.dumps(rec, default=str), flush=True); continue
        # rational data
        try:
            Dq, Nq, Sq = WC.rational_data(F, A, q, sysd['Ds'], sysd['Ns'], sysd['syz'], a.prime)
        except ValueError as ex:
            rec = dict(base, result=f'rational-data-failed: {ex}', seconds=round(time.time() - t0, 1))
            w2lib.log(rec); open(recfile, 'a').write(json.dumps(rec, default=str) + '\n'); print(json.dumps(rec, default=str), flush=True); continue
        EB = WC.ExactB(Dq, Nq, Sq)
        sector = cadd(cA, cadd(c, c))
        mons2d = monos_char(2 * a.d, sector)
        from hodge.core import shift, reduce_mod
        Fp = reduce_mod(F, a.prime)
        mons_low = monos_char(2 * a.d - 6, sector)
        red = w2_fast.MatReducer([shift(Fp, m, 1, a.prime) for m in mons_low], mons2d, a.prime)
        nt = sysd['nt']
        charts = list(range(nt)) if a.charts == 'all' else [int(x) for x in a.charts.split(',')]
        nprov = 0
        for ch in charts:
            t1 = time.time()
            res = WC.certify_chart(sysd, red, EB, F, mons_low, ch)
            ok = bool(res.get('verified_exact'))
            nprov += ok
            rec = dict(base, chart=ch, cert=res.get('cert'), steps=res.get('steps'), fail=res.get('fail'),
                       result='PROVEN-empty(exact certificate)' if ok else 'no-certificate(modp-only)', seconds=round(time.time() - t1, 1))
            w2lib.log(rec); open(recfile, 'a').write(json.dumps(rec, default=str) + '\n')
        print(json.dumps({'d': a.d, 'char': c, 'nt': nt, 'ns': sysd['ns'], 'charts': len(charts), 'proven': nprov, 'seconds': round(time.time() - t0, 1)}), flush=True)

if __name__ == '__main__':
    main()
