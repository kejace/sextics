"""Batch 1: exact E1-admissible families per (blocks, chi, d, pD); logs dims."""
import sys, time, json, argparse
sys.path.insert(0, '/tmp/hodge/results/symmetric-sextic-residual')
import w2lib as W
ap = argparse.ArgumentParser()
ap.add_argument('--blocks', default='S3S3'); ap.add_argument('--dmin', type=int, default=4); ap.add_argument('--dmax', type=int, default=12)
ap.add_argument('--chis', default='triv,sgnX,sgnY,sgnsgn'); ap.add_argument('--pDs', default='0,1')
args = ap.parse_args()
F, A, q = W.load()
rows = []
for d in range(args.dmin, args.dmax + 1):
    for chi in args.chis.split(','):
        for pD in [int(x) for x in args.pDs.split(',')]:
            t0 = time.time()
            fam = W.family(F, A, A and q, d, pD, chi, args.blocks, verbose=True)
            rec = {'batch': 'dims', 'family': f'{args.blocks}/{chi}', 'd': d, 'pD': pD, 'ker': fam['ker'], 'dimJ': fam['dimJ'],
                   'nD': fam['nD'], 'K': fam['ker'] + fam['dimJ'], 'ratrecon_ok': fam['ratrecon_ok'],
                   'kernel_verified_modp': fam.get('kernel_verified_modp'), 'primes_used': len(fam['primes_used']),
                   'rows_e1': fam['rows_e1'], 'nN_e1': fam['nN_e1'], 'seconds': round(time.time() - t0, 1)}
            W.wlog(rec); print(json.dumps(rec), flush=True)
