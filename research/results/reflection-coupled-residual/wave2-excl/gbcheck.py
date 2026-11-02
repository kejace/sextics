"""msolve GB (mod p) of an exported quadratic system + Singular dim/mult/Hilbert function => receipt.
Usage: python3 gbcheck.py <file.ms> <threads> [timeout_s] [maxdeg]
Verdict 'zero_dim_cone' (dim 0 => V(I)={0}): the ideal contains all forms of degree k0 = first zero of HF mod p,
hence over Q (rank argument with good reduction), hence V(I_Q) = {0}: only trivial solutions of the Hessian relaxation."""
import sys, json, time, subprocess, os
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl')
import w2lib as L
fn = sys.argv[1]; thr = int(sys.argv[2]); tmo = int(sys.argv[3]) if len(sys.argv) > 3 else 3600
maxdeg = int(sys.argv[4]) if len(sys.argv) > 4 else 10
meta = json.load(open(fn.replace('.ms', '.json')))
t0 = time.time()
gbf = fn.replace('.ms', '.gb')
try:
    r = subprocess.run(['msolve', '-f', fn, '-o', gbf, '-g', '2', '-t', str(thr), '-v', '1'], capture_output=True, text=True, timeout=tmo)
    msolve_log = r.stdout[-1500:] + r.stderr[-500:]
    tgb = round(time.time() - t0, 1)
except subprocess.TimeoutExpired:
    L.receipt('gbcheck', exp='gb-hilbert', file=os.path.basename(fn), d=meta['d'], prime=meta['p'], chis=meta['chis'], dim_W=meta['n'],
              n_quadrics=meta['n_quadrics'], rank_quadrics=meta['rank_quadrics'], result='msolve TIMEOUT', timeout_s=tmo)
    sys.exit(0)
vars_ = open(fn).readline().strip().split(',')
txt = open(gbf).read(); body = txt[txt.find('[') + 1: txt.rfind(']')]
polys = [t.strip() for t in body.replace('\n', '').split(',') if t.strip()]
sing = f"""ring r = {meta['p']}, ({','.join(vars_)}), dp;
ideal I = {','.join(polys)};
attrib(I,"isSB",1);
print("DIM"); print(dim(I)); print("MULT"); print(mult(I)); print("HF");
int k; for (k=0; k<={maxdeg}; k++) {{ print(size(kbase(I,k))); }}
quit;
"""
rs = subprocess.run(['Singular', '-q'], input=sing, capture_output=True, text=True, timeout=tmo)
toks = rs.stdout.split()
try:
    dimI = int(toks[toks.index('DIM') + 1]); multI = int(toks[toks.index('MULT') + 1])
    HF = [int(x) for x in toks[toks.index('HF') + 1:]]
except Exception:
    dimI, multI, HF = None, None, toks[-20:]
k0 = next((k for k, v in enumerate(HF) if v == 0), None) if isinstance(HF, list) and all(isinstance(v, int) for v in HF) else None
verdict = 'zero_dim_cone (V(I)={0}: only trivial solutions)' if dimI == 0 else f'POSITIVE-DIMENSIONAL solution cone (dim {dimI})'
L.receipt('gbcheck', exp='gb-hilbert', file=os.path.basename(fn), d=meta['d'], prime=meta['p'], chis=meta['chis'], dim_W=meta['n'],
          n_quadrics=meta['n_quadrics'], rank_quadrics=meta['rank_quadrics'], good_reduction=meta['good_reduction'],
          gb_size=len(polys), max_gb_degree=max(sum(int(x) for x in __import__('re').findall(r'\^(\d+)', pl.split('+')[0].split('-')[0])) or 1 for pl in polys[:1]) if polys else 0,
          dim=dimI, mult=multI, HF=HF, HF_zero_from=k0, msolve_seconds=tgb, seconds=round(time.time() - t0, 1), verdict=verdict,
          msolve_tail=msolve_log[-400:])
print('GBCHECK DONE', verdict, flush=True)
