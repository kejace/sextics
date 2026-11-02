"""Referee GB check: msolve grevlex GB mod p of an exported quadratic system, then Singular dim / vdim / Hilbert function
(kbase sizes) with the msolve basis declared as a standard basis. Usage: python3 ref_gb.py <file.ms> <threads> [timeout_s] [maxdeg]"""
import sys, json, time, subprocess, os, re
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee')
import reflib as R
fn = sys.argv[1]; thr = int(sys.argv[2]); tmo = int(sys.argv[3]) if len(sys.argv) > 3 else 3600
maxdeg = int(sys.argv[4]) if len(sys.argv) > 4 else 10
OUT = '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee/ref_gb.jsonl'
meta = json.load(open(fn.replace('.ms', '.json')))
t0 = time.time(); gbf = fn.replace('.ms', '.gb')
try:
    r = subprocess.run(['msolve', '-f', fn, '-o', gbf, '-g', '2', '-t', str(thr), '-v', '1'], capture_output=True, text=True, timeout=tmo)
except subprocess.TimeoutExpired:
    R.receipt(OUT, dict(exp='ref-gb', file=os.path.basename(fn), d=meta['d'], p=meta['p'], chis=meta['chis'], dim_W=meta['dim_W'], result='msolve TIMEOUT', timeout_s=tmo))
    sys.exit(0)
tgb = round(time.time() - t0, 1)
log = r.stdout + r.stderr
m = re.search(r'Dimension of quotient:\s*(\d+)', log); msolve_quot = int(m.group(1)) if m else None
vars_ = open(fn).readline().strip().split(',')
txt = open(gbf).read(); body = txt[txt.find('[') + 1: txt.rfind(']')]
polys = [t.strip() for t in body.replace('\n', '').split(',') if t.strip()]
sing = f"""ring r = {meta['p']}, ({','.join(vars_)}), dp;
ideal I = {','.join(polys)};
attrib(I,"isSB",1);
print("DIM"); print(dim(I)); print("MULT"); print(vdim(I)); print("HF");
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
k0 = next((k for k, v in enumerate(HF) if v == 0), None) if HF and all(isinstance(v, int) for v in HF) else None
R.receipt(OUT, dict(exp='ref-gb', file=os.path.basename(fn), d=meta['d'], p=meta['p'], chis=meta['chis'], dim_W=meta['dim_W'],
          n_quadrics=meta['n_quadrics'], rank_quadrics=meta['rank_quadrics'], good_reduction=meta['good_reduction'], gb_size=len(polys),
          msolve_quotient_dim=msolve_quot, singular_dim=dimI, singular_vdim=multI, HF=HF, HF_zero_from=k0, msolve_seconds=tgb,
          seconds=round(time.time() - t0, 1), verdict=('zero-dimensional cone: V(I)={0}' if dimI == 0 else f'dim {dimI}')))
