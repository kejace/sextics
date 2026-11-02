"""Shared helpers for the wave2 construction lens (quadric-contact-residual)."""
import sys, json, time
sys.path.insert(0,'/home/kejace/devel/sextics-work')
import sympy as sp
from hodge import core
X0,X1,X2,Y0,Y1,Y2 = V6 = sp.symbols('X0 X1 X2 Y0 Y1 Y2')
F,A,q,rec = core.load_instance('quadric-contact-residual')
Fs = sp.sympify(core.to_str(F)); As = sp.sympify(core.to_str(A))
LOG='/tmp/hodge/results/quadric-contact-residual/wave2-construct.log.jsonl'
def log(**kw):
    kw.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S')); kw.setdefault('lens','construct'); kw.setdefault('engine',2)
    with open(LOG,'a') as f: f.write(json.dumps(kw, default=str)+'\n')
