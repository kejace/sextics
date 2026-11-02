"""Shared helpers for the wave2-excl packet (diagonal-permutation)."""
import sys, json, time, os
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core
LANE = '/tmp/hodge/results/diagonal-permutation'
PACKET = os.path.join(LANE, 'wave2-excl')
LOG = os.path.join(LANE, 'wave2-excl.log.jsonl')
IID = 'diagonal-permutation'

def log(entry):
    entry.setdefault('lens', 'excl'); entry.setdefault('engine', 2)
    entry.setdefault('ts', time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(LOG, 'a') as f:
        f.write(json.dumps(entry, default=str) + '\n')
    print(json.dumps(entry, default=str)[:800], flush=True)

def receipt(name, obj):
    with open(os.path.join(PACKET, name), 'w') as f:
        json.dump(obj, f, indent=1, default=str)

def singular_poly(P):
    """Polynomial dict -> Singular expression string (variables X0..Y2, '^' powers)."""
    terms = []
    names = ['X0', 'X1', 'X2', 'Y0', 'Y1', 'Y2']
    for e, c in sorted(P.items()):
        c = core.Fraction(c)
        mon = '*'.join(f'{names[i]}^{k}' if k > 1 else names[i] for i, k in enumerate(e) if k)
        cs = f'({c.numerator})' if c.denominator == 1 else f'({c.numerator}/{c.denominator})'
        terms.append(cs + ('*' + mon if mon else ''))
    return '+'.join(terms) if terms else '0'
