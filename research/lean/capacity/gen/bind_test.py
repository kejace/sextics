#!/usr/bin/env python3
"""Test the witness-binding route of the Kronecker library:

  HW_D (exact output of scripts/proposals/check.py:polynomial on a balanced witness string)
     = denW wD            (by rfl)
     = toC tD             (by matchW_sound, decide +kernel)

Usage: python3 bind_test.py --d 8 --kD 300 --name Bind300
"""
from __future__ import annotations
import argparse, ast, json, os, sys, random
from fractions import Fraction
sys.path.insert(0, '/tmp/hodge/lean/capacity/gen'); sys.path.insert(0, '/home/kejace/devel/sextics-work')
sys.path.insert(0, '/home/kejace/devel/sextics/scripts/proposals')
from hodge import core
from hodge.core import NV, VARS
import kron_gen as K
from check import polynomial as translate

CAP = '/tmp/hodge/lean/capacity'


def term_str(e, c):
    mon = '*'.join(f'{VARS[i]}**{k}' if k > 1 else VARS[i] for i, k in enumerate(e) if k)
    c = int(c)
    if abs(c) == 1:
        s = mon if mon else '1'
        return ('-' + s) if c < 0 else s
    s = f'{abs(c)}*{mon}' if mon else f'{abs(c)}'
    return ('-' + s) if c < 0 else s


def witness_string(P, m):
    """Balanced parenthesized sum of the terms in the same (position-sorted) order as kron_gen.pt_lean."""
    items = sorted(((K.pos(e, m), e, c) for e, c in P.items() if c), key=lambda t: t[0])
    def build(lo, hi):
        if hi - lo == 1:
            return term_str(items[lo][1], items[lo][2])
        mid = (lo + hi) // 2
        return f'({build(lo, mid)} + {build(mid, hi)})'
    return build(0, len(items))


def w_ast(src):
    idx = {v: i for i, v in enumerate(VARS)}
    def t(n):
        if isinstance(n, ast.Name):
            return f'(.var {idx[n.id]})'
        if isinstance(n, ast.Constant):
            return f'(.nat {n.value})'
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            return f'(.neg {t(n.operand)})'
        if isinstance(n, ast.BinOp):
            if isinstance(n.op, ast.Pow):
                return f'(.pow {t(n.left)} {n.right.value})'
            op = {ast.Add: 'add', ast.Sub: 'sub', ast.Mult: 'mul'}[type(n.op)]
            return f'(.{op} {t(n.left)} {t(n.right)})'
        raise ValueError(ast.dump(n))
    return t(ast.parse(src, mode='eval').body)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, default=8)
    ap.add_argument('--kD', type=int, default=300)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--name', required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    mons = core.monomials(a.d)
    chosen = rng.sample(mons, min(a.kD, len(mons)))
    D = {}
    for e in chosen:
        c = 0
        while c == 0:
            c = rng.randint(-50, 50)
        D[e] = c
    m = a.d + 1
    s = witness_string(D, m)
    lean_term = translate(s)
    w = w_ast(s)
    t = K.pt_lean(D, m, a.d)
    ns = f'HodgeComputer.Cap.{a.name}'
    L = ['import Reflect.Kron\n', 'set_option maxHeartbeats 4000000\nset_option maxRecDepth 8000\n',
         f'namespace {ns}\nopen HodgeComputer.Kron\n',
         'noncomputable section\n',
         f'def HW_D : HodgeComputer.Poly := {lean_term}\n',
         f'noncomputable def wD : W := {w}\n', f'noncomputable def tD : PT := {t}\n',
         'theorem b1 : HW_D = denW wD := rfl\n',
         'theorem b2 : denW wD = toC tD := matchW_sound _ _ (by decide +kernel)\n',
         'theorem b3 : HW_D = toC tD := b1.trans b2\n',
         f'end\nend {ns}\n']
    src = '\n'.join(L)
    path = os.path.join(CAP, 'src', 'Cap', a.name + '.lean')
    open(path, 'w').write(src)
    json.dump({'D': len(D), 'd': a.d, 'bytes': len(src.encode()), 'route': 'bind', 'name': a.name,
               'witness_bytes': len(s)}, open(os.path.join(CAP, 'src', 'Cap', a.name + '.json'), 'w'))
    print(path, len(D), 'terms', len(src), 'bytes; witness string', len(s), 'bytes')


if __name__ == '__main__':
    main()
