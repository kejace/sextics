#!/usr/bin/env python3
"""Verdict for an `msolve -g 1` leading-ideal output of a HOMOGENEOUS ideal I over F_p.

V(I) = {0} over the algebraic closure of F_p  <=>  the leading (initial) ideal of I contains a pure power of every
variable (a monomial ideal is m-primary iff it contains a pure power of each variable, and in(I) and I share the
Hilbert function).  Usage: wave3_verdict.py <input.ms> <output.lead> <msolve-rc>   -> prints one JSON object.
Adapted from /tmp/hodge/results/quadric-contact-residual/wave2-excl/msverdict.py (wave 2)."""
import sys, re, os, json


def verdict(ms_file, lead_file, rc):
    res = {}
    try:
        vars_ = open(ms_file).readline().strip().split(',')
        res['n_vars'] = len(vars_)
        with open(ms_file) as f:
            f.readline(); res['prime'] = f.readline().strip(); res['n_polys'] = sum(1 for _ in f)
    except Exception as e:
        return {'verdict': 'input-unreadable', 'error': str(e)}
    if rc != 0 or not os.path.exists(lead_file) or os.path.getsize(lead_file) == 0:
        res['verdict'] = 'msolve-failed-or-timeout'; return res
    txt = open(lead_file).read()
    m = re.search(r'\[(.*)\]', txt, re.S)
    if not m:
        res['verdict'] = 'unparsed'; return res
    lead = [s.strip() for s in m.group(1).replace('\n', '').split(',') if s.strip()]
    res['n_lead'] = len(lead)
    if lead == ['1']:
        res['verdict'] = 'UNIT-IDEAL (impossible for a homogeneous ideal; check input)'; return res
    powers = {}
    for mono in lead:
        mm = re.fullmatch(r'([A-Za-z_]\w*)(?:\^(\d+))?', mono)
        if mm:
            v, k = mm.group(1), int(mm.group(2) or 1)
            powers[v] = min(powers.get(v, 10**9), k)
    missing = [v for v in vars_ if v not in powers]
    res['n_pure_powers'] = len(powers)
    res['max_pure_power'] = max(powers.values()) if powers else None
    res['n_missing_variables'] = len(missing)
    res['missing_variables_head'] = missing[:12]
    res['verdict'] = 'ONLY-ORIGIN (m-primary leading ideal)' if not missing else 'POSITIVE-DIMENSIONAL'
    return res


if __name__ == '__main__':
    print(json.dumps(verdict(sys.argv[1], sys.argv[2], int(sys.argv[3]))))
