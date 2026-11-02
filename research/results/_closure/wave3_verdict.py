#!/usr/bin/env python3
"""Verdict for an msolve Groebner output (`-g 1` leading ideal, or `-g 2` reduced basis) of a HOMOGENEOUS ideal I over F_p.

V(I) = {0} over the algebraic closure of F_p  <=>  the leading (initial) ideal of I contains a pure power of every
variable (a monomial ideal is m-primary iff it contains a pure power of each variable, and in(I) and I share the
Hilbert function).  Usage: wave3_verdict.py <input.ms> <output.lead|.gb> <msolve-rc>   -> prints one JSON object.
For `-g 2` files each basis element is printed with its grevlex leading term first ("1*w21^2+..."), so the leading
monomial is the first term with its coefficient stripped.  Adapted from
/tmp/hodge/results/quadric-contact-residual/wave2-excl/msverdict.py (wave 2)."""
import sys, re, os, json


def leading_monomials(txt):
    m = re.search(r'\[(.*)\]', txt, re.S)
    if not m:
        return None
    body = m.group(1).replace('\n', '')
    elems = [s.strip() for s in body.split(',') if s.strip()]
    leads = []
    for e in elems:
        first = e.split('+')[0].strip()           # leading term (coefficients mod p are printed positive)
        first = re.sub(r'^\d+\*', '', first)      # strip coefficient
        first = re.sub(r'^-', '', first)
        leads.append(first)
    return leads


def verdict(ms_file, out_file, rc):
    res = {}
    try:
        vars_ = open(ms_file).readline().strip().split(',')
        res['n_vars'] = len(vars_)
        with open(ms_file) as f:
            f.readline(); res['prime'] = f.readline().strip(); res['n_polys'] = sum(1 for _ in f)
    except Exception as e:
        return {'verdict': 'input-unreadable', 'error': str(e)}
    if rc != 0 or not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        res['verdict'] = 'msolve-failed-or-timeout'; return res
    txt = open(out_file).read()
    res['output_kind'] = 'leading-ideal' if 'Leading ideal' in txt[:200] else ('reduced-GB' if 'Groebner basis' in txt[:200] else 'unknown')
    lead = leading_monomials(txt)
    if lead is None:
        res['verdict'] = 'unparsed'; return res
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
