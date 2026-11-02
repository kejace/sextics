#!/usr/bin/env python3
"""Witness -> Lean package -> local compile check (only run after an exact, verified witness exists).
  python3 w2_make_lean.py [witness-candidate.json] [--no-check]
Writes /tmp/hodge/lean/off-paired/{Witness.lean,Solution.lean,proof.json,witness.json} and runs hodge.leancheck."""
import os, sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from hodge import core, leangen, jacobian
src = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('--') else os.path.join(os.path.dirname(HERE), 'witness-candidate.json')
w = json.load(open(src))
F, A, q, rec = core.load_instance('off-paired')
D = core.parse(w['D']); N = [core.parse(s) for s in w['N']]; Q = core.parse(w['Q'])
ok, info = core.verify_witness(F, A, q, D, N, Q)
print('exact identity verified:', ok, info)
assert ok, 'witness does not verify'
divF, _ = jacobian.divisible_by(F, D)
assert D and not divF, 'D must be nonzero and not divisible by F (F irreducible => coprime)'
print('D nonzero and F does not divide D: coprimality OK (F is irreducible)')
out = '/tmp/hodge/lean/off-paired'; os.makedirs(out, exist_ok=True)
files = leangen.generate('off-paired', {'d': w['d'], 'D': w['D'], 'N': w['N'], 'Q': w['Q']}, name='OffPaired')
open(f'{out}/Witness.lean', 'w').write(files['Witness.lean']); open(f'{out}/Solution.lean', 'w').write(files['Solution.lean'])
json.dump({'schemaVersion': 2, 'instanceId': 'off-paired', 'theorem': 'Submission.certificate', 'witnessTheorem': 'Submission.witness'},
          open(f'{out}/proof.json', 'w'), indent=2)
json.dump(leangen.witness_json(w['d'], D, N, Q), open(f'{out}/witness.json', 'w'))
print('Lean package written to', out, 'meta', files['meta'])
if '--no-check' not in sys.argv:
    from hodge import leancheck
    t0 = time.time()
    res = leancheck.run('off-paired', out, jobs=2)
    print('leancheck result:', res, f'({time.time()-t0:.0f}s)')
