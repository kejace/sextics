"""Verify witness-candidate.json exactly, check coprimality, generate the Lean package and compile it locally.
  python3 wave2/package_witness.py [path/to/witness-candidate.json] [--no-compile]
"""
import sys, os, json, subprocess
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, leangen, jacobian
IID = 'rank-one-biquartic'
path = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('--') else f'/tmp/hodge/results/{IID}/witness-candidate.json'
w = json.load(open(path))
F, A, q, rec = core.load_instance(IID)
D = core.parse(w['D']); N = [core.parse(s) for s in w['N']]; Q = core.parse(w['Q'])
ok, info = core.verify_witness(F, A, q, D, N, Q)
print('exact verify:', ok, info)
assert ok
divF, _ = jacobian.divisible_by(F, D)
print('F divides D (mod p)?', divF)
assert not divF
assert core.is_homogeneous(D) and core.parity_of(D) is not None
out = f'/tmp/hodge/lean/{IID}'
os.makedirs(out, exist_ok=True)
files = leangen.generate(IID, {'d': w['d'], 'D': w['D'], 'N': w['N'], 'Q': w['Q']}, name='RankOneBiquartic')
open(f'{out}/Witness.lean', 'w').write(files['Witness.lean'])
open(f'{out}/Solution.lean', 'w').write(files['Solution.lean'])
json.dump({'schemaVersion': 2, 'instanceId': IID, 'theorem': 'Submission.certificate', 'witnessTheorem': 'Submission.witness'}, open(f'{out}/proof.json', 'w'), indent=2)
json.dump({'d': int(w['d']), 'D': w['D'], 'N': w['N'], 'Q': w['Q']}, open(f'{out}/witness.json', 'w'), indent=1)
print('Lean package written to', out, files['meta'])
if '--no-compile' not in sys.argv:
    subprocess.run([sys.executable, '-m', 'hodge.leancheck', '--instance', IID, '--dir', out], cwd='/home/kejace/devel/sextics-work')
