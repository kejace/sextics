"""Run the official sandboxed verifier (scripts/proposals/sandbox.verify) on a folder of proposal files."""
import json, sys, tempfile
from pathlib import Path
import yaml
ROOT = Path('/home/kejace/devel/sextics')
sys.path.insert(0, str(ROOT / 'scripts/proposals'))
from check import module_order
from sandbox import verify

def run(instance_id, folder):
    folder = Path(folder)
    files = {p.name: p.read_bytes() for p in folder.iterdir() if p.is_file() and p.name in ('proof.json','solution.md','witness.json') or p.suffix == '.lean'}
    if 'solution.md' not in files:
        files['solution.md'] = b'# placeholder explanation\n'
    manifest = json.loads(files['proof.json'])
    instance = yaml.safe_load((ROOT / 'data' / instance_id / 'instance.rec').read_text())
    policy = json.loads((ROOT / '.github/solve-policy.json').read_text())
    order, imports = module_order(files)
    logs = Path(tempfile.mkdtemp(prefix='hodge-proof-check-'))
    print(f'Verification logs: {logs}', flush=True)
    print(json.dumps(verify(ROOT, files, instance, manifest, order, imports, policy, logs), indent=2))

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])
