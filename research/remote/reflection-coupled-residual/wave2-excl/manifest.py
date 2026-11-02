"""Write MANIFEST.json: sha256 + size of every file in the packet (receipts, scripts, exports)."""
import hashlib, os, json, time
D = '/tmp/hodge/results/reflection-coupled-residual/wave2-excl'
man = {}
for root, _, files in os.walk(D):
    for f in files:
        if f == 'MANIFEST.json': continue
        p = os.path.join(root, f)
        h = hashlib.sha256(open(p, 'rb').read()).hexdigest()
        man[os.path.relpath(p, D)] = {'sha256': h, 'bytes': os.path.getsize(p)}
json.dump({'generated': time.strftime('%Y-%m-%dT%H:%M:%S'), 'files': man}, open(f'{D}/MANIFEST.json', 'w'), indent=1)
print(len(man), 'files')
