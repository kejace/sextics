import sys, json
sys.path.insert(0, '/tmp/hodge/results/short-power-presentation')
from w2lib import *
F, A, q, rec = load()
fam = json.load(open('families_d8_p0_w2c5.json'))['families']['XM+YL']
reps = [el for el in fam['rational'] if el.get('nonjacobian')]
out = []
for i, el in enumerate(reps):
    cs = el['coords']
    keys = list(cs)
    nXM = sum(1 for k in keys if k.startswith('XM')); nYL = sum(1 for k in keys if k.startswith('YL'))
    D = core.parse(el['D'])
    okE1 = "by construction"
    hx = max(abs(Fraction(v).numerator).bit_length() for v in cs.values()); hd = max(Fraction(v).denominator.bit_length() for v in cs.values())
    Md = sorted({int(k.split()[4]) for k in keys if k.startswith('XM')}); Ld = sorted({int(k.split()[4]) for k in keys if k.startswith('YL')})
    print(f'rep {i}: #XM coords {nXM} (M-degrees {Md}), #YL coords {nYL} (L-degrees {Ld}), E1 {okE1}, heights num {hx} bits / den {hd} bits, terms {len(D)}', flush=True)
    out.append({'D': el['D'], 'nXM': nXM, 'nYL': nYL, 'E1': bool(okE1)})
json.dump(out, open('K8even_twocone_reps.json', 'w'))
print('saved K8even_twocone_reps.json', flush=True)
