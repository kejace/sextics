"""Full fixed-D test (engine v2, compressed + exact lift) of the unique non-Jacobian degree-8 E1-kernel element D8,
and of its (01)-image and (01)-symmetrisations.  Also reports the served target dimension."""
import sys, json, time; sys.path.insert(0, '/tmp/hodge/results/diagonal-permutation')
from hlib import *
from hodge.core import add, sub, scale, parse, to_str
D8 = parse(json.load(open('/tmp/hodge/results/diagonal-permutation/D8_e1kernel.json'))['D'])
im = json.load(open('/tmp/hodge/results/diagonal-permutation/D8_images.json'))
D01 = parse(im['D01'])
tests = [('D8 = unique non-Jacobian E1-kernel element (deg 8, canonical rep mod J_8^+)', D8, 'E1kernel-d8')]
for name, D, fam_ in tests:
    e = test_D_image(name, fam_, D=D)
    print(name, '-> d', e['d'], 'e1', e['e1'], 'cons', e.get('consistent'), 'served', e.get('served_dim'), e['seconds'], 's', flush=True)
# (01)-symmetrisations only if they still pass E1 (they do iff (01)D8 = +-D8 mod J)
for name, D in [('D8 + (01)D8', add(D8, D01)), ('D8 - (01)D8', sub(D8, D01))]:
    if not D: print(name, 'is zero'); continue
    e = test_D_image(name, 'E1kernel-d8-sym', D=D)
    print(name, '-> d', e['d'], 'e1', e['e1'], 'cons', e.get('consistent'), 'served', e.get('served_dim'), e['seconds'], 's', flush=True)
