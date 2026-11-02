"""Add a random affine linear equation l(a) - 1 = 0 to a homogeneous msolve system (so that msolve solves the
zero-dimensional dehomogenized system).  python3 wave2/dehom.py in.ms out.ms [seed]"""
import sys, random
src, dst = sys.argv[1], sys.argv[2]
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
lines = open(src).read().split('\n')
names = lines[0].split(','); p = int(lines[1])
rnd = random.Random(seed)
lin = '+'.join(f'{rnd.randrange(1, p)}*{n}' for n in names) + f'+{p-1}'
body = '\n'.join(lines[2:]).rstrip().rstrip(',')
open(dst, 'w').write(lines[0] + '\n' + lines[1] + '\n' + body + ',\n' + lin + '\n')
print('wrote', dst, 'with', len(names), 'vars')
