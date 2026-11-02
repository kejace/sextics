"""Convert an msolve -g 1 leading-ideal output into a Singular script computing dim / vdim / Hilbert series of the
monomial ideal (same as the ideal's, since leading ideals share the Hilbert function).  Usage: lead2sing.py lead.txt n p"""
import sys, re
lead, n, p = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
txt = open(lead).read()
mons = [m.strip() for m in re.split(r',\s*', txt.replace('[', '').replace(']', '').replace(':', '')) if m.strip()]
out = sys.argv[4] if len(sys.argv) > 4 else lead + '.sing'
with open(out, 'w') as f:
    f.write(f'ring R = {p}, (a1..a{n}), dp;\nideal I = ' + ',\n'.join(mons) + ';\n')
    f.write('print("size = " + string(size(I)));\nprint("dim = " + string(dim(I)));\nprint("vdim = " + string(vdim(I)));\n')
    f.write('intvec h = hilb(I, 1); print("hilb1 = " + string(h));\nquit;\n')
print('wrote', out, 'with', len(mons), 'leading monomials')
