#!/usr/bin/env python3
"""(i) covers creating a reflection, book §2.7 shape.  A §2.7 cover needs F = f(rest) + y^6 + l1^6 + l2^6
+ y*(l1*a1 + l2*a2) with y,l1,l2 coordinates (Fermat-curve automorphisms force coordinates, see notes);
then T^6 = l1^6+l2^6, S = y^2+T^2, P = yT and y<->T is the created reflection.  We enumerate every
ordered choice (y; {l1,l2}) among the six coordinates and list the monomials of F that violate the shape:
(a) monomials containing y with y-degree != 1 or 6, (b) monomials containing both l1,l2 (other than through
l1^6+l2^6), (c) monomials with l-degree not in {0,1,6}, (d) y-linear monomials not linear in exactly one l."""
import itertools, json
F = {(6,0,0,0,0,0):1,(0,6,0,0,0,0):1,(0,0,6,0,0,0):1,(0,0,0,6,0,0):1,(0,0,0,0,6,0):1,(0,0,0,0,0,6):1,
     (3,1,0,1,1,0):0.5,(0,3,1,0,1,1):0.5,(1,0,3,1,0,1):0.5}
names=["X0","X1","X2","Y0","Y1","Y2"]
def mono_str(m): return "*".join(f"{names[i]}^{e}" if e>1 else names[i] for i,e in enumerate(m) if e)
out=[]
for y in range(6):
    for l1,l2 in itertools.combinations([i for i in range(6) if i!=y],2):
        bad=[]
        for m in F:
            ey,e1,e2=m[y],m[l1],m[l2]
            if ey==6 and e1==e2==0: continue
            if (e1,e2) in ((6,0),(0,6)) and ey==0: continue
            if ey==0 and e1==0 and e2==0: continue          # f(rest)
            if ey==1 and ((e1,e2)==(1,0) or (e1,e2)==(0,1)): continue   # y*l*a(rest): OK
            bad.append(mono_str(m))
        out.append({"y":names[y],"l":[names[l1],names[l2]],"violating":bad})
ok=[o for o in out if not o["violating"]]
print("choices:",len(out),"admissible:",len(ok))
# minimal violation count
out.sort(key=lambda o:len(o["violating"]))
for o in out[:8]: print(o["y"],o["l"],"violations:",o["violating"])
json.dump(out,open("w2c_cover.json","w"),indent=0)
