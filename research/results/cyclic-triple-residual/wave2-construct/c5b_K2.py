# C5b: K2 = {u in R_12^chi : Def_G . u = 0 in R_18^chi}  (annihilator of the invariant deformations in the (2,2)-part of S_A)
import os, sys, json, time, itertools
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[v]='1'
sys.path.insert(0,'/tmp/hodge/results/cyclic-triple-residual'); sys.path.insert(0,'/home/kejace/devel/sextics-work')
exec(open('/tmp/hodge/results/cyclic-triple-residual/wave2-construct/c5_sector_vhs.py').read().split("Def = [e for e in mons(6)")[0])
Def = [e for e in mons(6) if chi(e)==(0,0,0)]
p = 1000003
t0 = time.time()
S12, idx12, g12, M12, r12 = sector_data(12, chiA, p)
S18, idx18, g18, M18, r18 = sector_data(18, chiA, p)
R12, R18 = len(S12)-r12, len(S18)-r18
# quotient coordinates for S18/J18 via rref of the transpose
MT = nmod_mat(M18.ncols(), len(S18), p)
for i in range(len(S18)):
    for j in range(M18.ncols()):
        if M18[i,j] != 0: MT[j,i] = M18[i,j]
MTr, rkT = MT.rref()
piv=[]; 
for i in range(rkT):
    j=0
    while j < len(S18) and MTr[i,j]==0: j+=1
    piv.append(j)
pivset=set(piv); nonpiv=[j for j in range(len(S18)) if j not in pivset]
assert len(nonpiv)==R18, (len(nonpiv), R18)
rows=[[int(MTr[i,j]) for j in range(len(S18))] for i in range(rkT)]
def nf(w):
    for i,pi in enumerate(piv):
        if w[pi]:
            c=w[pi]; row=rows[i]
            for j in range(len(S18)):
                if row[j]: w[j]=(w[j]-c*row[j])%p
    return [w[j] for j in nonpiv]
print('rref J18 done %.1fs' % (time.time()-t0), flush=True)
big = nmod_mat(len(Def)*R18, len(S12), p)
for jcol, b12 in enumerate(S12):
    for k, m in enumerate(Def):
        w=[0]*len(S18); e=tuple(x+y for x,y in zip(b12,m)); w[idx18[e]]=1
        v=nf(w)
        for i,val in enumerate(v):
            if val: big[k*R18+i, jcol]=val
kerdim = len(S12) - big.rank()
K2 = kerdim - (len(S12) - R12)
rec={'p':p,'R12':R12,'R18':R18,'dim_K2_annihilated_by_Def_in_R12':K2,'seconds':round(time.time()-t0,1)}
print(rec, flush=True)
json.dump(rec, open('/tmp/hodge/results/cyclic-triple-residual/wave2-construct/c5b_K2.json','w'), indent=1)
