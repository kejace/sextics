import sys, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import core, fixed_d
LOG='/tmp/hodge/results/polarized-k3-remainder/wave4/log.jsonl'
def log(**kw):
    kw['ts']=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()); open(LOG,'a').write(json.dumps(kw, default=str)+'\n')
F, A, q, rec = core.load_instance('polarized-k3-remainder')
Amix=core.parse(json.load(open('/tmp/hodge/results/polarized-k3-remainder/wave4/amix_pieces.json'))['A_mix'])
X=[core.parse(v) for v in ['X0','X1','X2','Y0','Y1','Y2']]
U=[core.add(X[i],X[i+3]) for i in range(3)]; V=[core.sub(X[i],X[i+3]) for i in range(3)]
def f_of(Z):
    z0,z1,z2=Z; m=core.mul; pw=core.power; from fractions import Fraction as Fr
    terms=[(1,pw(z0,6)),(Fr(5,7),m(pw(z0,5),z1)),(Fr(10,7),m(pw(z1,5),z2)),(Fr(15,7),m(pw(z2,5),z0)),(Fr(5,7),m(m(pw(z0,2),pw(z1,3)),z2)),(1,pw(z1,6)),(1,pw(z2,6))]
    out={}
    for cf,P in terms: out=core.add(out, core.scale(P,cf))
    return out
fU=f_of(U); fV=f_of(V)
assert core.sub(core.scale(core.add(fU,fV),0.5) if False else core.add(fU,fV), core.scale(F,2))=={}, 'F != (f(U)+f(V))/2'
cands={'f(V)':fV,'f(U)':fU,'f(U)-f(V)':core.sub(fU,fV)}
# control: joining-sector target U0^3 V0^3 with D = f(V) (wave-2 B2 says strict certificate exists)
ctrl=core.mul(core.power(U[0],3),core.power(V[0],3))
for nm,D in cands.items():
    for tn,Tg in [('A_mix',Amix),('control U0^3V0^3',ctrl)]:
        if tn.startswith('control') and nm!='f(V)': continue
        t0=time.time()
        r=fixed_d.solve_fixed_D_fast(F,Tg,1,D,want_exact=True,verbose=False)
        info={k:(v if k not in ('N','Q','D') else '...') for k,v in (r.items() if isinstance(r,dict) else [('raw',str(r)[:200])])}
        ok=None
        if isinstance(r,dict) and r.get('N') is not None:
            try:
                ok=core.verify_witness(F,Tg,1,D,r['N'],r.get('Q'))
            except Exception as e: ok=('verify-error',str(e)[:100])
        print(nm,'|',tn,'|',{k:info[k] for k in list(info)[:8]},'| verified:',ok, '| %.1fs'%(time.time()-t0), flush=True)
        log(what='e2-fixed-D-degree6-calibration', D=nm, target=tn, result={'solver':{k:str(v)[:120] for k,v in info.items()},'verified':str(ok)}, prime=core.PRIMES[0])
