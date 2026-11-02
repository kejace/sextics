import json, time
done=set()
for l in open('/tmp/hodge/notes/wave2-verification.log.jsonl'):
    try: e=json.loads(l)
    except: continue
    if e.get('check','').startswith('qcr_d7_pair_'): done.add(e['check'])
ts=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()); n=0
for l in open('qcr7_results.jsonl'):
    r=json.loads(l); name='qcr_d7_pair_'+r['pair']
    if name in done: continue
    a,b=r['pair'].split('_'); diff=sum(x!=y for x,y in zip(a,b))
    kind={1:'edge',2:'face diagonal',3:'space diagonal'}[diff]
    e={"ts":ts,"check":name,"lane":"quadric-contact-residual","kind":"relaxation","claim":f"excl items 3-4: degree-7 denominators supported on the parity-class pair {{{a},{b}}} ({kind}) are excluded (msolve -g 1 m-primary, two primes)","method":"own code qcr7.py at 998244353: D8* recomputed, Euler slices, Hessian relaxation mod J_8 (+D8*), msolve -g 1 leading ideal","prime":998244353,"expected":{"m_primary":True},"got":{k:r.get(k) for k in ('n_vars','independent_quadrics','m_primary','msolve','seconds','max_pure_power')},"agree":bool(r.get('m_primary')) if r.get('msolve')=='ok' else None,"ok":r.get('msolve')=='ok',"error":None if r.get('msolve')=='ok' else r.get('msolve')}
    open('/tmp/hodge/notes/wave2-verification.log.jsonl','a').write(json.dumps(e)+'\n'); n+=1
print('synced',n)
