import json
rows=[json.loads(l) for l in open('/tmp/hodge/results/quadratic-cubic-pullback/log.jsonl')]
print(f"{'family':18s} {'d':>2s} {'E1':>5s} {'cons':>5s} {'rank':>5s} {'null':>4s} {'sec':>6s}  D / note")
for r in rows:
    D=r['D']; D=D if len(D)<40 else D[:37]+'...'
    print(f"{r['family']:18s} {str(r['d']):>2s} {str(r['e1']):>5s} {str(r['consistent']):>5s} {str(r['rank']):>5s} {str(r['nullity']):>4s} {str(r['seconds']):>6s}  {D}  {r.get('note','')[:60]}")
print(len(rows),'rows')
