"""Referee: exact rational identity A*K8 = sum N_i dF/dz_i (python Fractions), K8 not in J mod p, K8 in E1 kernel mod p."""
import sys, json, time
sys.path.insert(0, '/tmp/hodge/results/reflection-coupled-residual/wave2-excl/referee')
import reflib as R
p = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000021
PK = '/tmp/hodge/results/reflection-coupled-residual/wave2-excl'
e = json.load(open(f'{PK}/K8_E1_exact.json'))
K8 = R.parse(e['D']); N = [R.parse(s) for s in e['N']]
F, A, q = R.load_instance()
t0 = time.time()
lhs = A * K8
rhs = R.Poly()
for i in range(6): rhs = rhs + N[i] * F.deriv(i)
diff = lhs - rhs
print('deg K8', K8.degree(), 'terms', len(K8.t), 'chars', K8.chars(), 'deg N_i', [n.degree() for n in N], 'N terms', [len(n.t) for n in N])
print('EXACT A*K8 - sum N_i F_i == 0 over Q:', diff.is_zero(), f'({time.time()-t0:.1f}s)')
# consistency with the other copies
k8b = json.load(open(f'{PK}/../kernel_d8_X03Y03.json'))[0]; k8c = json.load(open(f'{PK}/e1reps_d8_chi00_Q.json'))[0]
print('same D string in kernel_d8_X03Y03.json / e1reps_d8_chi00_Q.json:', (R.parse(k8b) - K8).is_zero(), (R.parse(k8c) - K8).is_zero())
S = R.Setup(p)
R8 = R.Reducer(S, 8, (0, 0)); print('good reduction deg 8 (0,0):', R8.good, R8.rank, R8.dimRJ, R8.fermat)
print('K8 in J_8 mod p?', R8.rank_mod_J([K8.modp(p)]) == 0)
R14 = R.Reducer(S, 14, (0, 1)); print('good reduction deg 14 (0,1):', R14.good)
print('A*K8 in J_14 mod p?', R14.rank_mod_J([(A * K8).modp(p)]) == 0)
R.receipt(f'{PK}/referee/ref_k8.jsonl', dict(exp='ref-k8', p=p, exact_identity=diff.is_zero(), K8_terms=len(K8.t), K8_chars=[list(c) for c in K8.chars()],
          K8_in_J_modp=(R8.rank_mod_J([K8.modp(p)]) == 0), good8=R8.good, good14=R14.good, AK8_in_J_modp=(R14.rank_mod_J([(A * K8).modp(p)]) == 0)))
