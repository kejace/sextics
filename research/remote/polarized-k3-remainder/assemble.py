import subprocess, os
os.chdir('/tmp/hodge/results/polarized-k3-remainder')
tables = subprocess.run(['/usr/bin/python3', 'mktables.py'], capture_output=True, text=True).stdout
nexts = """1. **Decide ρ and test loophole (B)** — certified periods of `w² = f` (Lairez–Sertöz `PeriodSuite` or
   Sertöz's K3 period code), LLL for NS, then test whether the period vector of η = [z0²z1²z2²] lies in
   the NS_prim span. CPU: hours on one core; needs software not present on this machine (Sage/Magma +
   PeriodSuite). Outcome: either closes the instance conditionally (ρ = 1 ⟹ equivalent to the open GHC
   case) or exhibits the only cheap route. Alternative, cheaper first step: Picard rank bound at a second
   prime via p-adic point counting (Costa–Harvey controlledreduction) — if ρ(F̄_p) = 2 at some p with a
   discriminant class different from the p = 3 one, van Luijk gives ρ = 1 unconditionally.
2. **Degree-8 exclusion for irreducible D** using the sum principle (§1.3 ii) and the dichotomy (§4.3):
   restrict the discussion's 315-parameter Koszul relaxation to the complement of the reducible locus
   {f_U h, (f_U − f_V) h, P_i·cubic}, and add the 18 further K3 targets as simultaneous right-hand sides
   (same 315 denominator parameters). msolve / Singular modStd, hours–days, several GB. This would raise the
   unconditional bound to degree ≥ 9.
3. **Served-space diagnostics for every future candidate** (`leanimage.py`, ≈ 10 min at d = 8, ≈ 1 h at
   d = 10 with ≈ 6 GB): by the dichotomy a single nonzero K3 projection is decisive, so this is the right
   first test for any D anybody proposes, before any Lean work.
4. **Degree-10 structural candidates** (polar pairing Σ ∂_if(U)∂_if(V), twisted pairings, products of
   the two excluded sextics f(U)−f(V) and quartics are reducible hence dead): ≈ 20k × 20k compressed
   solves, 30–60 min and ≈ 5 GB each on an idle machine (`candidates_d10.json` is ready); expected
   inconsistent by §1; low priority.
5. **Discussion post** (not done in this session, per rules): the conditional equivalence theorem (§1.2),
   the dichotomy and sum principle (§1.3), the maximal-rank table through d = 12 (§1.5) and the negative
   tables (§2) sharpen 1337isnot1337's "degree ≥ 8" to "an irreducible factor of degree ≥ 8, and in fact
   a divisor whose resolution carries so(T)(1)". Credit: 1337isnot1337 (degree ≤ 7 exclusion, End_Hdg(T) = Q,
   Picard bound), the book's Ch. 3A sector decomposition, DiamondJM's degree-6 sweeps.
"""
src = open('note_draft.md').read().replace('<<TABLES>>', tables).replace('<<NEXTSTEPS>>', nexts)
open('/tmp/hodge/notes/attack1-polarized-k3-remainder.md', 'w').write(src)
print('assembled', len(src), 'bytes')
