import re, sys
src = open('/tmp/hodge/lean/capacity/src/Reflect/Kron.lean').read()
head, sep, body = src.partition('namespace HodgeComputer.Kron\n')
assert sep, 'namespace line not found'
assert body.rstrip().endswith('end HodgeComputer.Kron'), body[-80:]
body = body.rstrip()[: -len('end HodgeComputer.Kron')].rstrip('\n') + '\n'
header = '''import HodgeComputer.Verify

/-!
# Kronecker-substitution checker for explicit polynomial identities

Large exact identities between explicit polynomials in `MvPolynomial (Fin 6) ℂ`, such as the
cleared certificate equations of a witness, are verified by the kernel from one big-integer
comparison instead of by `ring` normalization.

An explicit polynomial with integer coefficients is a balanced binary tree of monomials (`PT`)
whose leaves are sorted by Kronecker position. An identity `L = R` between expressions (`E`)
built from such trees by products, sums, differences, natural scalars, powers, and partial
derivatives is certified by three decidable checks:

* `homE`, `degE`: both sides are homogeneous of the same degree `n < m`;
* `l1E`: the ℓ¹ norms bound every coefficient of `L - R` strictly below `2 ^ W`;
* `psiEq`: the Kronecker values agree, where `psiE m W` evaluates at `X i ↦ 2 ^ (W * m ^ i)`
  for `i < 5` and at `X 5 ↦ 1`.

`sound` derives `denZ L = denZ R` over `ℤ` from these checks by the uniqueness of bounded signed
digits in a mixed-radix expansion (`eq_zero_of_eval_kpt`); `soundC` and `soundN` transfer the
identity to `ℂ`. The Kronecker values are computed with kernel-accelerated `Nat` shifts and
multiplications on pairs of natural numbers (`psiEN`); the divide-and-conquer evaluation
`psiRelN` keeps every intermediate number proportional to the span of positions in its subtree.
`sound4` is the variant sliced by the degree in `X 4`, for high-degree witnesses whose unsliced
Kronecker numbers would exceed the verifier's memory.

`W`, `denW`, and `matchW_sound` bind the verifier's generated `HodgeWitness` terms to the trees.
Every check is discharged by `decide +kernel`, and the soundness theorems depend only on the
standard axioms `propext`, `Classical.choice`, and `Quot.sound`. The directory `README.md`
describes the application pattern, the verification budget, and pitfalls.

Developed in the Lean capacity study prepared with Claude Code assistance for `kejace`.
-/

set_option autoImplicit false

namespace HodgeComputer.Kronecker
'''
out = header + body + '\nend HodgeComputer.Kronecker\n'
# collapse 3+ blank lines to a single blank line
out = re.sub(r'\n{3,}', '\n\n', out)
open('/home/kejace/devel/sextics-kron/formal/lean/HodgeComputer/Reflection/Kronecker.lean', 'w').write(out)
print(len(out.splitlines()), 'lines')
