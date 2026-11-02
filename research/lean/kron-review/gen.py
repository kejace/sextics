import ast
def polynomial(source):
    variables = dict(zip(('X0','X1','X2','Y0','Y1','Y2'), range(6)))
    def translate(n):
        if isinstance(n, ast.Name) and n.id in variables:
            return f'(MvPolynomial.X (R := ℂ) ({variables[n.id]} : Fin 6))'
        if isinstance(n, ast.Constant) and type(n.value) is int:
            return f'({n.value} : HodgeComputer.Poly)'
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub): return f'(-{translate(n.operand)})'
        if isinstance(n, ast.BinOp):
            if isinstance(n.op, ast.Pow):
                return f'({translate(n.left)} ^ ({n.right.value} : ℕ))'
            op = {ast.Add:'+', ast.Sub:'-', ast.Mult:'*'}.get(type(n.op))
            if op: return f'({translate(n.left)} {op} {translate(n.right)})'
        raise ValueError('Unsupported')
    return translate(ast.parse(source, mode='eval').body)
s = "-5*X0**6 + (3*X0**4*X1**2 + -X1**5*Y2)"
print(polynomial(s))
