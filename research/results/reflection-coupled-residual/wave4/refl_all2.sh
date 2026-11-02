#!/bin/bash
# wave4 Exp B (lambda-eliminated): Fermat-5 control, six coordinate pencils, full X. p = 2^31-1 (msolve < 2^31).
cd /tmp/hodge/results/reflection-coupled-residual/wave4
ulimit -v 6000000
export MSOLVE_THREADS=${MSOLVE_THREADS:-3}
P=2147483647
FEXPR='X0**6 + 5*X0**3*X1*Y0*Y1/7 + 5*X0*X1*Y0**3*Y1/7 + 5*X0*X1*Y0**2*Y2**2/7 + 10*X0*X1*Y1**2*Y2**2/7 + 15*X0*X2**3*Y0*Y1/7 + X1**6 + 10*X1**3*X2*Y0*Y1/7 + 5*X1*X2*Y0*Y1**3/7 + X2**6 + Y0**6 + Y1**6 + Y2**6'
run() { name=$1; shift; echo "=== $name $(date -u)"; timeout 7200 python3 reflsearch2.py "$name" $P 7 "$@" > "$name.out" 2>&1; grep -v "^  reflection" "$name.out" | tail -1 | cut -c1-600;
  python3 -c "
import json,sys
try:
    S=json.load(open('$name.json'))
    rec={'t':__import__('time').strftime('%Y-%m-%dT%H:%MZ',__import__('time').gmtime()),'what':'reflsearch2 (lambda-eliminated reflection search)','name':'$name','p':S['p'],'n':S['n'],'subst':S.get('subst'),'n_eqs':S['n_eqs'],'dim':S.get('dim'),'degree':S.get('degree'),'factor_degrees':S.get('factor_degrees'),'genuine_Fp':S.get('genuine_Fp_reflections'),'spurious_Fp':S.get('spurious_Fp'),'msolve_s':S.get('msolve_seconds'),'file':'wave4/$name.json'}
    open('log.jsonl','a').write(json.dumps(rec)+'\n')
except Exception as e: print('log fail',e)
"; }
run ctrlF5 'X0**6+X1**6+X2**6+Y0**6+Y1**6'
run fib_Y0Y1 "$FEXPR" 'Y0=l*Y1'
run fib_Y0Y2 "$FEXPR" 'Y0=l*Y2'
run fib_Y1Y2 "$FEXPR" 'Y1=l*Y2'
run fib_X0X1 "$FEXPR" 'X0=l*X1'
run fib_X0X2 "$FEXPR" 'X0=l*X2'
run fib_X1X2 "$FEXPR" 'X1=l*X2'
run full_X "$FEXPR"
echo "=== ALL DONE $(date -u)"
