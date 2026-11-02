import sys, runpy
sys.path.insert(0, '/tmp/hodge/results/four-plus-two-residual/wave2-excl')
import w2lib as W
W.REMOTE_OUT = '/tmp/hodge/work/referee42/out/d6jac2-replay'   # keep the agent's files untouched
import os; os.makedirs(W.REMOTE_OUT, exist_ok=True)
fam = sys.argv[1]
sys.argv = ['w2_d6jac2.py', 'referee-d6jac2-replay', fam]
runpy.run_path('/tmp/hodge/results/four-plus-two-residual/wave2-excl/w2_d6jac2.py', run_name='__main__')
