"""Sampled replay of redsys_lift.py (3 exact degree-1 polynomial Farkas vectors on the reduced trivial-sector system).
Outputs are redirected to /tmp/hodge/referee/cycle-q2/packet_shadow/; the agent's log is not touched."""
import sys, os
sys.path.insert(0, '/tmp/hodge/results/cycle-q2/wave2-excl')
import w2lib
w2lib.PACKET = '/tmp/hodge/referee/cycle-q2/packet_shadow'
w2lib.log = lambda **kw: None
sys.argv = ['redsys_lift.py', '--nvec', '3', '--tag', 'referee_red_trivial6']
import runpy
runpy.run_path('/tmp/hodge/results/cycle-q2/wave2-excl/redsys_lift.py', run_name='__main__')
