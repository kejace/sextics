import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
DATA = '/home/kejace/devel/sextics/data'
SCRATCH = '/tmp/hodge/results/_tests'
SMALL = ['collision-h3', 'fermat-pair-rational-d4-strict', 'three-block-rational-d4-strict']

def scratch_dir(prefix):
    import tempfile
    os.makedirs(SCRATCH, exist_ok=True)
    return tempfile.mkdtemp(prefix=prefix, dir=SCRATCH)
