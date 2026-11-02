#!/usr/bin/env bash
# Re-apply the VM-side engine patches (idempotent).  Run after a push from another machine reverted hodge/*.py.
cd /home/kejace/devel/sextics-work || exit 1
for s in tools/patch_fixed_d_rankfix.py tools/patch_fixed_d_sparse.py tools/patch_fixed_d_fastauto.py tools/patch_sweep_schema.py tools/patch_rjob_local.py; do python3 "$s"; done
HODGE_QUIET=1 python3 -c "
import sys; sys.path.insert(0, '.')
from hodge import fixed_d, rjob
print('engine', fixed_d.ENGINE_VERSION, 'solve_sparse', hasattr(fixed_d, 'solve_sparse'), 'rjob.is_local', hasattr(rjob, 'is_local'))"
