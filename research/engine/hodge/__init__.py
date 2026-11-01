"""hodge: scratch engine for the hodge.computer sextic challenge (see README.md)."""
import os as _os
# BLAS threading hurts badly on a loaded machine and the dense solves here are small;
# default to single-threaded BLAS unless the caller set the variables (must happen before numpy import).
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')
