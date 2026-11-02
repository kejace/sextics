#!/bin/bash
cd /tmp/hodge/results/symmetric-sextic-residual
python3 w2x_relax_job.py --kind full --d 6 --par mixed --prime 2147483647 --threads 16 --timeout 3600 --singular
python3 w2x_relax_job.py --kind full --d 7 --par 0 --prime 2147483647 --threads 32 --timeout 14400 --singular
python3 w2x_relax_job.py --kind full --d 7 --par 1 --prime 2147483647 --threads 32 --timeout 14400 --singular
python3 w2x_relax_job.py --kind full --d 6 --par mixed --prime 2147483629 --threads 16 --timeout 3600
python3 w2x_relax_job.py --kind full --d 7 --par 0 --prime 2147483629 --threads 32 --timeout 14400
python3 w2x_relax_job.py --kind full --d 7 --par 1 --prime 2147483629 --threads 32 --timeout 14400
