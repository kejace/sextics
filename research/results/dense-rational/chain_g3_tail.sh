for i in 0 1 2; do timeout 7200 python3 e2norm.py e1kernel_d10_par0_reps.json $i >> e2norm_d10.log 2>&1; done; echo done > chain_g3.done
