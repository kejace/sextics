cd /tmp/hodge/results/short-power-presentation
python3 w2_kernel.py 9 0 --primes 2 --maxprimes 4 --threads 2 --tag chain9 > w2_kernel_d9_p0.out 2>&1
python3 w2_kernel.py 9 1 --primes 2 --maxprimes 4 --threads 2 --tag chain9 > w2_kernel_d9_p1.out 2>&1
python3 w2_families.py 9 0 --plist 536870909,536870879 --maxr 5 --tag w2 > w2_families_d9_p0.out 2>&1
python3 w2_families.py 9 1 --plist 536870909,536870879 --maxr 5 --tag w2 > w2_families_d9_p1.out 2>&1
python3 w2_cands.py 9 0 --plist 536870909,536870879 --tag w2 > w2_cands_d9_p0.out 2>&1
python3 w2_cands.py 9 1 --plist 536870909,536870879 --tag w2 > w2_cands_d9_p1.out 2>&1
echo "chain9 done $(date)"
