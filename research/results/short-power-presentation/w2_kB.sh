cd /tmp/hodge/results/short-power-presentation
for spec in "9 0" "9 1"; do
  set -- $spec
  echo "== kernel d=$1 par=$2 $(date)"
  python3 w2_kernel.py $1 $2 --primes 3 --maxprimes 6 --threads 8 --numerators --tag kB > w2_kernel_d$1_p$2.out 2>&1
  tail -1 w2_kernel_d$1_p$2.out
done
echo "== done $(date)"
